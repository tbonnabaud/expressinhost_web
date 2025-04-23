import asyncio
import re
import time
from dataclasses import dataclass

from aiohttp import ClientSession
from rq import get_current_job
from selectolax.parser import HTMLParser, Node

from ..crud.codon_tables import CodonTableRepository
from ..crud.codon_translations import CodonTranslationRepository
from ..crud.last_web_scraping import LastWebScrapingRepository
from ..database import IntegrityError, LocalSession, SQLAlchemyError
from ..job_manager import Job, update_job_meta, web_scraping_queue
from ..logger import scraping_logger
from ..routes.codon_tables import assign_codon_table_id
from ..schemas import CodonTableFormWithTranslations, CodonTranslation
from .mappings import AMINO_ACID_MAPPING, WOBBLE_MAPPING

SOURCE = "Lowe Lab"
BASE_URL = "https://gtrnadb.ucsc.edu"
GENOME_LIST_URL = f"{BASE_URL}/browse.html"

# Third group to handle cell containing X/Y
CELL_REGEX = re.compile(r"([A-Z]+)\s+(\d*)/?(\d*)")


@dataclass(slots=True)
class TRNACell:
    anticodon: str
    count: int


@dataclass(slots=True)
class GenomePageMetadata:
    organism: str
    link: str


def get_codon_from_anticodon(amino_acid: str, anticodon: str) -> str | None:
    anticodon_mapping = AMINO_ACID_MAPPING.get(amino_acid)

    if anticodon_mapping:
        return anticodon_mapping.get(anticodon)


def get_potential_wobble_codons(amino_acid: str, codon: str) -> list[str]:
    return WOBBLE_MAPPING[amino_acid].get(codon)


def parse_cell(td_node: Node) -> TRNACell | None:
    """Parse the content of the cell to extract the anticodon and its count."""
    cell_text = td_node.text()
    m = CELL_REGEX.match(cell_text)

    if m:
        trna = m[2] if m[3] == "" else m[3]
        anticodon = m[1].replace("T", "U")

        if trna == "":
            return TRNACell(anticodon, 0)

        return TRNACell(anticodon, int(trna))

    return None


def potential_wobble_codon_has_trna(
    amino_acid: str, potential_wobble_codon: str, trna_sublist: list[TRNACell]
) -> bool:
    for t in trna_sublist:
        codon = AMINO_ACID_MAPPING[amino_acid][t.anticodon]

        if potential_wobble_codon == codon and t.count > 0:
            return True

    return False


def get_amino_acid_translations(amino_acid: str, trna_sublist: list[TRNACell]):
    for tnra_count in trna_sublist:
        anticodon = tnra_count.anticodon
        default_codon = get_codon_from_anticodon(amino_acid, anticodon)

        if tnra_count.count == 0:
            potential_wobble_codons = get_potential_wobble_codons(
                amino_acid, default_codon
            )

            if potential_wobble_codons:
                wobble_codon = "---"

                if potential_wobble_codon_has_trna(
                    amino_acid, potential_wobble_codons[0], trna_sublist
                ):
                    wobble_codon = potential_wobble_codons[0]
                elif potential_wobble_codon_has_trna(
                    amino_acid, potential_wobble_codons[1], trna_sublist
                ):
                    wobble_codon = potential_wobble_codons[1]

                yield CodonTranslation(
                    codon=default_codon,
                    anticodon=anticodon,
                    amino_acid=amino_acid,
                    trna_gcn=0.0,
                    wobble_codon=wobble_codon,
                    wobble_rate=0.35,
                )

            else:
                scraping_logger.error(
                    f"No potential wobble codons for: {amino_acid=}, {anticodon=}, {default_codon=}"
                )

        else:
            yield CodonTranslation(
                codon=default_codon,
                anticodon=anticodon,
                amino_acid=amino_acid,
                trna_gcn=float(tnra_count.count),
                wobble_codon="---",
                wobble_rate=0.0,
            )


def extract_amino_acid_table(table_node: Node):
    """Nodes of tbody tag inside tRNA-box class."""
    node_iter = table_node.iter()
    # Skip the first node (table header)
    next(node_iter)

    for row_node in node_iter:
        row_iter = row_node.iter()
        amino_acid = next(row_iter).text()

        # Ignore stop codon
        if amino_acid not in ["Supres", "SelCys"]:
            if amino_acid in ["fMet/Met", "iMet/Met"]:
                amino_acid = "Met"

            trna_sublist = [
                count
                for count in map(parse_cell, row_iter)
                # We ignore CAU anticodon for isoleucine
                if count and not (amino_acid == "Ile" and count.anticodon == "CAU")
            ]

            for translation in get_amino_acid_translations(amino_acid, trna_sublist):
                yield translation


def parse_trna_gene_summary(html_content: str):
    tree = HTMLParser(html_content)
    # organism = tree.css_first("#page-header h5").text()
    # scraping_logger.debug(organism)
    tables = tree.css(".tRNA-box tbody")

    for table_node in tables:
        for translation in extract_amino_acid_table(table_node):
            yield translation


def parse_genome_list(html_content: str):
    tree = HTMLParser(html_content)
    link_list = tree.css(".panel td a")

    for link_node in link_list:
        href = link_node.attributes.get("href")
        # Remove the first part containing ".."
        page_link = f"{BASE_URL}/{href}"

        yield GenomePageMetadata(link_node.text(), page_link)


async def get_url_content(session: ClientSession, url: str):
    async with session.get(url, ssl=False) as response:
        if response.status == 200:
            return await response.text()

        else:
            scraping_logger.error(
                f"{response.status} {response.reason}: {response.url}"
            )


async def task(job: Job, session: ClientSession, genome_page: GenomePageMetadata):
    summary_page = await get_url_content(session, genome_page.link)
    codon_table = None
    step = job.get_meta()["step"]
    update_job_meta(job, genome_page.organism, step + 1)

    if summary_page:
        translations = parse_trna_gene_summary(summary_page)
        codon_table = CodonTableFormWithTranslations(
            organism=genome_page.organism,
            name="Default",
            source=SOURCE,
            translations=sorted(translations, key=lambda t: t.amino_acid),
        )

    return codon_table


async def get_last_release():
    async with ClientSession() as session:
        home_page = await get_url_content(session, BASE_URL)
        tree = HTMLParser(home_page)
        release = tree.css_first("#homepage-section h5").text()

        return release.replace("Data ", "")


async def run_scraping():
    job = get_current_job()
    start_time = time.time()
    scraping_logger.info("Fetching and parsing HTML data...")
    update_job_meta(job, "Starting...", 0, 0)
    last_release = await get_last_release()

    async with ClientSession() as session:
        genome_list_page = await get_url_content(session, GENOME_LIST_URL)
        genome_list = list(parse_genome_list(genome_list_page))

        update_job_meta(job, "Fetching list of genomes...", 0, 2 * len(genome_list))

        results = await asyncio.gather(
            *[task(job, session, genome_page) for genome_page in genome_list]
        )

    insert_message = "Insertion of the extracted codon tables in the database..."
    scraping_logger.info(insert_message)

    with LocalSession() as db_session:
        codon_table_repo = CodonTableRepository(db_session)
        codon_translation_repo = CodonTranslationRepository(db_session)
        step = job.get_meta()["step"] + 1

        # To avoid to try to insert duplicates from the page listing available genomes
        inserted_organisms = set()

        for result in filter(lambda r: r is not None, results):
            update_job_meta(job, insert_message, step)

            if result.organism not in inserted_organisms:
                try:
                    meta_dict = result.model_dump(exclude={"translations"})
                    meta_dict["user_id"] = None
                    codon_table_id = codon_table_repo.add(meta_dict)
                    codon_translation_repo.add_batch(
                        [
                            assign_codon_table_id(codon_table_id, x)
                            for x in result.translations
                        ]
                    )

                    db_session.commit()
                    inserted_organisms.add(result.organism)

                except IntegrityError:
                    scraping_logger.warning(
                        f"{result.organism} - {result.name} already exists."
                    )
                    db_session.rollback()

                except SQLAlchemyError as exc:
                    scraping_logger.error(exc)
                    db_session.rollback()

                step += 1

        LastWebScrapingRepository(db_session).upsert(SOURCE, last_release)

    end_time = time.time()
    elapsed_time = end_time - start_time
    end_message = f"Elapsed time for the web scraping: {elapsed_time:.0f} seconds. {len(inserted_organisms)} new table(s) inserted."
    scraping_logger.info(end_message)


async def periodic_web_scraping():
    """Example background task that runs every 24 hours."""
    while True:
        try:
            scraping_logger.info("Check if database is up to date.")
            last_release = await get_last_release()
            scraping_in_db = None
            scraping_logger.info(f"Last release is {last_release}")

            with LocalSession() as db_session:
                scraping_in_db = LastWebScrapingRepository(db_session).get_last_from(
                    SOURCE
                )

            if not scraping_in_db or scraping_in_db.release != last_release:
                scraping_logger.info("Database is not up to date. Run web scraping...")
                web_scraping_queue.enqueue(run_scraping)

            else:
                scraping_logger.info("Database is up to date.")

            await asyncio.sleep(60 * 60 * 24)

        except asyncio.CancelledError:
            scraping_logger.info("Periodic web scraping was cancelled.")
            break

        except Exception as e:
            scraping_logger.error(f"Error in periodic web scraping: {str(e)}")
            await asyncio.sleep(60 * 60)  # Wait before retrying


if __name__ == "__main__":
    asyncio.run(periodic_web_scraping())
