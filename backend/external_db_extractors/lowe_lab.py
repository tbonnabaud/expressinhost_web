import asyncio
import re
import time
from dataclasses import dataclass

from aiohttp import ClientSession
from selectolax.parser import HTMLParser, Node

from ..schemas import CodonTableFormWithTranslations, CodonTranslation
from .mappings import AMINO_ACID_MAPPING, WOBBLE_MAPPING

# from ..crud.codon_tables import CodonTableRepository
# from ..crud.codon_translations import CodonTranslationRepository
# from ..database import LocalSession, Session

BASE_URL = "https://gtrnadb.ucsc.edu"
GENOME_LIST_URL = f"{BASE_URL}/cgi-bin/trna_chooseorg?org="

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
                print(f"ERROR: {amino_acid=}, {anticodon=}, {default_codon=}")

        else:
            yield CodonTranslation(
                codon=default_codon,
                anticodon=anticodon,
                amino_acid=amino_acid,
                trna_gcn=float(tnra_count.count),
                wobble_codon="---",
                wobble_rate=0.0,
            )


def extract_amino_acid_table(table_node: Node) -> list[CodonTranslation]:
    """Nodes of tbody tag inside tRNA-box class."""
    node_iter = table_node.iter()
    # Skip the first node (table header)
    next(node_iter)

    translations = []

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
            translations.extend(get_amino_acid_translations(amino_acid, trna_sublist))

    return translations


def parse_trna_gene_summary(html_content: str) -> list[CodonTranslation]:
    tree = HTMLParser(html_content)
    # organism = tree.css_first("#page-header h5").text()
    # print(organism)
    translations = []

    tables = tree.css(".tRNA-box tbody")

    for table_node in tables:
        amino_acid_translations = extract_amino_acid_table(table_node)
        translations.extend(amino_acid_translations)

    return translations


def parse_genome_list(html_content: str):
    tree = HTMLParser(html_content)
    link_list = tree.css("li a")

    for link_node in link_list:
        href = link_node.attributes.get("href")
        # Remove the first part containing ".."
        page_link = f"{BASE_URL}{href[2:]}"

        yield GenomePageMetadata(link_node.text(), page_link)


async def get_url_content(session: ClientSession, url: str):
    async with session.get(url, ssl=False) as response:
        if response.status == 200:
            return await response.text()

        else:
            print("Error:", response.url, response.status, response.reason)


async def task(session: ClientSession, genome_page: GenomePageMetadata):
    summary_page = await get_url_content(session, genome_page.link)

    if summary_page:
        translations = parse_trna_gene_summary(summary_page)
        codon_table = CodonTableFormWithTranslations(
            organism=genome_page.organism,
            name="Default",
            source="Lowe Lab",
            translations=sorted(translations, key=lambda t: t.amino_acid),
        )

        return codon_table


async def run_scraping():
    start_time = time.time()

    async with ClientSession() as session:
        genome_list_page = await get_url_content(session, GENOME_LIST_URL)

        genome_list = parse_genome_list(genome_list_page)

        results = await asyncio.gather(
            *[task(session, genome_page) for genome_page in genome_list]
        )

        print(len(results))
        print(results[0])

        # gene_summary_page = await get_url_content(
        #     session,
        #     "https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/bacteria/Esch_coli/",
        # )

        # if gene_summary_page:
        #     result = parse_trna_gene_summary(gene_summary_page)
        #     pprint(result)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")


if __name__ == "__main__":
    asyncio.run(run_scraping())
