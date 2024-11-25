import csv
import re
from pathlib import Path

from .core.constantes import TABLE_BASE_PATH
from .crud.codon_tables import CodonTableRepository
from .crud.codon_translations import CodonTranslationRepository
from .database import LocalSession
from .routes.codon_tables import assign_codon_table_id
from .schemas import CodonTableFormWithTranslations, CodonTranslation


def get_csv_file_list(dir_path: Path):
    for path in dir_path.iterdir():
        if path.is_file() and path.name.endswith(".csv"):
            yield path


def read_rows_from_csv(path: Path) -> list[dict]:
    with path.open() as file:
        reader = csv.DictReader(file, delimiter="\t")

        return list(reader)


def main():
    with LocalSession() as session:
        codon_table_repo = CodonTableRepository(session)

        if codon_table_repo.count_all() == 0:
            codon_translation_repo = CodonTranslationRepository(session)

            for path in get_csv_file_list(TABLE_BASE_PATH):
                csv_rows = read_rows_from_csv(path)
                translations = [CodonTranslation(**row) for row in csv_rows]

                # Use binomial convention
                organism_name = re.sub(
                    r"[\s_\-]+", " ", path.name.removesuffix(".csv")
                ).capitalize()
                codon_table_form = CodonTableFormWithTranslations(
                    organism=organism_name, name="Example", translations=translations
                )
                meta_dict = codon_table_form.model_dump(exclude={"translations"})
                meta_dict["user_id"] = None
                codon_table_id = codon_table_repo.add(meta_dict)

                codon_translation_repo.add_batch(
                    list(
                        map(
                            lambda x: assign_codon_table_id(codon_table_id, x),
                            translations,
                        )
                    )
                )


if __name__ == "__main__":
    main()
