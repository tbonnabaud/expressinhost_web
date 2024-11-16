import csv
from pathlib import Path

from .core.constantes import TABLE_BASE_PATH
from .crud.codon_tables import CodonTableRepository
from .crud.codon_translations import CodonTranslationRepository
from .database import LocalSession
from .schemas import CodonTableForm, CodonTranslation

# import re


# def format_to_binomial(value: str):
#     """Return value with binomial nomenclature."""
#     return re.sub(r"[\s_\-]+", " ", value).capitalize()


def get_csv_file_list(dir_path: Path):
    for path in dir_path.iterdir():
        if path.is_file() and path.name.endswith(".csv"):
            yield path


def create_codon_table(organism: str, name: str):
    data = CodonTableForm(organism=organism, name=name).model_dump()
    data["user_id"] = None

    return data


def main():
    with LocalSession() as session:
        codon_table_repo = CodonTableRepository(session)

        if codon_table_repo.count_all() == 0:
            codon_translation_repo = CodonTranslationRepository(session)

            for path in get_csv_file_list(TABLE_BASE_PATH):
                with path.open() as file:
                    organism_name = path.name.removesuffix(".csv")
                    codon_table = create_codon_table(
                        organism_name, organism_name + "_example"
                    )
                    codon_table_id = codon_table_repo.add(codon_table)

                    reader = csv.DictReader(file, delimiter="\t")

                    row_list = [
                        CodonTranslation(
                            codon_table_id=codon_table_id, **row
                        ).model_dump()
                        for row in reader
                    ]

                    codon_translation_repo.add_batch(row_list)


if __name__ == "__main__":
    main()
