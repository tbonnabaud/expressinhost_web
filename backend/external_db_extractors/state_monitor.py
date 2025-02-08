from dataclasses import dataclass
from enum import Enum


class Status(str, Enum):
    idle = "Idle"
    fetching_genome_list = "Fetching list of genomes"
    parsing_codon_tables = "Parsing codon tables"
    database_insertion = "Insertion in the database"
    error = "Error"
    success = "Success"


@dataclass
class StateMonitor:
    status: Status = Status.idle
    done: int | None = None
    total: int | None = None

    def reset(self):
        self.status = Status.idle
        self.done = None
        self.total = None

    def start_with_total(self, total: int):
        self.status = Status.fetching_genome_list
        self.done = 0
        self.total = total

    def progress(self):
        self.done += 1

        if self.done == self.total:
            self.status = Status.success

    def error(self):
        self.status = Status.error
