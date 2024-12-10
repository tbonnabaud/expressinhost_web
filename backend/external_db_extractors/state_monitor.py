from dataclasses import dataclass
from enum import Enum


class State(str, Enum):
    idle = "Idle"
    fetching_genome_list = "Fetching list of genomes"
    parsing_codon_tables = "Parsing codon tables"
    database_insertion = "Insertion in the database"
    error = "Error"
    success = "Success"


@dataclass
class StateMonitor:
    state: State = State.idle
    done: int | None = None
    total: int | None = None

    def reset(self):
        self.state = State.idle
        self.done = None
        self.total = None

    def start_with_total(self, total: int):
        self.state = State.fetching_genome_list
        self.done = 0
        self.total = total

    def progress(self):
        self.done += 1

        if self.done == self.total:
            self.state = State.success

    def error(self):
        self.state = State.error
