from enum import Enum

from pydantic import BaseModel


class Status(str, Enum):
    IDLE = "Idle"
    RUNNING = "Running"
    ERROR = "Error"
    SUCCESS = "Success"


class StateMonitor(BaseModel):
    status: Status = Status.IDLE
    message: str = ""
    done: int | None = None
    total: int | None = None

    def reset(self):
        self.status = Status.IDLE
        self.done = None
        self.total = None
        self.message = ""

    def start(self, message: str = ""):
        self.status = Status.RUNNING
        self.done = None
        self.total = None
        self.message = message

    def set_total(self, total: int):
        self.status = Status.RUNNING
        self.done = 0
        self.total = total

    def progress(self):
        self.done += 1

        if self.done == self.total:
            self.status = Status.SUCCESS

    def error(self):
        self.status = Status.ERROR
