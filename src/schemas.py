from pydantic import BaseModel


class BadSequence(BaseModel):
    name: str
    msg: str


class MismatchingSequences(BaseModel):
    name: str
    sequence1: str
    sequence2: str


class TuningParameters(BaseModel):
    slow_speed_threshold: float
    conservation_threshold: float
