class ExpressInHostError(Exception):
    pass


class FastaFormatError(ExpressInHostError):
    pass


class ClustalFormatError(ExpressInHostError):
    pass


class NoIdenticalSequencesError(ExpressInHostError):
    pass


class NoAminoAcidConservation(ExpressInHostError):
    pass
