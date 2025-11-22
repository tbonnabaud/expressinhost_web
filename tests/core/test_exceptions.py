from backend.core.exceptions import (
    ClustalFormatError,
    ExpressInHostError,
    FastaFormatError,
    NoAminoAcidConservation,
    NoIdenticalSequencesError,
)


def test_exception_hierarchy():
    """Test the exception hierarchy structure."""
    # All custom exceptions should be catchable as ExpressInHostError
    try:
        raise FastaFormatError("test")
    except ExpressInHostError:
        pass

    try:
        raise ClustalFormatError("test")
    except ExpressInHostError:
        pass

    try:
        raise NoIdenticalSequencesError("test")
    except ExpressInHostError:
        pass

    try:
        raise NoAminoAcidConservation("test")
    except ExpressInHostError:
        pass
