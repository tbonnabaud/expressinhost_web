import pytest

from backend.core.exceptions import (
    ClustalFormatError,
    ExpressInHostError,
    FastaFormatError,
    NoAminoAcidConservation,
    NoIdenticalSequencesError,
)


def test_express_in_host_error_is_exception():
    """Test that ExpressInHostError is a subclass of Exception."""
    assert issubclass(ExpressInHostError, Exception)


def test_express_in_host_error_can_be_raised():
    """Test that ExpressInHostError can be raised with a message."""
    with pytest.raises(ExpressInHostError) as exc_info:
        raise ExpressInHostError("Test error message")
    assert str(exc_info.value) == "Test error message"


def test_fasta_format_error_is_express_in_host_error():
    """Test that FastaFormatError is a subclass of ExpressInHostError."""
    assert issubclass(FastaFormatError, ExpressInHostError)


def test_fasta_format_error_can_be_raised():
    """Test that FastaFormatError can be raised with a message."""
    with pytest.raises(FastaFormatError) as exc_info:
        raise FastaFormatError("Invalid FASTA format")
    assert str(exc_info.value) == "Invalid FASTA format"


def test_clustal_format_error_is_express_in_host_error():
    """Test that ClustalFormatError is a subclass of ExpressInHostError."""
    assert issubclass(ClustalFormatError, ExpressInHostError)


def test_clustal_format_error_can_be_raised():
    """Test that ClustalFormatError can be raised with a message."""
    with pytest.raises(ClustalFormatError) as exc_info:
        raise ClustalFormatError("Invalid CLUSTAL format")
    assert str(exc_info.value) == "Invalid CLUSTAL format"


def test_no_identical_sequences_error_is_express_in_host_error():
    """Test that NoIdenticalSequencesError is a subclass of ExpressInHostError."""
    assert issubclass(NoIdenticalSequencesError, ExpressInHostError)


def test_no_identical_sequences_error_can_be_raised():
    """Test that NoIdenticalSequencesError can be raised with a message."""
    with pytest.raises(NoIdenticalSequencesError) as exc_info:
        raise NoIdenticalSequencesError("Sequences do not match")
    assert str(exc_info.value) == "Sequences do not match"


def test_no_amino_acid_conservation_is_express_in_host_error():
    """Test that NoAminoAcidConservation is a subclass of ExpressInHostError."""
    assert issubclass(NoAminoAcidConservation, ExpressInHostError)


def test_no_amino_acid_conservation_can_be_raised():
    """Test that NoAminoAcidConservation can be raised with a message."""
    with pytest.raises(NoAminoAcidConservation) as exc_info:
        raise NoAminoAcidConservation("Amino acids not conserved")
    assert str(exc_info.value) == "Amino acids not conserved"


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
