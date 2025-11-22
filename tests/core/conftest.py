import csv
from io import StringIO

import pytest


@pytest.fixture
def raw_codon_table_data():
    """Fixture providing raw codon table data for testing."""
    csv_data = """amino_acid,anticodon,codon,trna_gcn,wobble_codon,wobble_rate
Ala,IGC,GCU,2.5,GCU,0.0
Ala,GGC,GCC,3.0,GCC,0.0
Ala,CGC,GCG,1.5,GCG,0.0
Ala,UGC,GCA,2.0,GCA,0.0
Arg,ICG,CGU,1.0,CGU,0.0
Arg,GCG,CGC,2.0,CGC,0.0
Arg,CCG,CGG,0.5,CGG,0.0
Arg,UCG,CGA,0.8,CGA,0.0
Arg,ICU,AGA,0.0,CGU,0.3
Arg,GCU,AGG,0.0,CGC,0.3
Asn,IGU,AAU,3.0,AAU,0.0
Asn,GUU,AAC,4.0,AAC,0.0
Asp,IGU,GAU,2.5,GAU,0.0
Asp,GUC,GAC,3.5,GAC,0.0
Cys,IGC,UGU,1.5,UGU,0.0
Cys,GCA,UGC,2.5,UGC,0.0
Gln,IUG,CAG,3.0,CAG,0.0
Gln,CUG,CAA,2.0,CAA,0.0
Glu,IUC,GAG,3.5,GAG,0.0
Glu,CUC,GAA,2.5,GAA,0.0
Gly,ICC,GGU,2.0,GGU,0.0
Gly,GCC,GGC,3.0,GGC,0.0
Gly,CCC,GGG,1.5,GGG,0.0
Gly,UCC,GGA,2.5,GGA,0.0
His,IGU,CAU,2.0,CAU,0.0
His,GUG,CAC,3.0,CAC,0.0
Ile,IAU,AUU,3.0,AUU,0.0
Ile,GAU,AUC,4.0,AUC,0.0
Ile,CAU,AUA,1.5,AUA,0.0
Leu,IAG,CUG,4.0,CUG,0.0
Leu,GAG,CUC,3.0,CUC,0.0
Leu,CAG,CUU,2.0,CUU,0.0
Leu,UAG,CUA,1.0,CUA,0.0
Leu,IAA,UUG,2.5,UUG,0.0
Leu,CAA,UUA,1.5,UUA,0.0
Lys,IUU,AAG,3.5,AAG,0.0
Lys,CUU,AAA,4.5,AAA,0.0
Met,IAU,AUG,5.0,AUG,0.0
Phe,IAA,UUU,3.0,UUU,0.0
Phe,GAA,UUC,4.0,UUC,0.0
Pro,IGG,CCU,2.0,CCU,0.0
Pro,GGG,CCC,3.0,CCC,0.0
Pro,CGG,CCG,1.5,CCG,0.0
Pro,UGG,CCA,2.5,CCA,0.0
Ser,IGA,UCU,2.5,UCU,0.0
Ser,GGA,UCC,3.5,UCC,0.0
Ser,CGA,UCG,1.5,UCG,0.0
Ser,UGA,UCA,2.0,UCA,0.0
Ser,ICU,AGU,0.0,UCU,0.3
Ser,GCU,AGC,0.0,UCC,0.3
Thr,IGU,ACU,2.5,ACU,0.0
Thr,GGU,ACC,3.5,ACC,0.0
Thr,CGU,ACG,1.5,ACG,0.0
Thr,UGU,ACA,2.0,ACA,0.0
Trp,ICA,UGG,4.0,UGG,0.0
Tyr,IUA,UAU,2.5,UAU,0.0
Tyr,GUA,UAC,3.5,UAC,0.0
Val,IAC,GUU,3.0,GUU,0.0
Val,GAC,GUC,4.0,GUC,0.0
Val,CAC,GUG,2.0,GUG,0.0
Val,UAC,GUA,2.5,GUA,0.0"""

    reader = csv.DictReader(StringIO(csv_data))
    return list(reader)
