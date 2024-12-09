# Amino-acid to anticodon to codon
AMINO_ACID_MAPPING = {
    "Ala": {"AGC": "GCU", "GGC": "GCC", "CGC": "GCG", "UGC": "GCA"},
    "Arg": {
        "ACG": "CGU",
        "GCG": "CGC",
        "CCG": "CGG",
        "UCG": "CGA",
        "CCU": "AGG",
        "UCU": "AGA",
    },
    "Asn": {"AUU": "AAU", "GUU": "AAC"},
    "Asp": {
        "AUC": "GAU",
        "GUC": "GAC",
    },
    "Cys": {"ACA": "UGU", "GCA": "UGC"},
    "Gln": {"CUG": "CAG", "UUG": "CAA"},
    "Glu": {"CUC": "GAG", "UUC": "GAA"},
    "Gly": {"ACC": "GGU", "GCC": "GGC", "CCC": "GGG", "UCC": "GGA"},
    "His": {"AUG": "CAU", "GUG": "CAC"},
    "Ile": {"AAU": "AUU", "GAU": "AUC", "UAU": "AUA"},
    "Leu": {
        "AAG": "CUU",
        "GAG": "CUC",
        "CAG": "CUG",
        "UAG": "CUA",
        "CAA": "UUG",
        "UAA": "UUA",
    },
    "Lys": {"CUU": "AAG", "UUU": "AAA"},
    "Met": {"CAU": "AUG"},
    "Phe": {"AAA": "UUU", "GAA": "UUC"},
    "Pro": {"AGG": "CCU", "GGG": "CCC", "CGG": "CCG", "UGG": "CCA"},
    "Ser": {
        "AGA": "UCU",
        "GGA": "UCC",
        "CGA": "UCG",
        "UGA": "UCA",
        "ACU": "AGU",
        "GCU": "AGC",
    },
    "Thr": {"AGU": "ACU", "GGU": "ACC", "CGU": "ACG", "UGU": "ACA"},
    "Trp": {"CCA": "UGG"},
    "Tyr": {"AUA": "UAU", "GUA": "UAC"},
    "Val": {"AAC": "GUU", "GAC": "GUC", "CAC": "GUG", "UAC": "GUA"},
}

# Amino-acid to codon to list of potential wobble codons
WOBBLE_MAPPING = {
    "Ala": {
        "GCA": ["GCG", "GCC"],
        "GCG": ["GCA", "GCC"],
        "GCC": ["GCU", "GCG"],
        "GCU": ["GCC", "GCG"],
    },
    "Arg": {
        "CGA": ["CGG", "CGC"],
        "CGG": ["CGA", "CGC"],
        "CGU": ["CGC", "CGG"],
        "CGC": ["CGU", "CGG"],
        "AGA": ["AGG", "CGG"],
        "AGG": ["AGA", "CGG"],
    },
    "Asn": {"AAU": ["AAC", "AAC"], "AAC": ["AAU", "AAU"]},
    "Asp": {"GAU": ["GAC", "GAC"], "GAC": ["GAU", "GAU"]},
    "Cys": {"UGU": ["UGC", "UGC"], "UGC": ["UGU", "UGU"]},
    "Gln": {"CAA": ["CAG", "CAG"], "CAG": ["CAA", "CAA"]},
    "Glu": {"GAA": ["GAG", "GAG"], "GAG": ["GAA", "GAA"]},
    "Gly": {
        "GGA": ["GGG", "GGC"],
        "GGG": ["GGA", "GGC"],
        "GGU": ["GGC", "GGG"],
        "GGC": ["GGU", "GGG"],
    },
    "His": {"CAU": ["CAC", "CAC"], "CAC": ["CAU", "CAU"]},
    "Ile": {"AUA": ["AUC", "AUU"], "AUU": ["AUC", "AUA"], "AUC": ["AUU", "AUA"]},
    "Leu": {
        "CUA": ["CUG", "CUC"],
        "CUG": ["CUA", "CUC"],
        "CUU": ["CUC", "CUG"],
        "CUC": ["CUU", "CUG"],
        "UUA": ["UUG", "CUG"],
        "UUG": ["UUA", "CUG"],
    },
    "Lys": {"AAA": ["AAG", "AAG"], "AAG": ["AAA", "AAA"]},
    "Met": {"AUG": ["AUG", "AUG"]},
    "Phe": {"UUU": ["UUC", "UUC"], "UUC": ["UUU", "UUU"]},
    "Pro": {
        "CCA": ["CCG", "CCC"],
        "CCG": ["CCA", "CCC"],
        "CCU": ["CCC", "CCG"],
        "CCC": ["CCU", "CCG"],
    },
    "Ser": {
        "UCA": ["UCG", "UCC"],
        "UCG": ["UCA", "UCC"],
        "UCU": ["UCC", "UCG"],
        "UCC": ["UCU", "UCG"],
        "AGU": ["AGC", "UCC"],
        "AGC": ["AGU", "UCC"],
    },
    "Thr": {
        "ACA": ["ACG", "ACC"],
        "ACG": ["ACA", "ACC"],
        "ACU": ["ACC", "ACG"],
        "ACC": ["ACU", "ACG"],
    },
    "Trp": {"UGG": ["UGG", "UGG"]},
    "Tyr": {"UAU": ["UAC", "UAC"], "UAC": ["UAU", "UAU"]},
    "Val": {
        "GUA": ["GUG", "GUC"],
        "GUG": ["GUA", "GUC"],
        "GUU": ["GUC", "GUG"],
        "GUC": ["GUU", "GUG"],
    },
    "Stop": {"UAA": ["UAG", "UGA"], "UAG": ["UAA", "UGA"], "UGA": ["UAA", "UAG"]},
}

if __name__ == "__main__":
    import csv
    from collections import defaultdict

    with open("backend/external_db_extractors/wobble_table.csv") as f:
        reader = csv.DictReader(f, delimiter="\t")
        output = defaultdict(dict)

        for row in reader:
            output[row["amino_acid"]][row["codon"]] = [row["first"], row["second"]]

        print(output)
