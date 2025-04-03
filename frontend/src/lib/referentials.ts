import type { CodonTranslation, RestrictionSite } from './interfaces'
import { groupByAminoAcid } from './helpers'

export const MODE_LABEL_MAPPING: Record<string, string> = {
  direct_mapping: 'Direct mapping',
  optimisation_and_conservation_1: 'Optimisation and conservation 1',
  optimisation_and_conservation_2: 'Optimisation and conservation 2',
}

export const FIVE_PRIME_REGION_TUNING_LABEL_MAPPING: Record<string, string> = {
  partial_untuning: 'Partial untuning',
  fine_tuning: 'Fine tuning',
}

const BASE_CODON_TABLE = [
  { amino_acid: 'Ala', anticodon: 'AGC', codon: 'GCU' },
  { amino_acid: 'Ala', anticodon: 'GGC', codon: 'GCC' },
  { amino_acid: 'Ala', anticodon: 'CGC', codon: 'GCG' },
  { amino_acid: 'Ala', anticodon: 'UGC', codon: 'GCA' },
  { amino_acid: 'Arg', anticodon: 'ACG', codon: 'CGU' },
  { amino_acid: 'Arg', anticodon: 'GCG', codon: 'CGC' },
  { amino_acid: 'Arg', anticodon: 'CCG', codon: 'CGG' },
  { amino_acid: 'Arg', anticodon: 'UCG', codon: 'CGA' },
  { amino_acid: 'Arg', anticodon: 'CCU', codon: 'AGG' },
  { amino_acid: 'Arg', anticodon: 'UCU', codon: 'AGA' },
  { amino_acid: 'Asn', anticodon: 'AUU', codon: 'AAU' },
  { amino_acid: 'Asn', anticodon: 'GUU', codon: 'AAC' },
  { amino_acid: 'Asp', anticodon: 'AUC', codon: 'GAU' },
  { amino_acid: 'Asp', anticodon: 'GUC', codon: 'GAC' },
  { amino_acid: 'Cys', anticodon: 'ACA', codon: 'UGU' },
  { amino_acid: 'Cys', anticodon: 'GCA', codon: 'UGC' },
  { amino_acid: 'Gln', anticodon: 'CUG', codon: 'CAG' },
  { amino_acid: 'Gln', anticodon: 'UUG', codon: 'CAA' },
  { amino_acid: 'Glu', anticodon: 'CUC', codon: 'GAG' },
  { amino_acid: 'Glu', anticodon: 'UUC', codon: 'GAA' },
  { amino_acid: 'Gly', anticodon: 'ACC', codon: 'GGU' },
  { amino_acid: 'Gly', anticodon: 'GCC', codon: 'GGC' },
  { amino_acid: 'Gly', anticodon: 'CCC', codon: 'GGG' },
  { amino_acid: 'Gly', anticodon: 'UCC', codon: 'GGA' },
  { amino_acid: 'His', anticodon: 'AUG', codon: 'CAU' },
  { amino_acid: 'His', anticodon: 'GUG', codon: 'CAC' },
  { amino_acid: 'Ile', anticodon: 'AAU', codon: 'AUU' },
  { amino_acid: 'Ile', anticodon: 'GAU', codon: 'AUC' },
  { amino_acid: 'Ile', anticodon: 'UAU', codon: 'AUA' },
  { amino_acid: 'Leu', anticodon: 'AAG', codon: 'CUU' },
  { amino_acid: 'Leu', anticodon: 'GAG', codon: 'CUC' },
  { amino_acid: 'Leu', anticodon: 'CAG', codon: 'CUG' },
  { amino_acid: 'Leu', anticodon: 'UAG', codon: 'CUA' },
  { amino_acid: 'Leu', anticodon: 'CAA', codon: 'UUG' },
  { amino_acid: 'Leu', anticodon: 'UAA', codon: 'UUA' },
  { amino_acid: 'Lys', anticodon: 'CUU', codon: 'AAG' },
  { amino_acid: 'Lys', anticodon: 'UUU', codon: 'AAA' },
  { amino_acid: 'Met', anticodon: 'CAU', codon: 'AUG' },
  { amino_acid: 'Phe', anticodon: 'AAA', codon: 'UUU' },
  { amino_acid: 'Phe', anticodon: 'GAA', codon: 'UUC' },
  { amino_acid: 'Pro', anticodon: 'AGG', codon: 'CCU' },
  { amino_acid: 'Pro', anticodon: 'GGG', codon: 'CCC' },
  { amino_acid: 'Pro', anticodon: 'CGG', codon: 'CCG' },
  { amino_acid: 'Pro', anticodon: 'UGG', codon: 'CCA' },
  { amino_acid: 'Ser', anticodon: 'AGA', codon: 'UCU' },
  { amino_acid: 'Ser', anticodon: 'GGA', codon: 'UCC' },
  { amino_acid: 'Ser', anticodon: 'CGA', codon: 'UCG' },
  { amino_acid: 'Ser', anticodon: 'UGA', codon: 'UCA' },
  { amino_acid: 'Ser', anticodon: 'ACU', codon: 'AGU' },
  { amino_acid: 'Ser', anticodon: 'GCU', codon: 'AGC' },
  { amino_acid: 'Thr', anticodon: 'AGU', codon: 'ACU' },
  { amino_acid: 'Thr', anticodon: 'GGU', codon: 'ACC' },
  { amino_acid: 'Thr', anticodon: 'CGU', codon: 'ACG' },
  { amino_acid: 'Thr', anticodon: 'UGU', codon: 'ACA' },
  { amino_acid: 'Trp', anticodon: 'CCA', codon: 'UGG' },
  { amino_acid: 'Tyr', anticodon: 'AUA', codon: 'UAU' },
  { amino_acid: 'Tyr', anticodon: 'GUA', codon: 'UAC' },
  { amino_acid: 'Val', anticodon: 'AAC', codon: 'GUU' },
  { amino_acid: 'Val', anticodon: 'GAC', codon: 'GUC' },
  { amino_acid: 'Val', anticodon: 'CAC', codon: 'GUG' },
  { amino_acid: 'Val', anticodon: 'UAC', codon: 'GUA' },
]

export const DEFAULT_CODON_TABLE: Array<CodonTranslation> =
  BASE_CODON_TABLE.map(e => ({
    ...e,
    trna_gcn: 1,
    wobble_codon: '---',
    wobble_rate: 0,
  }))

export const CODON_LIST = BASE_CODON_TABLE.map(e => e.codon).sort()

export const CODON_MAPPING = Object.fromEntries(
  BASE_CODON_TABLE.map(e => [e.codon, e]),
)

export const AMINO_ACID_MAPPING = groupByAminoAcid(DEFAULT_CODON_TABLE)

export const UTR_EXAMPLE =
  'ACCCGGCGCTCCATTAAATAGCCGTAGACGGAACTTCGCCTTTCTCTCGGCCTTAGCGCCATTTTTTTGGGTGAGTGTTTTTTGGTTCCTGCGTTGGGATTCCGTGTACAATCCATAGACATCTGACCTCGGCACTTAGCATCATCACAGCAAACTAACTGTAGCCTTTCTCTCTTTCCCTGTAGAAACCTCTGCGCC'

// Restriction enzyme recognition sites
export const RESTRICTION_SITES: RestrictionSite[] = [
  {
    enzyme: 'AatI',
    sequence: 'AGGCCT',
  },
  {
    enzyme: 'AatII',
    sequence: 'GACGTC',
  },
  {
    enzyme: 'AbrI',
    sequence: 'CTCGAG',
  },
  {
    enzyme: 'AccII',
    sequence: 'CGCG',
  },
  {
    enzyme: 'AccIII',
    sequence: 'TCCGGA',
  },
  {
    enzyme: 'Acc65I',
    sequence: 'GGTACC',
  },
  {
    enzyme: 'AciI',
    sequence: 'CCGC',
  },
  {
    enzyme: 'AclI',
    sequence: 'AACGTT',
  },
  {
    enzyme: 'Afa22MI',
    sequence: 'CGATCG',
  },
  {
    enzyme: 'AflII',
    sequence: 'CTTAAG',
  },
  {
    enzyme: 'AgeI',
    sequence: 'ACCGGT',
  },
  {
    enzyme: 'AluI',
    sequence: 'AGCT',
  },
  {
    enzyme: 'Alw26I',
    sequence: 'GTCTC',
  },
  {
    enzyme: 'Aor13HI',
    sequence: 'TCCGGA',
  },
  {
    enzyme: 'ApaI',
    sequence: 'GGGCCC',
  },
  {
    enzyme: 'AplI',
    sequence: 'CTGCAG',
  },
  {
    enzyme: 'AscI',
    sequence: 'GGCGCGCC',
  },
  {
    enzyme: 'AseI',
    sequence: 'ATTAAT',
  },
  {
    enzyme: 'AvaIII',
    sequence: 'ATGCAT',
  },
  {
    enzyme: 'BalI',
    sequence: 'TGGCCA',
  },
  {
    enzyme: 'BamFI',
    sequence: 'GGATCC',
  },
  {
    enzyme: 'BamHI',
    sequence: 'GGATCC',
  },
  {
    enzyme: 'Bbr02I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R1.BbrUI',
    sequence: 'GGCGCC',
  },
  {
    enzyme: 'BbrUII',
    sequence: 'GTCGAC',
  },
  {
    enzyme: 'BbrUIII',
    sequence: 'CTGCAG',
  },
  {
    enzyme: 'BbvI',
    sequence: 'GCAGC',
  },
  {
    enzyme: 'BbvCI (BbvCIA)',
    sequence: 'CCTCAGC',
  },
  {
    enzyme: 'BbvCI (BbvCIB)',
    sequence: 'CCTCAGC',
  },
  {
    enzyme: 'BceSIII',
    sequence: 'ACGGC',
  },
  {
    enzyme: 'BceSIV',
    sequence: 'GCAGC',
  },
  {
    enzyme: 'BclI',
    sequence: 'TGATCA',
  },
  {
    enzyme: 'BcoDI',
    sequence: 'GTCTC',
  },
  {
    enzyme: 'BfaI (BfaIA)',
    sequence: 'CTAG',
  },
  {
    enzyme: 'BfaI (BfaIB)',
    sequence: 'CTAG',
  },
  {
    enzyme: 'BfiI',
    sequence: 'ACTGGG',
  },
  {
    enzyme: 'BfuAI',
    sequence: 'ACCTGC',
  },
  {
    enzyme: 'BfuAII',
    sequence: 'GCATGC',
  },
  {
    enzyme: 'BfuCI',
    sequence: 'GATC',
  },
  {
    enzyme: 'BglII',
    sequence: 'AGATCT',
  },
  {
    enzyme: 'BhaI',
    sequence: 'GCATC',
  },
  {
    enzyme: 'BhaII',
    sequence: 'GGCC',
  },
  {
    enzyme: 'BmrI',
    sequence: 'ACTGGG',
  },
  {
    enzyme: 'BmtI',
    sequence: 'GCTAGC',
  },
  {
    enzyme: 'BpeH640I',
    sequence: 'AAGCTT',
  },
  {
    enzyme: 'BpuJI',
    sequence: 'CCCGT',
  },
  {
    enzyme: 'BsaI',
    sequence: 'GGTCTC',
  },
  {
    enzyme: 'BscXI',
    sequence: 'GCAGGC',
  },
  {
    enzyme: 'BscXII',
    sequence: 'GATC',
  },
  {
    enzyme: 'BseYI (BseYIA)',
    sequence: 'CCCAGC',
  },
  {
    enzyme: 'BseYI (BseYIB)',
    sequence: 'CCCAGC',
  },
  {
    enzyme: 'BsmI',
    sequence: 'GAATGC',
  },
  {
    enzyme: 'BsmAI',
    sequence: 'GTCTC',
  },
  {
    enzyme: 'BsmBI',
    sequence: 'CGTCTC',
  },
  {
    enzyme: 'Bsp98I',
    sequence: 'GGATCC',
  },
  {
    enzyme: 'BspD6I (BspD6IA)',
    sequence: 'GAGTC',
  },
  {
    enzyme: 'BspD6I (BspD6IB)',
    sequence: 'GAGTC',
  },
  {
    enzyme: 'BspEI',
    sequence: 'TCCGGA',
  },
  {
    enzyme: 'BspHI',
    sequence: 'TCATGA',
  },
  {
    enzyme: 'BspNCI',
    sequence: 'CCAGA',
  },
  {
    enzyme: 'BspRI',
    sequence: 'GGCC',
  },
  {
    enzyme: 'BsrI',
    sequence: 'ACTGG',
  },
  {
    enzyme: 'BsrDI (BsrDIA)',
    sequence: 'GCAATG',
  },
  {
    enzyme: 'BsrDI (BsrDIB)',
    sequence: 'GCAATG',
  },
  {
    enzyme: 'BstFI',
    sequence: 'AAGCTT',
  },
  {
    enzyme: 'BstVI',
    sequence: 'CTCGAG',
  },
  {
    enzyme: 'BstXII',
    sequence: 'GATC',
  },
  {
    enzyme: 'BstZ1II',
    sequence: 'AAGCTT',
  },
  {
    enzyme: 'BsuBI',
    sequence: 'CTGCAG',
  },
  {
    enzyme: 'BsuFI',
    sequence: 'CCGG',
  },
  {
    enzyme: 'R1.BsuMI',
    sequence: 'CTCGAG',
  },
  {
    enzyme: 'R2.BsuMI',
    sequence: 'CTCGAG',
  },
  {
    enzyme: 'R3.BsuMI',
    sequence: 'CTCGAG',
  },
  {
    enzyme: 'BsuRI',
    sequence: 'GGCC',
  },
  {
    enzyme: 'BtsI (BtsIA)',
    sequence: 'GCAGTG',
  },
  {
    enzyme: 'BtsI (BtsIB)',
    sequence: 'GCAGTG',
  },
  {
    enzyme: 'BtsCI',
    sequence: 'GGATG',
  },
  {
    enzyme: 'CatHI',
    sequence: 'CTCTTC',
  },
  {
    enzyme: 'CbeI',
    sequence: 'GGCC',
  },
  {
    enzyme: 'CceI',
    sequence: 'CCGG',
  },
  {
    enzyme: 'Cce743I',
    sequence: 'GACGC',
  },
  {
    enzyme: 'CchI',
    sequence: 'CTAG',
  },
  {
    enzyme: 'CcoLI',
    sequence: 'GATC',
  },
  {
    enzyme: 'Cfr9I',
    sequence: 'CCCGGG',
  },
  {
    enzyme: 'Cfr42I',
    sequence: 'CCGCGG',
  },
  {
    enzyme: 'ClaI',
    sequence: 'ATCGAT',
  },
  {
    enzyme: 'CpaAI',
    sequence: 'CGCG',
  },
  {
    enzyme: 'CphBI',
    sequence: 'CCCGGG',
  },
  {
    enzyme: 'Csp231I',
    sequence: 'AAGCTT',
  },
  {
    enzyme: 'Csp12AI',
    sequence: 'GGATG',
  },
  {
    enzyme: 'Csp12AORFA',
    sequence: 'GGATG',
  },
  {
    enzyme: 'Csp104CI',
    sequence: 'GCATC',
  },
  {
    enzyme: 'Csp68KII',
    sequence: 'TTCGAA',
  },
  {
    enzyme: 'Csp68KIII',
    sequence: 'ATGCAT',
  },
  {
    enzyme: 'Csp68KVI',
    sequence: 'CGCG',
  },
  {
    enzyme: 'CviAI',
    sequence: 'GATC',
  },
  {
    enzyme: 'CviAII',
    sequence: 'CATG',
  },
  {
    enzyme: 'CviQI',
    sequence: 'GTAC',
  },
  {
    enzyme: 'DdsI',
    sequence: 'GGATCC',
  },
  {
    enzyme: 'DpnII',
    sequence: 'GATC',
  },
  {
    enzyme: 'DraI',
    sequence: 'TTTAAA',
  },
  {
    enzyme: 'EacI',
    sequence: 'GGATC',
  },
  {
    enzyme: 'EagI',
    sequence: 'CGGCCG',
  },
  {
    enzyme: 'Eco31I',
    sequence: 'GGTCTC',
  },
  {
    enzyme: 'Eco1524I',
    sequence: 'AGGCCT',
  },
  {
    enzyme: 'EcoGIII',
    sequence: 'CTGCAG',
  },
  {
    enzyme: 'EcoRI',
    sequence: 'GAATTC',
  },
  {
    enzyme: 'EcoRV',
    sequence: 'GATATC',
  },
  {
    enzyme: 'EcoVIII',
    sequence: 'AAGCTT',
  },
  {
    enzyme: 'Eco29kI',
    sequence: 'CCGCGG',
  },
  {
    enzyme: 'EsaWC1I',
    sequence: 'GGCC',
  },
  {
    enzyme: 'Esp3I',
    sequence: 'CGTCTC',
  },
  {
    enzyme: 'EspCI',
    sequence: 'ACCTGC',
  },
  {
    enzyme: 'FnuDI',
    sequence: 'GGCC',
  },
  {
    enzyme: 'FokI',
    sequence: 'GGATG',
  },
  {
    enzyme: 'FpsJI',
    sequence: 'CCGG',
  },
  {
    enzyme: 'FseI',
    sequence: 'GGCCGGCC',
  },
  {
    enzyme: 'FspI',
    sequence: 'TGCGCA',
  },
  {
    enzyme: 'Fsp291I',
    sequence: 'ACCTGC',
  },
  {
    enzyme: 'Gba686I',
    sequence: 'GGATCC',
  },
  {
    enzyme: 'HaeIII',
    sequence: 'GGCC',
  },
  {
    enzyme: 'HgiDII',
    sequence: 'GTCGAC',
  },
  {
    enzyme: 'HhaI',
    sequence: 'GCGC',
  },
  {
    enzyme: 'Hin4II',
    sequence: 'CCTTC',
  },
  {
    enzyme: 'HinP1I',
    sequence: 'GCGC',
  },
  {
    enzyme: 'HindIII',
    sequence: 'AAGCTT',
  },
  {
    enzyme: 'HpaI',
    sequence: 'GTTAAC',
  },
  {
    enzyme: 'HpaII',
    sequence: 'CCGG',
  },
  {
    enzyme: 'HphI',
    sequence: 'GGTGA',
  },
  {
    enzyme: 'Hpy30XI',
    sequence: 'CCATC',
  },
  {
    enzyme: 'Hpy99III',
    sequence: 'GCGC',
  },
  {
    enzyme: 'Hpy299IX',
    sequence: 'CCGG',
  },
  {
    enzyme: 'Hpy303I',
    sequence: 'GTAC',
  },
  {
    enzyme: 'Hpy312I',
    sequence: 'GTAC',
  },
  {
    enzyme: 'Hpy401I',
    sequence: 'GTAC',
  },
  {
    enzyme: 'Hpy421I',
    sequence: 'GTAC',
  },
  {
    enzyme: 'Hpy423I',
    sequence: 'GTAC',
  },
  {
    enzyme: 'Hpy471I',
    sequence: 'GTAC',
  },
  {
    enzyme: 'Hpy501I',
    sequence: 'GTAC',
  },
  {
    enzyme: 'HpyAII',
    sequence: 'GAAGA',
  },
  {
    enzyme: 'HpyAIII',
    sequence: 'GATC',
  },
  {
    enzyme: 'HpyAV',
    sequence: 'CCTTC',
  },
  {
    enzyme: 'HpyC1I',
    sequence: 'CCATC',
  },
  {
    enzyme: 'HpyCH4I',
    sequence: 'CATG',
  },
  {
    enzyme: 'HpyCH4V',
    sequence: 'TGCA',
  },
  {
    enzyme: 'HpyGII',
    sequence: 'TGCA',
  },
  {
    enzyme: 'HpyHII',
    sequence: 'GTAC',
  },
  {
    enzyme: 'HpyNSH57I',
    sequence: 'GTAC',
  },
  {
    enzyme: 'HsoI',
    sequence: 'GCGC',
  },
  {
    enzyme: 'KasI',
    sequence: 'GGCGCC',
  },
  {
    enzyme: 'KpnI',
    sequence: 'GGTACC',
  },
  {
    enzyme: 'Kpn2I',
    sequence: 'TCCGGA',
  },
  {
    enzyme: 'LlaAI',
    sequence: 'GATC',
  },
  {
    enzyme: 'LlaCI',
    sequence: 'AAGCTT',
  },
  {
    enzyme: 'LlaDI',
    sequence: 'AGTACT',
  },
  {
    enzyme: 'LlaDCHI',
    sequence: 'GATC',
  },
  {
    enzyme: 'LlaKR2I',
    sequence: 'GATC',
  },
  {
    enzyme: 'LraI',
    sequence: 'GAATTC',
  },
  {
    enzyme: 'Lsp1109I',
    sequence: 'GCAGC',
  },
  {
    enzyme: 'MboI',
    sequence: 'GATC',
  },
  {
    enzyme: 'MboII',
    sequence: 'GAAGA',
  },
  {
    enzyme: 'McaCI',
    sequence: 'CCATC',
  },
  {
    enzyme: 'McaTI',
    sequence: 'GCGCGC',
  },
  {
    enzyme: 'MhyGDL1III',
    sequence: 'GATC',
  },
  {
    enzyme: 'Mis1I',
    sequence: 'GATC',
  },
  {
    enzyme: 'MjaI',
    sequence: 'CTAG',
  },
  {
    enzyme: 'MjaIII',
    sequence: 'GATC',
  },
  {
    enzyme: 'MjaV',
    sequence: 'GTAC',
  },
  {
    enzyme: 'MluI',
    sequence: 'ACGCGT',
  },
  {
    enzyme: 'MmeII',
    sequence: 'GATC',
  },
  {
    enzyme: 'MmyCVI',
    sequence: 'CCATC',
  },
  {
    enzyme: 'MnlI',
    sequence: 'CCTC',
  },
  {
    enzyme: 'MscI',
    sequence: 'TGGCCA',
  },
  {
    enzyme: 'MseI',
    sequence: 'TTAA',
  },
  {
    enzyme: 'MspI',
    sequence: 'CCGG',
  },
  {
    enzyme: 'MthTI',
    sequence: 'GGCC',
  },
  {
    enzyme: 'MthZI',
    sequence: 'CTAG',
  },
  {
    enzyme: 'MunI',
    sequence: 'CAATTG',
  },
  {
    enzyme: 'Mva1269I',
    sequence: 'GAATGC',
  },
  {
    enzyme: 'Mva1312II',
    sequence: 'GCATC',
  },
  {
    enzyme: 'NaeI',
    sequence: 'GCCGGC',
  },
  {
    enzyme: 'NcoI',
    sequence: 'CCATGG',
  },
  {
    enzyme: 'NcuI',
    sequence: 'GAAGA',
  },
  {
    enzyme: 'NflHII',
    sequence: 'CCGG',
  },
  {
    enzyme: 'NgoAII',
    sequence: 'GGCC',
  },
  {
    enzyme: 'NgoAIII',
    sequence: 'CCGCGG',
  },
  {
    enzyme: 'NgoAIV',
    sequence: 'GCCGGC',
  },
  {
    enzyme: 'NgoAVII',
    sequence: 'GCCGC',
  },
  {
    enzyme: 'NgoAXVI',
    sequence: 'GGTGA',
  },
  {
    enzyme: 'NgoBVIII',
    sequence: 'GGTGA',
  },
  {
    enzyme: 'NgoMIII',
    sequence: 'CCGCGG',
  },
  {
    enzyme: 'NgoMIV',
    sequence: 'GCCGGC',
  },
  {
    enzyme: 'NgoMVIII',
    sequence: 'GGTGA',
  },
  {
    enzyme: 'NgoPII',
    sequence: 'GGCC',
  },
  {
    enzyme: 'NgoSII',
    sequence: 'GGCC',
  },
  {
    enzyme: 'NheI',
    sequence: 'GCTAGC',
  },
  {
    enzyme: 'NmeBI',
    sequence: 'GACGC',
  },
  {
    enzyme: 'NmeSI',
    sequence: 'AGTACT',
  },
  {
    enzyme: 'NsoJS138I',
    sequence: 'CAGCTG',
  },
  {
    enzyme: 'NspV',
    sequence: 'TTCGAA',
  },
  {
    enzyme: 'Nti1539I',
    sequence: 'ACCTGC',
  },
  {
    enzyme: 'PabI',
    sequence: 'GTAC',
  },
  {
    enzyme: 'PacI',
    sequence: 'TTAATTAA',
  },
  {
    enzyme: 'Pac25I',
    sequence: 'CCCGGG',
  },
  {
    enzyme: 'PaeR7I',
    sequence: 'CTCGAG',
  },
  {
    enzyme: 'Pam7686I',
    sequence: 'CCATGG',
  },
  {
    enzyme: 'PciI',
    sequence: 'ACATGT',
  },
  {
    enzyme: 'PenI',
    sequence: 'GCAGT',
  },
  {
    enzyme: 'PhoI',
    sequence: 'GGCC',
  },
  {
    enzyme: 'PluTI',
    sequence: 'GGCGCC',
  },
  {
    enzyme: 'PmeI',
    sequence: 'GTTTAAAC',
  },
  {
    enzyme: 'PmlI',
    sequence: 'CACGTG',
  },
  {
    enzyme: 'PspOMI',
    sequence: 'GGGCCC',
  },
  {
    enzyme: 'PvuI',
    sequence: 'CGATCG',
  },
  {
    enzyme: 'RflFI',
    sequence: 'GTCGAC',
  },
  {
    enzyme: 'RrhJ1I',
    sequence: 'GCCGGC',
  },
  {
    enzyme: 'RsaI',
    sequence: 'GTAC',
  },
  {
    enzyme: 'RshI',
    sequence: 'CGATCG',
  },
  {
    enzyme: 'Rsp241I',
    sequence: 'CGATCG',
  },
  {
    enzyme: 'RsrI',
    sequence: 'GAATTC',
  },
  {
    enzyme: 'SalI',
    sequence: 'GTCGAC',
  },
  {
    enzyme: 'SbaI',
    sequence: 'CAGCTG',
  },
  {
    enzyme: 'ScaI',
    sequence: 'AGTACT',
  },
  {
    enzyme: 'SdaI',
    sequence: 'CCTGCAGG',
  },
  {
    enzyme: 'SenpCI',
    sequence: 'CCGCGG',
  },
  {
    enzyme: 'SfaNI',
    sequence: 'GCATC',
  },
  {
    enzyme: 'SghWII',
    sequence: 'GCCGGC',
  },
  {
    enzyme: 'Sgr13350I',
    sequence: 'GAGCTC',
  },
  {
    enzyme: 'Sgr13350III',
    sequence: 'GCTCTTC',
  },
  {
    enzyme: 'SmaI',
    sequence: 'CCCGGG',
  },
  {
    enzyme: 'SmoLI',
    sequence: 'CCGG',
  },
  {
    enzyme: 'SnaBI',
    sequence: 'TACGTA',
  },
  {
    enzyme: 'SonI',
    sequence: 'ATCGAT',
  },
  {
    enzyme: 'SphI',
    sequence: 'GCATGC',
  },
  {
    enzyme: 'Spn23FI',
    sequence: 'GATC',
  },
  {
    enzyme: 'Sse9I',
    sequence: 'AATT',
  },
  {
    enzyme: 'SsoI',
    sequence: 'GAATTC',
  },
  {
    enzyme: 'SspI',
    sequence: 'AATATT',
  },
  {
    enzyme: 'Ssu211I',
    sequence: 'GATC',
  },
  {
    enzyme: 'Ssu212I',
    sequence: 'GATC',
  },
  {
    enzyme: 'Ssu220I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R1.Ssu2479I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R2.Ssu2479I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R1.Ssu4109I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R2.Ssu4109I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R1.Ssu4961I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R2.Ssu4961I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R1.Ssu8074I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R2.Ssu8074I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R1.Ssu11318I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R2.Ssu11318I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R1.SsuDAT1I',
    sequence: 'GATC',
  },
  {
    enzyme: 'R2.SsuDAT1I',
    sequence: 'GATC',
  },
  {
    enzyme: 'Sth368I',
    sequence: 'GATC',
  },
  {
    enzyme: 'StsI',
    sequence: 'GGATG',
  },
  {
    enzyme: 'SuaI',
    sequence: 'GGCC',
  },
  {
    enzyme: 'Sxa3845I',
    sequence: 'GCCGGC',
  },
  {
    enzyme: 'TdeII',
    sequence: 'CTCTTC',
  },
  {
    enzyme: 'TfiTok6A1I',
    sequence: 'TCGA',
  },
  {
    enzyme: 'TflI',
    sequence: 'TCGA',
  },
  {
    enzyme: 'ThaI',
    sequence: 'CGCG',
  },
  {
    enzyme: 'TliI',
    sequence: 'CTCGAG',
  },
  {
    enzyme: 'TmaI',
    sequence: 'CGCG',
  },
  {
    enzyme: 'TneDI',
    sequence: 'CGCG',
  },
  {
    enzyme: 'Tsp32I',
    sequence: 'TCGA',
  },
  {
    enzyme: 'Tsp509I',
    sequence: 'AATT',
  },
  {
    enzyme: 'TspDTI',
    sequence: 'ATGAA',
  },
  {
    enzyme: 'TspMI',
    sequence: 'CCCGGG',
  },
  {
    enzyme: 'TthHB8I',
    sequence: 'TCGA',
  },
  {
    enzyme: 'Tvu2HI',
    sequence: 'GGCC',
  },
  {
    enzyme: 'UpaP162I',
    sequence: 'CATG',
  },
  {
    enzyme: 'Van91II',
    sequence: 'GAATTC',
  },
  {
    enzyme: 'XamI',
    sequence: 'GTCGAC',
  },
  {
    enzyme: 'XbaI',
    sequence: 'TCTAGA',
  },
  {
    enzyme: 'XcyI',
    sequence: 'CCCGGG',
  },
  {
    enzyme: 'XhoI',
    sequence: 'CTCGAG',
  },
  {
    enzyme: 'XmaIII',
    sequence: 'CGGCCG',
  },
  {
    enzyme: 'XorKI',
    sequence: 'CGATCG',
  },
  {
    enzyme: 'XorKII',
    sequence: 'CTGCAG',
  },
  {
    enzyme: 'XphI',
    sequence: 'CTGCAG',
  },
  {
    enzyme: 'XveI',
    sequence: 'CTGCAG',
  },
  {
    enzyme: 'NotI',
    sequence: 'GCGGCCGC',
  },
  {
    enzyme: 'SacI',
    sequence: 'GAGCTC',
  },
  {
    enzyme: 'PstI',
    sequence: 'CTGCAG',
  },
  {
    enzyme: 'SpeI',
    sequence: 'ACTAGT',
  },
  {
    enzyme: 'PvuII',
    sequence: 'CAGCTG',
  },
  {
    enzyme: 'NdeI',
    sequence: 'CATATG',
  },
  {
    enzyme: 'MfeI',
    sequence: 'CAATTG',
  },
  {
    enzyme: 'NarI',
    sequence: 'GGCGCC',
  },
  {
    enzyme: 'NsiI',
    sequence: 'ATGCAT',
  },
  {
    enzyme: 'SacII',
    sequence: 'CCGCGG',
  },
  {
    enzyme: 'SbfI',
    sequence: 'CCTGCAGG',
  },
  {
    enzyme: 'StuI',
    sequence: 'AGGCCT',
  },
  {
    enzyme: 'TaqI',
    sequence: 'TCGA',
  },
  {
    enzyme: 'XmaI',
    sequence: 'CCCGGG',
  },
  {
    enzyme: 'AvrII',
    sequence: 'CCTAGG',
  },
  {
    enzyme: 'BbsI',
    sequence: 'GAAGAC',
  },
  {
    enzyme: 'BmgBI',
    sequence: 'CACGTC',
  },
  {
    enzyme: 'BseRI',
    sequence: 'GAGGAG',
  },
  {
    enzyme: 'BseYI',
    sequence: 'CCCAGC',
  },
  {
    enzyme: 'BsmFI',
    sequence: 'GGGAC',
  },
  {
    enzyme: 'BspDI',
    sequence: 'ATCGAT',
  },
  {
    enzyme: 'BspMI',
    sequence: 'ACCTGC',
  },
  {
    enzyme: 'BspQI',
    sequence: 'GCTCTTC',
  },
  {
    enzyme: 'BsrBI',
    sequence: 'CCGCTC',
  },
  {
    enzyme: 'BsrDI',
    sequence: 'GCAATG',
  },
  {
    enzyme: 'BsrGI',
    sequence: 'TGTACA',
  },
  {
    enzyme: 'BssHII',
    sequence: 'GCGCGC',
  },
  {
    enzyme: 'BssSI',
    sequence: 'CACGAG',
  },
  {
    enzyme: 'BstBI',
    sequence: 'TTTCGAA',
  },
  {
    enzyme: 'BstZ17I',
    sequence: 'GTATAC',
  },
  {
    enzyme: 'EarI',
    sequence: 'CTCTTC',
  },
  {
    enzyme: 'EciI',
    sequence: 'GGCGGA',
  },
  {
    enzyme: 'EcoP15I',
    sequence: 'CAGCAG',
  },
  {
    enzyme: 'FauI',
    sequence: 'CCCGC',
  },
  {
    enzyme: 'HgaI',
    sequence: 'GACGC',
  },
  {
    enzyme: 'MlyI',
    sequence: 'GAGTC',
  },
  {
    enzyme: 'NlaIII',
    sequence: 'CATG',
  },
  {
    enzyme: 'NmeAIII',
    sequence: 'GCCGGC',
  },
  {
    enzyme: 'PleI',
    sequence: 'GAGTC',
  },
  {
    enzyme: 'Sau3AI',
    sequence: 'GATC',
  },
  {
    enzyme: 'SgrDI',
    sequence: 'CGTCGACG',
  },
  {
    enzyme: 'SnaI',
    sequence: 'GTATAC',
  },
  {
    enzyme: 'ZraI',
    sequence: 'GACGTC',
  },
  {
    enzyme: 'NruI',
    sequence: 'TCGCGA',
  },
]
