import type { CodonTranslation } from './interfaces'
import { groupByAminoAcid } from './helpers'

export const MODE_LABEL_MAPPING: Record<string, string> = {
  direct_mapping: 'Direct mapping',
  optimisation_and_conservation_1: 'Optimisation and conservation 1',
  optimisation_and_conservation_2: 'Optimisation and conservation 2',
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

export const CODON_LIST = BASE_CODON_TABLE.map(e => e.codon)

export const CODON_MAPPING = Object.fromEntries(
  BASE_CODON_TABLE.map(e => [e.codon, e]),
)

export const AMINO_ACID_MAPPING = groupByAminoAcid(DEFAULT_CODON_TABLE)
