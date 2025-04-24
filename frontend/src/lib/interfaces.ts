import type { Component } from 'vue'

export interface User {
  id: string
  creation_date: string
  email: string
  role: string
  full_name: string
  contact_consent: boolean
}

export interface UserForm {
  email: string
  password: string
  full_name: string
  contact_consent: boolean
}

export interface UserProfileForm {
  email: string
  full_name: string
  contact_consent: boolean
}

export interface UserRoleForm {
  role: string
}

export interface UserPasswordForm {
  reset_token: string
  password: string
}

export interface UserLogin {
  username: string
  password: string
}

export interface Token {
  access_token: string
  token_type: string
}

export enum TuningModeName {
  DIRECT_MAPPING = 'direct_mapping',
  OPTIMISATION_AND_CONSERVATION_1 = 'optimisation_and_conservation_1',
  OPTIMISATION_AND_CONSERVATION_2 = 'optimisation_and_conservation_2',
  PROTEIN_STRUCTURE_ANALYSIS = 'protein_structure_analysis',
}

export interface PartialUntuningMode {
  mode: 'partial_untuning'
  untuned_codon_number: number
}

// To use OSTIR
export interface FineTuningMode {
  mode: 'fine_tuning'
  codon_window_size: number
  utr: string
}

export interface SlowedDownMode {
  mode: 'slowed_down'
  slowed_down_codon_number: number
}

export type FivePrimeRegionTuningMode =
  | PartialUntuningMode
  | FineTuningMode
  | SlowedDownMode

export interface RestrictionSite {
  enzyme: string
  sequence: string
}

export interface RunTrainingForm {
  name: string
  nucleotide_file_content: string
  pdb_file_content: string | null
  clustal_file_content: string | null
  host_codon_table_id: string
  sequences_native_codon_tables: Record<string, string>
  mode: TuningModeName
  slow_speed_threshold: number
  conservation_threshold: number | null
  five_prime_region_tuning: FivePrimeRegionTuningMode | null
  restriction_sites: RestrictionSite[]
  send_email: boolean
}

export interface CodonTable {
  id: string
  user_id: string | null
  creation_date: string
  name: string
  organism: string
  source: string | null
}

export interface CodonTranslation {
  // codon_table_id: string
  codon: string
  anticodon: string
  amino_acid: string
  trna_gcn: number
  wobble_codon: string
  wobble_rate: number
}

export interface CodonTableForm {
  name: string
  organism: string
  translations: Array<CodonTranslation>
}

export interface Result {
  id: string | null
  user_id: string | null
  creation_date: string
  name: string
  host_codon_table_id: string
  sequences_native_codon_tables: Record<string, string>
  mode: string
  slow_speed_threshold: number
  conservation_threshold: number | null
  host_codon_table: CodonTable
  five_prime_region_tuning: FivePrimeRegionTuningMode | null
  restriction_sites: RestrictionSite[] | null
}

export interface ResultWithId {
  id: string
  user_id: string | null
  creation_date: string
  name: string
  host_codon_table_id: string
  sequences_native_codon_tables: Record<string, string>
  mode: string
  slow_speed_threshold: number
  conservation_threshold: number | null
  host_codon_table: CodonTable
  five_prime_region_tuning: FivePrimeRegionTuningMode | null
  restriction_sites: RestrictionSite[] | null
}

export interface Profiles {
  speed: Array<number>
  rank: Array<number>
}

export interface TunedSequence {
  id: string | null
  result_id: string | null
  name: string
  input: string
  output: string
  identity_percentage: number
  input_profiles: Profiles | null
  output_profiles: Profiles
}

export interface TuningOutput {
  result: Result
  tuned_sequences: Array<TunedSequence>
}

export interface ComponentMeta {
  name: string
  component: Component
}

export interface LastWebScraping {
  source: string
  release: string
  scraping_date: string
}

export interface RunInfo {
  id: string
  creation_date: string
  duration: string
  sequence_number: number
  mode: string
  slow_speed_threshold: number
  conservation_threshold: number | null
}

export interface RunInfoDurationStats {
  min_duration: string
  avg_duration: string
  max_duration: string
}

export interface RunInfoSequenceNumberStats {
  min_sequence_number: number
  avg_sequence_number: number
  max_sequence_number: number
}

export interface SequenceComparatorForm {
  sequence1: string
  sequence2: string
  host_codon_table_id: string
  slow_speed_threshold: number
}
