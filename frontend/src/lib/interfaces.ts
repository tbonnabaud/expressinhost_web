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

export interface RunTrainingForm {
  name: string
  nucleotide_file_content: string
  clustal_file_content: string | null
  host_codon_table_id: string
  sequences_native_codon_tables: Record<string, string>
  mode: string
  slow_speed_threshold: number
  conservation_threshold?: number
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
  input_profiles: Profiles
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
