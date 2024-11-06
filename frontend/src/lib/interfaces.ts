export interface User {
  id: string
  creation_date: string
  email: string
  role: string
  full_name: string
}

export interface UserForm {
  email: string
  password: string
  full_name: string
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
  nucleotide_file_content: string
  clustal_file_content: string | null
  host_codon_table_name: string
  sequences_native_codon_tables: Record<string, string>
  mode: string
  slow_speed_threshold: number
  conservation_threshold?: number
}

export interface CodonTable {
  name: string
  organism: string
  custom: boolean
}

export interface Result {
  id: string | null
  user_id: string | null
  creation_date: string
  host_codon_table_name: string
  sequences_native_codon_tables: Record<string, string>
  mode: string
  slow_speed_threshold: number
  conservation_threshold: number | null
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
