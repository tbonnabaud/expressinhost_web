export interface RunTrainingForm {
  nucleotide_file_content: string
  clustal_file_content?: string
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
  creation_date: string
  host_codon_table_name: string
  sequences_native_codon_tables: Record<string, string>
  mode: string
  slow_speed_threshold: number
  conservation_threshold?: number
}

export interface TunedSequence {
  name: string
  input: string
  output: string
  identity_percentage: number
}

export interface TuningOutput {
  result: Result
  tuned_sequences: TunedSequence
}
