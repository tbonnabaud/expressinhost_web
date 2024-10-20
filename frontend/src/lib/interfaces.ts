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
