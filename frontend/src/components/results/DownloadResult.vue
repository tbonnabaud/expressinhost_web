<script setup lang="ts">
import type { TuningOutput, TunedSequence } from '@/lib/interfaces'
import { downloadFile } from '@/lib/helpers'
import JSZip from 'jszip'

const props = defineProps<TuningOutput>()

function foldSequence(sequence: string) {
  return sequence.replace(/(.{80})/g, '$1\n')
}

function formatToFasta(tunedSequences: Array<TunedSequence>) {
  return tunedSequences
    .map(e => `>${e.name}\n${foldSequence(e.output)}`)
    .join('\n\n')
}

/**
 * Returns a CSV formatted string with columns for codon_index, input, and output.
 *
 * @param input - Array of input numbers.
 * @param output - Array of output numbers.
 * @returns {string} The CSV formatted string.
 */
function formatProfileCSV(input: number[], output: number[]): string {
  if (input.length !== output.length) {
    throw new Error('Input and output arrays must have the same length')
  }

  // Create the CSV header
  const csvRows = ['codon_index,input,output']

  // Create the CSV rows
  for (let i = 0; i < input.length; i++) {
    const row = `${i},${input[i]},${output[i]}`
    csvRows.push(row)
  }

  // Join all rows into a single CSV string
  return csvRows.join('\n')
}

async function downloadZip() {
  const zip = new JSZip()
  // FASTA
  const fastaContent = formatToFasta(props.tuned_sequences)
  const fastaFileName = `tuned_sequences_${props.result.mode}.fasta`
  // ZIP file name
  const zipFileName = `${props.result.name}.zip`

  // Add the files into the archive
  zip.file(fastaFileName, fastaContent)

  for (const tunedSequence of props.tuned_sequences) {
    // Match only ID to avoid a too long name
    const matchId = tunedSequence.name.match(/^\S+/)

    if (matchId) {
      const seqId = matchId[0]
      // Add speed profiles
      zip.file(
        `speeds/${seqId}_speed_profiles.csv`,
        formatProfileCSV(
          tunedSequence.input_profiles.speed,
          tunedSequence.output_profiles.speed,
        ),
      )
      // Add rank profiles
      zip.file(
        `ranks/${seqId}_rank_profiles.csv`,
        formatProfileCSV(
          tunedSequence.input_profiles.rank,
          tunedSequence.output_profiles.rank,
        ),
      )
    } else {
      console.warn(`ID not found for sequence: ${tunedSequence.name}`)
    }
  }

  const zipBlob = await zip.generateAsync({ type: 'blob' })
  downloadFile(zipBlob, zipFileName)
}
</script>

<template>
  <button @click="downloadZip">Download zipped output</button>
</template>

<style scoped>
button {
  width: 100%;
  margin: 1em 0 1.5em 0;
}
</style>
