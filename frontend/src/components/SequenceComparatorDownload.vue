<script setup lang="ts">
import type { TunedSequence } from '@/lib/interfaces'
import { downloadFile } from '@/lib/helpers'
import { API } from '@/lib/api'
import JSZip from 'jszip'
import { ref } from 'vue'

const props = defineProps<{
  comparedSequences: TunedSequence
}>()

const isLoading = ref(false)

function formatToFasta(comparedSequences: TunedSequence) {
  const seq1 = `>Sequence1 ${comparedSequences.name}\n${comparedSequences.input}`
  const seq2 = `>Sequence2 ${comparedSequences.name}\n${comparedSequences.output}`

  return `${seq1}\n\n${seq2}`
}

/**
 * Generates a CSV formatted string from input and output codon arrays and their corresponding profiles.
 *
 * @param inputCodonArray - An array of input codon strings.
 * @param inputProfile - An array of numerical values representing the profile of the input codons.
 * @param outputCodonArray - An array of output codon strings.
 * @param outputProfile - An array of numerical values representing the profile of the output codons.
 * @param profileType - The type of profile being represented, either 'speed' or 'rank'.
 * @returns A CSV formatted string with columns for index, input codon, input profile, output codon, and output profile.
 * @throws Will throw an error if the inputProfile and outputProfile arrays do not have the same length.
 */
function formatProfileCSV(
  inputCodonArray: string[],
  inputProfile: number[] | null | undefined,
  outputCodonArray: string[],
  outputProfile: number[],
  profileType: 'speed' | 'rank',
): string {
  const csvRows: string[] = []

  if (inputProfile) {
    if (inputProfile.length !== outputProfile.length) {
      alert('Input and output arrays must have the same length')
      throw new Error('Input and output arrays must have the same length')
    }

    // Create the CSV header
    csvRows.push(
      `index,input_codon,input_${profileType},output_codon,output_${profileType}`,
    )

    // Create the CSV rows
    for (let i = 0; i < inputProfile.length; i++) {
      const row = `${i},${inputCodonArray[i]},${inputProfile[i]},${outputCodonArray[i]},${outputProfile[i]}`
      csvRows.push(row)
    }
  } else {
    // Create the CSV header
    csvRows.push(`index,output_codon,output_${profileType}`)

    // Create the CSV rows
    for (let i = 0; i < outputProfile.length; i++) {
      const row = `${i},${outputCodonArray[i]},${outputProfile[i]}`
      csvRows.push(row)
    }
  }

  // Join all rows into a single CSV string
  return csvRows.join('\n')
}

async function fetchCodonTableFullName(id: string) {
  const [data, error] = await API.codonTables.get(id)

  if (!error && data) {
    return `${data.organism} - ${data.name}`
  } else {
    return null
  }
}

async function downloadZip() {
  isLoading.value = true
  const zip = new JSZip()
  // FASTA
  const fastaContent = formatToFasta(props.comparedSequences)
  const fastaFileName = `compared_sequences_${props.comparedSequences.name}.fasta`
  // ZIP file name
  const zipFileName = `${props.comparedSequences.name}.zip`

  // Add the files into the archive
  zip.file(fastaFileName, fastaContent)

  // Match only ID to avoid a too long name
  const matchId = props.comparedSequences.name.match(/^\S+/)

  if (matchId) {
    const seqId = matchId[0]
    const inputCodonArray = props.comparedSequences.input.match(/.{3}/g) || []
    const outputCodonArray = props.comparedSequences.output.match(/.{3}/g) || []
    // Add speed profiles
    zip.file(
      `speeds/${seqId}_speed_profiles.csv`,
      formatProfileCSV(
        inputCodonArray,
        props.comparedSequences.input_profiles?.speed,
        outputCodonArray,
        props.comparedSequences.output_profiles.speed,
        'speed',
      ),
    )
    // Add rank profiles
    zip.file(
      `ranks/${seqId}_rank_profiles.csv`,
      formatProfileCSV(
        inputCodonArray,
        props.comparedSequences.input_profiles?.rank,
        outputCodonArray,
        props.comparedSequences.output_profiles.rank,
        'rank',
      ),
    )
  } else {
    console.warn(`ID not found for sequence: ${props.comparedSequences.name}`)
  }

  const zipBlob = await zip.generateAsync({ type: 'blob' })
  downloadFile(zipBlob, zipFileName)
  isLoading.value = false
}
</script>

<template>
  <button @click="downloadZip" :aria-busy="isLoading" :disabled="isLoading">
    Download zipped output
  </button>
</template>

<style scoped>
button {
  width: 100%;
  margin: 1em 0 1.5em 0;
}
</style>
