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

async function downloadFasta() {
  const zip = new JSZip()
  const fastaContent = formatToFasta(props.tuned_sequences)
  const fastaFileName = `tuned_sequences_${props.result.mode}.fasta`
  const zipFileName = `${props.result.name}.zip`

  zip.file(fastaFileName, fastaContent)

  const zipBlob = await zip.generateAsync({ type: 'blob' })
  downloadFile(zipBlob, zipFileName)
}
</script>

<template>
  <button @click="downloadFasta">Download zipped output</button>
</template>

<style scoped>
button {
  width: 100%;
  margin: 1em 0 1.5em 0;
}
</style>
