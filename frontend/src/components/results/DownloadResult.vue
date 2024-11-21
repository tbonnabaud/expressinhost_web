<script setup lang="ts">
import type { TuningOutput, TunedSequence } from '@/lib/interfaces'

const props = defineProps<TuningOutput>()

function foldSequence(sequence: string) {
  return sequence.replace(/(.{80})/g, '$1\n')
}

function formatToFasta(tunedSequences: Array<TunedSequence>) {
  return tunedSequences
    .map(e => `>${e.name}\n${foldSequence(e.output)}`)
    .join('\n\n')
}

function downloadFile() {
  // Creating a Blob from the data
  const blob = new Blob([formatToFasta(props.tuned_sequences)], {
    type: 'text/plain',
  })
  const url = URL.createObjectURL(blob)

  // Creating a temporary link element
  const link = document.createElement('a')
  link.href = url
  link.download = `tuned_sequences_${props.result.mode}.fasta`

  // Append to the body, click and remove it
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)

  // Revoke the object URL after the download
  URL.revokeObjectURL(url)
}
</script>

<template>
  <button @click="downloadFile">Download output FASTA file</button>
</template>

<style scoped>
button {
  width: 100%;
  margin: 1em 0 1.5em 0;
}
</style>
