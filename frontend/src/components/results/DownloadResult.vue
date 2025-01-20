<script setup lang="ts">
import type { TuningOutput, TunedSequence } from '@/lib/interfaces'
import { downloadFile } from '@/lib/helpers'

const props = defineProps<TuningOutput>()

function foldSequence(sequence: string) {
  return sequence.replace(/(.{80})/g, '$1\n')
}

function formatToFasta(tunedSequences: Array<TunedSequence>) {
  return tunedSequences
    .map(e => `>${e.name}\n${foldSequence(e.output)}`)
    .join('\n\n')
}

function downloadFasta() {
  downloadFile(
    formatToFasta(props.tuned_sequences),
    `tuned_sequences_${props.result.mode}.fasta`,
  )
}
</script>

<template>
  <button @click="downloadFasta">Download output FASTA file</button>
</template>

<style scoped>
button {
  width: 100%;
  margin: 1em 0 1.5em 0;
}
</style>
