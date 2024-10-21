<script setup lang="ts">
import { computed } from 'vue'
import type { TunedSequence } from '@/lib/interfaces'

const props = defineProps<{
  tunedSequence: TunedSequence
  open?: boolean
}>()

const seqComparison = computed(() =>
  colorSequences(props.tunedSequence.input, props.tunedSequence.output),
)

function colorSequences(inputSequence: string, outputSequence: string) {
  const inputCodonArray = inputSequence.match(/.{3}/g) || []
  const outputCodonArray = outputSequence.match(/.{3}/g) || []

  return [
    inputCodonArray
      .map((item, i) => {
        if (item != outputCodonArray[i]) {
          return `<span class="in-seq">${item}</span>`
        }
        return item
      })
      .join(''),
    outputCodonArray
      .map((item, i) => {
        if (item != inputCodonArray[i]) {
          return `<span class="out-seq">${item}</span>`
        }
        return item
      })
      .join(''),
  ]
}
</script>

<template>
  <details :open="open">
    <summary>{{ tunedSequence.name }}</summary>

    <p>&rarr; Similarity:{{ tunedSequence.identity_percentage.toFixed(2) }}%</p>

    <div class="flex-container">
      <div class="sequence-group-label">
        <label>Input: </label>
        <label>Output: </label>
      </div>

      <div class="sequence-group">
        <div class="sequence" v-html="seqComparison[0]"></div>
        <div class="sequence" v-html="seqComparison[1]"></div>
      </div>
    </div>
  </details>

  <hr />
</template>

<style scoped>
.sequence {
  height: 3em;
  white-space: nowrap;
  font-family: monospace;
  font-size: 1.5em;
}

.sequence-group {
  overflow-x: scroll;
  margin: 2em 0 0 1em;
}

.sequence-group-label {
  margin: 2em 0;
}

.sequence-group-label label {
  height: 3em;
  font-weight: bold;
  text-decoration: underline;
}

details summary {
  font-size: 1.2em;
}
</style>
