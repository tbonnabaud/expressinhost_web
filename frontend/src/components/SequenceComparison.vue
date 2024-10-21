<script setup lang="ts">
import { computed } from 'vue'
import type { TunedSequence } from '@/lib/interfaces'

const props = defineProps<TunedSequence>()

const seqComparison = computed(() => colorSequences(props.input, props.output))

function colorSequences(inputSequence: string, outputSequence: string) {
  const inputCodonArray = inputSequence.match(/.{3}/g) || []
  const outputCodonArray = outputSequence.match(/.{3}/g) || []

  return [
    inputCodonArray
      .map((item, i) => {
        if (item == outputCodonArray[i]) {
          return `<span class="in-seq">${item}</span>`
        }
        return item
      })
      .join(''),
    outputCodonArray
      .map((item, i) => {
        if (item == inputCodonArray[i]) {
          return `<span class="out-seq">${item}</span>`
        }
        return item
      })
      .join(''),
  ]
}
</script>

<template>
  <h3>{{ name }}</h3>

  <p><strong>-> Similarity:</strong> {{ identity_percentage.toFixed(2) }}%</p>

  <div class="sequence-group">
    <div class="sequence" v-html="seqComparison[0]"></div>
    <div class="sequence" v-html="seqComparison[1]"></div>
  </div>

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
  margin: 2em 0;
}
</style>
