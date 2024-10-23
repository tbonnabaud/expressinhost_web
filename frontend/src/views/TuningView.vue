<script setup lang="ts">
import { ref } from 'vue'
import type { TuningOutput } from '@/lib/interfaces'
import TuningForm from '@/components/TuningForm.vue'
import ResultContent from '@/components/results/ResultContent.vue'

const tuningOutput = ref({} as TuningOutput)

function handleSubmit(output: TuningOutput) {
  if (output) {
    tuningOutput.value = output
  }
}
</script>

<template>
  <main class="container">
    <h1>Sequence tuning</h1>

    <hr />

    <ResultContent
      v-if="tuningOutput.result"
      :result="tuningOutput.result"
      :tuned_sequences="tuningOutput.tuned_sequences"
    />

    <TuningForm v-else @submit="handleSubmit" />
  </main>
</template>

<style scoped>
@media (min-width: 1024px) {
  .container {
    padding-left: 15em;
    padding-right: 15em;
  }
}
</style>
