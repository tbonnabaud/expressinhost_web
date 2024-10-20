<script setup lang="ts">
import { ref } from 'vue'
import TuningForm from '@/components/TuningForm.vue'
import ResultContent from '@/components/ResultContent.vue'
import type { TuningOutput } from '@/lib/interfaces'

const tuningOutput = ref({} as TuningOutput)

function handleSubmit(ouput: TuningOutput) {
  tuningOutput.value = ouput
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
