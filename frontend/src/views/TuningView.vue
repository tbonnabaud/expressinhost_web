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

function goBackToForm() {
  tuningOutput.value = {} as TuningOutput
}
</script>

<template>
  <main class="container">
    <h1>Sequence tuning</h1>

    <hr />

    <div v-if="tuningOutput.result">
      <button id="goBackToFormButton" @click="goBackToForm">
        Go back to form
      </button>

      <ResultContent
        :result="tuningOutput.result"
        :tuned_sequences="tuningOutput.tuned_sequences"
      />
    </div>

    <TuningForm v-else @submit="handleSubmit" />
  </main>
</template>

<style scoped>
#goBackToFormButton {
  float: right;
}
</style>
