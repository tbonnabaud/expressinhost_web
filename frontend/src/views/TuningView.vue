<script setup lang="ts">
import { ref } from 'vue'
import type { TuningOutput } from '@/lib/interfaces'
import TuningForm from '@/components/TuningForm.vue'
import ResultContent from '@/components/results/ResultContent.vue'
import BaseModal from '@/components/BaseModal.vue'

const tuningOutput = ref({} as TuningOutput)
const showBackToFormModal = ref(false)

function handleSubmit(output: TuningOutput) {
  if (output) {
    tuningOutput.value = output
    window.scrollTo(0, 0)
  }
}

function goBackToForm() {
  localStorage.removeItem('tuningJobId')
  tuningOutput.value = {} as TuningOutput
}
</script>

<template>
  <main class="container">
    <h1>Tuning</h1>

    <hr />

    <div v-if="tuningOutput.result">
      <BaseModal
        :open="showBackToFormModal"
        title="Confirmation"
        @close="showBackToFormModal = false"
      >
        Are you sure to want to go back to a new form?

        <footer>
          <button class="secondary" @click="showBackToFormModal = false">
            Cancel
          </button>
          <button id="goBackToFormButton" @click="goBackToForm">Confirm</button>
        </footer>
      </BaseModal>

      <button
        id="showGoBackToFormModalButton"
        @click="showBackToFormModal = true"
      >
        Go back to a new form
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
#showGoBackToFormModalButton {
  float: right;
}
</style>
