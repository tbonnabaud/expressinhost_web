<script setup lang="ts">
import { ref, watch } from 'vue'
import { UTR_EXAMPLE } from '@/lib/referentials'
import { checkUtrSequence } from '@/lib/checkers'
import WithAlertError from '@/components/WithAlertError.vue'

const model = defineModel<string>()
const utr = ref('')
const errors = ref<string[]>([])

watch(utr, value => {
  errors.value = checkUtrSequence(value)
  model.value = value
})

function addExampleUTR() {
  utr.value = UTR_EXAMPLE
}
</script>

<template>
  <WithAlertError :errors="errors">
    <label id="utrSequenceLabel" for="utrSequence">5' UTR sequence</label>
    <textarea
      id="utrSequence"
      v-model="utr"
      placeholder="Put 5' UTR sequence."
      rows="3"
      spellcheck="false"
    ></textarea>
  </WithAlertError>

  <button
    id="addUtrExampleButton"
    @click.prevent="addExampleUTR"
    class="outline secondary"
  >
    Fill with an example
  </button>
</template>

<style scoped>
#addUtrExampleButton {
  font-size: 90%;
  padding: 5px 10px;
}
</style>
