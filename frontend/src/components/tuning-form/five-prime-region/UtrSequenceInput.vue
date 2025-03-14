<script setup lang="ts">
import { ref, watch } from 'vue'
import { UTR_EXAMPLE } from '@/lib/referentials'

const model = defineModel<string>()
const utr = ref('')

watch(utr, value => {
  if (value) {
    model.value = value
  }
})

function addExampleUTR() {
  utr.value = UTR_EXAMPLE
}
</script>

<template>
  <label id="utrSequenceLabel" for="utrSequence">UTR sequence</label>
  <textarea
    id="utrSequence"
    v-model="model"
    placeholder="Put UTR sequence"
    rows="3"
    spellcheck="false"
  ></textarea>
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
