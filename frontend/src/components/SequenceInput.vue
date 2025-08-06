<script setup lang="ts">
import WithAlertError from './WithAlertError.vue'

const { rows = 5, required } = defineProps<{
  errors?: string[]
  rows?: number
  required?: boolean
}>()
const model = defineModel<string>()
</script>

<template>
  <WithAlertError :errors="errors">
    <textarea
      name=""
      class="sequence-input"
      :rows="rows"
      spellcheck="false"
      v-model="model"
      :required="required"
    ></textarea>
  </WithAlertError>
</template>

<style scoped>
textarea {
  font-family: 'Courier New', Courier, monospace;
}
</style>
