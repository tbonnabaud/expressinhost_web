<script setup lang="ts">
import { computed, ref, useTemplateRef } from 'vue'

const props = defineProps<{
  options: Array<string>
}>()
const model = defineModel()
const filter = ref('')
const dropdownRef = useTemplateRef('dropdown')

const filteredOptions = computed(() => {
  const lowerCaseFilter = filter.value.toLowerCase()
  return props.options.filter(e => e.toLowerCase().includes(lowerCaseFilter))
})

function closeDropdown() {
  if (dropdownRef.value) {
    dropdownRef.value.open = false
    filter.value = ''
  }
}
</script>

<template>
  <details class="dropdown" ref="dropdown">
    <summary>{{ model }}</summary>
    <ul class="option-list">
      <li>
        <input
          type="search"
          placeholder="Filter..."
          ref="filter-input"
          v-model="filter"
        />
      </li>
      <li v-for="option in filteredOptions" :key="option">
        <label @click="closeDropdown">
          <input type="radio" v-model="model" :value="option" />
          {{ option }}
        </label>
      </li>
    </ul>
  </details>
</template>

<style scoped></style>
