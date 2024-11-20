<script setup lang="ts">
import { computed, ref, useTemplateRef } from 'vue'
import type { CodonTable } from '@/lib/interfaces'

const props = defineProps<{
  options: Array<CodonTable>
}>()
const model = defineModel<CodonTable | null>()
const filter = ref('')
const dropdownRef = useTemplateRef('dropdown')

const filteredOptions = computed(() => {
  const lowerCaseFilter = filter.value.toLowerCase()
  return props.options.filter(
    e =>
      e.name.toLowerCase().includes(lowerCaseFilter) ||
      e.organism.toLowerCase().includes(lowerCaseFilter),
  )
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
    <summary class="summary-select">
      <template v-if="model">
        <i>{{ model.organism }}</i> - {{ model.name }}
      </template>

      <template v-else>Select one...</template>
    </summary>

    <ul>
      <li>
        <input
          type="search"
          placeholder="Filter..."
          ref="filter-input"
          v-model="filter"
        />
      </li>
      <div class="option-list">
        <li>
          <label @click="closeDropdown">
            <input type="radio" v-model="model" :value="null" />
            (None)
          </label>
        </li>
        <li v-for="option in filteredOptions" :key="option.id">
          <label @click="closeDropdown">
            <input type="radio" v-model="model" :value="option" />
            <i>{{ option.organism }}</i> - {{ option.name }}
          </label>
        </li>
      </div>
    </ul>
  </details>
</template>

<style scoped>
.option-list {
  max-height: 30em;
  overflow-y: scroll;
}

input[type='search'] {
  margin-bottom: 0;
}

input[type='radio'] {
  display: none;
}

.dropdown .summary-select {
  height: fit-content;
}
</style>
