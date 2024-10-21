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
    <summary>{{ model || 'Select one...' }}</summary>

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
        <li v-for="option in filteredOptions" :key="option">
          <label @click="closeDropdown">
            <input type="radio" v-model="model" :value="option" />
            {{ option }}
          </label>
        </li>
      </div>
    </ul>
  </details>
</template>

<style scoped>
.option-list {
  max-height: 30em;
  bottom: 0;
  overflow-y: scroll;
}

input {
  margin-bottom: 0;
  display: none;
}
</style>
