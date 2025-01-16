<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import type { CodonTable } from '@/lib/interfaces'

const DEFAULT_NUMBER_TO_SHOW = 100

const props = defineProps<{
  options: Array<CodonTable>
}>()
const model = defineModel<CodonTable | null>()
const filter = ref('')
const collapseDropdown = ref(true)
const optionsToShow = ref(DEFAULT_NUMBER_TO_SHOW)

const filteredOptions = computed(() => {
  const lowerCaseFilter = filter.value.toLowerCase()
  return props.options
    .filter(
      e =>
        e.name.toLowerCase().includes(lowerCaseFilter) ||
        e.organism.toLowerCase().includes(lowerCaseFilter),
    )
    .slice(0, optionsToShow.value)
})

watch(filter, () => {
  optionsToShow.value = DEFAULT_NUMBER_TO_SHOW
})

watch(collapseDropdown, () => {
  optionsToShow.value = DEFAULT_NUMBER_TO_SHOW
})

function handleScroll(event: Event) {
  const target = event.target as HTMLElement
  const { scrollTop, scrollHeight, clientHeight } = target

  if (scrollTop + clientHeight >= scrollHeight * 0.9) {
    optionsToShow.value += DEFAULT_NUMBER_TO_SHOW
  }
}

function closeDropdown() {
  collapseDropdown.value = true
  filter.value = ''
}
</script>

<template>
  <details class="dropdown" :open="!collapseDropdown">
    <summary
      class="summary-select"
      @click.prevent="collapseDropdown = !collapseDropdown"
    >
      <template v-if="model">
        <i>{{ model.organism }}</i> - {{ model.name }}
      </template>

      <template v-else>Select one...</template>
    </summary>

    <ul v-if="!collapseDropdown">
      <li>
        <input type="search" placeholder="Filter..." v-model="filter" />
      </li>
      <div class="option-list" @scroll="handleScroll">
        <li>
          <label>
            <input
              type="radio"
              v-model="model"
              :value="null"
              @change="closeDropdown"
            />
            (None)
          </label>
        </li>
        <li v-for="option in filteredOptions" :key="option.id">
          <label>
            <input
              type="radio"
              v-model="model"
              :value="option"
              @change="closeDropdown"
            />
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
