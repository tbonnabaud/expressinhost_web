<script setup lang="ts">
import { computed, ref, watch, onUnmounted } from 'vue'
import type { CodonTable } from '@/lib/interfaces'

const DEFAULT_NUMBER_TO_SHOW = 100

const props = defineProps<{
  options: Array<CodonTable>
}>()
const model = defineModel<CodonTable | null>()
const filter = ref('')
const collapseDropdown = ref(true)
const optionsToShow = ref(DEFAULT_NUMBER_TO_SHOW)
const summaryRef = ref<HTMLElement | null>(null)
const dropdownStyle = ref({
  top: '0px',
  left: '0px',
  width: '300px',
})

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

function updateDropdownPosition() {
  if (summaryRef.value) {
    const rect = summaryRef.value.getBoundingClientRect()
    dropdownStyle.value = {
      top: `${rect.bottom + 4}px`,
      left: `${rect.left}px`,
      width: `${Math.max(rect.width, 300)}px`,
    }
  }
}

function handleWindowScroll() {
  if (!collapseDropdown.value) {
    updateDropdownPosition()
  }
}

watch(collapseDropdown, value => {
  optionsToShow.value = DEFAULT_NUMBER_TO_SHOW

  if (!value) {
    updateDropdownPosition()
    window.addEventListener('scroll', handleWindowScroll, true)
  } else {
    window.removeEventListener('scroll', handleWindowScroll, true)
  }
})

onUnmounted(() => {
  window.removeEventListener('scroll', handleWindowScroll, true)
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
    <summary ref="summaryRef" class="summary-select" @click.prevent="collapseDropdown = !collapseDropdown">
      <template v-if="model">
        <i>{{ model.organism }}</i> - {{ model.name }}
      </template>

      <template v-else>Select one...</template>
    </summary>

    <ul v-if="!collapseDropdown" class="option-dropdown"
      :style="{ top: dropdownStyle.top, left: dropdownStyle.left, width: dropdownStyle.width }">
      <li>
        <input type="search" placeholder="Filter..." v-model="filter" />
      </li>
      <div class="option-list" @scroll="handleScroll">
        <li>
          <label>
            <input type="radio" v-model="model" :value="null" @change="closeDropdown" />
            (None)
          </label>
        </li>
        <li v-for="option in filteredOptions" :key="option.id">
          <label>
            <input type="radio" v-model="model" :value="option" @change="closeDropdown" />
            <i>{{ option.organism }}</i> - {{ option.name }}
          </label>
        </li>
      </div>
    </ul>
  </details>
</template>

<style scoped>
.dropdown {
  position: relative;
}

.option-dropdown {
  position: fixed;
  z-index: 9999;
  padding: 0;
}

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

@media (max-width: 768px) {
  .option-dropdown {
    min-width: 90vw !important;
    max-width: 90vw !important;
  }

  .option-list {
    width: 100%;
  }
}
</style>
