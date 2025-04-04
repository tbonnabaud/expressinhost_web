<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import { RESTRICTION_SITES } from '@/lib/referentials'
import type { RestrictionSite } from '@/lib/interfaces'
import RestrictionSiteTag from '../RestrictionSiteTag.vue'

const DEFAULT_NUMBER_TO_SHOW = 100

const model = defineModel<RestrictionSite[]>()
const siteToAdd = ref<RestrictionSite | null>(null)
const filter = ref('')
const collapseDropdown = ref(true)
const optionsToShow = ref(DEFAULT_NUMBER_TO_SHOW)

const filteredOptions = computed(() => {
  const lowerCaseFilter = filter.value.toLowerCase()
  return RESTRICTION_SITES.filter(site =>
    site.sequence.toLowerCase().includes(lowerCaseFilter),
  ).slice(0, optionsToShow.value)
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

function addSiteSequence() {
  if (model.value == null) {
    model.value = []
  }

  const site = siteToAdd.value

  if (site !== null && !model.value.includes(site)) {
    model.value.push(site)
  }
}

function removeSite(siteSequence: string) {
  if (model.value) {
    model.value = model.value.filter(e => e.sequence !== siteSequence)
  }
}

function closeDropdown() {
  collapseDropdown.value = true
  filter.value = ''
}
</script>

<template>
  <div>
    <div v-if="model?.length" class="restriction-site-tags">
      <RestrictionSiteTag
        v-for="site in model"
        :site
        :key="site.sequence"
        :deletable="true"
        @delete="removeSite"
      />
    </div>

    <p v-else>No site selected.</p>
  </div>

  <div class="restriction-site-selector">
    <details class="dropdown" :open="!collapseDropdown">
      <summary
        class="summary-select"
        @click.prevent="collapseDropdown = !collapseDropdown"
      >
        <template v-if="siteToAdd">
          <i>{{ siteToAdd.enzyme }}</i> - {{ siteToAdd.sequence }}
        </template>

        <template v-else>Select one...</template>
      </summary>

      <ul v-if="!collapseDropdown" class="option-dropdown">
        <li>
          <input type="search" placeholder="Filter..." v-model="filter" />
        </li>
        <div class="option-list" @scroll="handleScroll">
          <li>
            <label>
              <input
                type="radio"
                v-model="siteToAdd"
                :value="null"
                @change="closeDropdown"
              />
              (None)
            </label>
          </li>
          <li v-for="option in filteredOptions" :key="option.sequence">
            <label>
              <input
                type="radio"
                v-model="siteToAdd"
                :value="option"
                @change="closeDropdown"
              />
              <i>{{ option.enzyme }}</i> - {{ option.sequence }}
            </label>
          </li>
        </div>
      </ul>
    </details>

    <button type="button" @click="addSiteSequence">Add</button>
  </div>
</template>

<style scoped>
.restriction-site-tags {
  text-align: center;
  margin-bottom: 1em;
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

.restriction-site-selector {
  display: flex;
}

.restriction-site-selector > details {
  flex: 1;
}

.restriction-site-selector > button {
  margin-left: 0.5em;
}

@media (max-width: 768px) {
  .option-list {
    width: 90vw;
  }
}
</style>
