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
  return RESTRICTION_SITES.filter(
    site =>
      site.enzyme.toLowerCase().includes(lowerCaseFilter) ||
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

  // Maximum of 3 sites
  if (model.value.length < 3) {
    const site = siteToAdd.value

    if (site !== null && !model.value.includes(site)) {
      model.value.push(site)
    }
  } else {
    alert('Maximum number of sites reached.')
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
  <div class="restriction-site-tags">
    <div v-if="model?.length">
      <RestrictionSiteTag
        v-for="site in model"
        :site
        :key="site.enzyme"
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
          <li v-for="option in filteredOptions" :key="option.enzyme">
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
  gap: 0.5em;
  flex-wrap: wrap;
}

.restriction-site-selector > details {
  flex: 1 1 auto;
  min-width: 200px;
}

.restriction-site-selector > button {
  flex: 0 0 auto;
  min-width: 80px;
}

@media (max-width: 768px) {
  .option-list {
    width: calc(100vw - 2rem);
    max-width: 500px;
    position: relative;
    left: 50%;
    transform: translateX(-50%);
  }

  .restriction-site-selector {
    flex-direction: column;
  }

  .restriction-site-selector > details,
  .restriction-site-selector > button {
    width: 100%;
  }
}

/* Very small phone optimization */
@media (max-width: 380px) {
  .option-list {
    width: calc(100vw - 1rem);
  }
}
</style>
