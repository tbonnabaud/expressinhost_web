<script setup lang="ts">
import { onMounted, ref, watch } from 'vue'
import ProgressBar from '../ProgressBar.vue'
import { API } from '@/lib/api'
import { Status, useStreamState } from '@/lib/streamedState'
import type { LastWebScraping } from '@/lib/interfaces'
import { formatToLocaleDateString } from '@/lib/helpers'

const isLoading = ref(false)
const lastWebScraping = ref<LastWebScraping | null>(null)
const { state: scrapingState, startStream: startScrapingStream } =
  useStreamState(localStorage.getItem('accessToken') || undefined)

onMounted(fetchLastRelease)
onMounted(async () => {
  const jobId = localStorage.getItem('webScrapingJobId')

  if (jobId) {
    await startScrapingStream(
      `/api/admin/external-db/web-scraping/state/${jobId}`,
    )
  }
})

watch(scrapingState, state => {
  if (state) {
    isLoading.value = [Status.STARTED, Status.QUEUED].includes(state.status)

    if (state.status == Status.NOT_FOUND) {
      localStorage.removeItem('webScrapingJobId')
    }
  }
})

async function runWebScraping() {
  isLoading.value = true
  const [jobId, error] = await API.admin.runWebScraping()

  if (!error) {
    await startScrapingStream(
      `/api/admin/external-db/web-scraping/state/${jobId}`,
    )
    localStorage.setItem('webScrapingJobId', jobId)
  } else {
    console.error('No job ID.')
  }
}

async function fetchLastRelease() {
  const [data, error] = await API.admin.getWebScrapingLastRelease()

  if (!error) {
    lastWebScraping.value = data
  }
}
</script>

<template>
  <div>
    <h2>Web scraping of codon tables</h2>

    <div class="grid">
      <button
        id="scrapingButton"
        type="button"
        @click="runWebScraping"
        :aria-busy="isLoading"
        :disabled="isLoading"
      >
        Run scraping of Lowe Lab database
      </button>

      <p v-if="lastWebScraping" id="scrapingInfos">
        <strong>Last release in database:</strong>
        {{ lastWebScraping.release }}, scraped on
        {{ formatToLocaleDateString(lastWebScraping.scraping_date) }}
      </p>
    </div>

    <p>{{ scrapingState?.status }}</p>

    <ProgressBar
      v-if="scrapingState"
      id="progressBar"
      :value="scrapingState.step || 0"
      :max="scrapingState.total || 0"
    />

    <p id="stateMessage">{{ scrapingState?.message }}</p>
  </div>
</template>

<style scoped>
#progressBar {
  height: 3em;
}

#stateMessage {
  margin-top: 1em;
}

#scrapingButton {
  height: 3em;
}

#scrapingInfos {
  text-align: center;
  border: 1px dashed grey;
  border-radius: 0.25rem;
  padding: 0 1em;
}
</style>
