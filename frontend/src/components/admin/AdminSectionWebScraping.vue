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
  useStreamState(
    '/api/admin/external-db/web-scraping/state',
    'GET',
    localStorage.getItem('accessToken') || undefined,
  )

onMounted(startScrapingStream)
onMounted(fetchLastRelease)

watch(scrapingState, value => {
  if (value) {
    isLoading.value = value.status == Status.RUNNING
  }
})

async function runWebScraping() {
  isLoading.value = true
  await API.admin.runWebScraping()
  await startScrapingStream()
}

async function fetchLastRelease() {
  const [data, error] = await API.admin.getWebScrapingLastRelease()

  if (!error) {
    lastWebScraping.value = data
  }
}
</script>

<template>
  <section>
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
        Last release: {{ lastWebScraping.release }}, scraped on
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
  </section>
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
  line-height: 3em;
  border: 1px dashed grey;
  border-radius: 0.25rem;
}
</style>
