<script setup lang="ts">
import { API } from '@/lib/api'
import { Status, useStreamState } from '@/lib/streamedState'
import { onMounted, ref, watch } from 'vue'
import ProgressBar from '../ProgressBar.vue'

const { state: scrapingState, startStream: startScrapingStream } =
  useStreamState(
    '/api/admin/external-db/web-scraping/state',
    'GET',
    localStorage.getItem('accessToken') || undefined,
  )

const isLoading = ref(false)

onMounted(startScrapingStream)

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
</script>

<template>
  <section>
    <h2>Web scraping of codon tables</h2>

    <button
      type="button"
      @click="runWebScraping"
      :aria-busy="isLoading"
      :disabled="isLoading"
    >
      Run scraping of Lowe Lab database
    </button>

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
</style>
