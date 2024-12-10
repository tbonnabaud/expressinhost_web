<script setup lang="ts">
import { API } from '@/lib/api'
import { onMounted, onUnmounted, ref } from 'vue'

interface ScrapingState {
  state: string
  done: number
  total: number
}

const scrapingState = ref(null as ScrapingState | null)
const requestInterval = ref(0)

onMounted(() => {
  requestInterval.value = setInterval(async () => {
    // console.log('test')
    await getWebScrapingState()
  }, 2000)
})

onUnmounted(() => clearInterval(requestInterval.value))

async function runWebScraping() {
  await API.admin.runWebScraping()
}

async function getWebScrapingState() {
  const [data, error] = await API.admin.getWebScrapingState()

  if (!error && data) {
    scrapingState.value = data
  }
}
</script>

<template>
  <button type="button" @click="runWebScraping">
    Run scraping of Lowe Lab database
  </button>

  <div v-if="scrapingState && scrapingState.done && scrapingState.total">
    <progress :value="scrapingState.done" :max="scrapingState.total" />
  </div>
</template>
