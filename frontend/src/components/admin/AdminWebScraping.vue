<script setup lang="ts">
import { API } from '@/lib/api'
import { onMounted, onUnmounted, ref, watch } from 'vue'

interface ScrapingState {
  state: string
  done: number
  total: number
}

const scrapingState = ref(null as ScrapingState | null)
const requestInterval = ref(0)

const isLoading = ref(false)

watch(scrapingState, value => {
  if (value && value.done && value.total) {
    isLoading.value = value.done !== value.total
  }
})

onMounted(async () => {
  await getWebScrapingState()
  requestInterval.value = setInterval(async () => {
    await getWebScrapingState()
  }, 2000)
})

onUnmounted(() => clearInterval(requestInterval.value))

async function runWebScraping() {
  isLoading.value = true
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
  <button :aria-busy="isLoading" type="button" @click="runWebScraping">
    Run scraping of Lowe Lab database
  </button>

  <div
    id="scrapingProgressDiv"
    v-if="scrapingState && scrapingState.done && scrapingState.total"
  >
    <span
      >{{
        ((scrapingState.done / scrapingState.total) * 100).toFixed(0)
      }}%</span
    >
    <progress :value="scrapingState.done" :max="scrapingState.total"></progress>
  </div>
</template>

<style scoped>
progress {
  height: 1.5em;
  background-color: #727a8d;
}

#scrapingProgressDiv span {
  position: absolute;
  display: inline-block;
  text-align: center;
  margin-left: 50%;
  font-weight: bold;
  color: #fff;
}

#scrapingProgressDiv {
  display: block;
  position: relative;
  width: 100%;
}
</style>
