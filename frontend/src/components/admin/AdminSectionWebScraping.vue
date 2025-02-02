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
  <section>
    <h2>Web scraping of codon tables</h2>

    <button :aria-busy="isLoading" type="button" @click="runWebScraping">
      Run scraping of Lowe Lab database
    </button>

    <div id="scrapingProgressWrapper">
      <div
        id="scrapingProgress"
        v-if="scrapingState && scrapingState.done && scrapingState.total"
      >
        <span>
          {{ ((scrapingState.done / scrapingState.total) * 100).toFixed(0) }}%
        </span>
        <progress
          id="progressBar"
          :value="scrapingState.done"
          :max="scrapingState.total"
        ></progress>
      </div>
    </div>
  </section>
</template>

<style scoped>
#progressBar {
  background-color: #727a8d;
  height: 1.5em;
}

#scrapingProgress span {
  position: absolute;
  display: inline-block;
  text-align: center;
  margin-left: 50%;
  font-weight: bold;
  color: #fff;
}

#scrapingProgress {
  display: block;
  position: relative;
  width: 100%;
}

#scrapingProgressWrapper {
  height: 1.5em;
}
</style>
