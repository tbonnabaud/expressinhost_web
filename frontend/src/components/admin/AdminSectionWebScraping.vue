<script setup lang="ts">
import { API } from '@/lib/api'
import { Status, useStreamState } from '@/lib/streamedState'
import { ref, watch } from 'vue'

const { state: scrapingState } = useStreamState(
  '/api/admin/external-db/web-scraping/state',
  localStorage.getItem('accessToken') || undefined,
)

const isLoading = ref(false)

watch(scrapingState, value => {
  if (value) {
    isLoading.value = value.status == Status.RUNNING
  }
})

async function runWebScraping() {
  isLoading.value = true
  await API.admin.runWebScraping()
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
    <p>{{ scrapingState?.message }}</p>

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
