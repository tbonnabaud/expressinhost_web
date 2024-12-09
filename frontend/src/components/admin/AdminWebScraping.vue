<script setup lang="ts">
import { API } from '@/lib/api'

async function runWebScraping() {
  await API.admin.runWebScraping()
}
</script>

<template>
  <button type="button" @click="runWebScraping">
    Run scraping of Lowe Lab database
  </button>
</template>
