<script setup lang="ts">
import { API } from '@/lib/api'
import { downloadFile } from '@/lib/helpers'

async function downloadLogFile() {
  const [data, error] = await API.admin.getLogFileContent()

  if (!error) {
    downloadFile(data, 'expressinhost.log')
  }
}
</script>

<template>
  <button @click="downloadLogFile">Download log file</button>
</template>
