<script setup lang="ts">
import { onMounted, ref } from 'vue'
import { API } from '@/lib/api'
import { parseISODuration } from '@/lib/helpers'
import type { RunInfoDurationStats } from '@/lib/interfaces'

const durationStats = ref(null as RunInfoDurationStats | null)

onMounted(async () => await fetchDurationStats())

async function fetchDurationStats() {
  const [data, error] = await API.runInfos.durationStats()

  if (!error) {
    durationStats.value = data
  }
}
</script>

<template>
  <div id="durationStats">
    <p>
      <strong>Minimum duration: </strong>
      <ins>
        {{ parseISODuration(durationStats?.min_duration || '0') }}
      </ins>
      seconds
    </p>

    <p>
      <strong>Average duration: </strong>
      <ins>
        {{ parseISODuration(durationStats?.avg_duration || '0') }}
      </ins>
      seconds
    </p>

    <p>
      <strong>Maximum duration: </strong>
      <ins>
        {{ parseISODuration(durationStats?.max_duration || '0') }}
      </ins>
      seconds
    </p>
  </div>
</template>

<style scoped>
#durationStats {
  min-height: 30vh;
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
}

#durationStats p {
  font-size: large;
}

#durationStats p ins {
  font-weight: bold;
}
</style>
