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

function formatDuration(duration: string): string {
  const { hours, minutes, seconds } = parseISODuration(duration)
  console.log(seconds)

  const parts = []

  if (hours > 0) parts.push(`${hours}h`)
  if (minutes > 0) parts.push(`${minutes}m`)
  if (seconds > 0 || parts.length === 0) parts.push(`${seconds.toFixed(2)}s`)

  return parts.join(' ')
}
</script>

<template>
  <div class="stats-container">
    <div class="stat">
      <div class="label">Minimum</div>
      <div class="number">
        {{ formatDuration(durationStats?.min_duration || '0') }}
      </div>
    </div>
    <div class="stat">
      <div class="label">Average</div>
      <div class="number">
        {{ formatDuration(durationStats?.avg_duration || '0') }}
      </div>
    </div>
    <div class="stat">
      <div class="label">Maximum</div>
      <div class="number">
        {{ formatDuration(durationStats?.max_duration || '0') }}
      </div>
    </div>
  </div>
</template>

<style scoped>
.stats-container {
  display: flex;
  gap: 20px;
  justify-content: center;
}

.stat {
  text-align: center;
  background-color: #fcfcfc;
  padding: 20px;
  border-radius: 8px;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

.stat .number {
  font-size: 16px;
  font-weight: bold;
  color: #1d6a54;
}

.stat .label {
  font-size: 16px;
}

@media (max-width: 768px) {
  .stats-container {
    flex-direction: column;
  }
}

@media (prefers-color-scheme: dark) {
  .stat {
    background-color: #13171f;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.5);
  }
}
</style>
