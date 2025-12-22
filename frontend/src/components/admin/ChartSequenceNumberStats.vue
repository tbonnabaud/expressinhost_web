<script setup lang="ts">
import { onMounted, ref } from 'vue'
import { API } from '@/lib/api'
import type { RunInfoSequenceNumberStats } from '@/lib/interfaces'

const sequenceNumberStats = ref(null as RunInfoSequenceNumberStats | null)

onMounted(async () => await fetchDurationStats())

async function fetchDurationStats() {
  const [data, error] = await API.runInfos.sequenceNumberStats()

  if (!error) {
    sequenceNumberStats.value = data
  }
}
</script>

<template>
  <div class="stats-container">
    <div class="stat">
      <div class="label">Minimum</div>
      <div class="number">
        {{ sequenceNumberStats?.min_sequence_number || 0 }}
      </div>
    </div>
    <div class="stat">
      <div class="label">Average</div>
      <div class="number">
        {{ Math.round(sequenceNumberStats?.avg_sequence_number || 0) }}
      </div>
    </div>
    <div class="stat">
      <div class="label">Maximum</div>
      <div class="number">
        {{ sequenceNumberStats?.max_sequence_number || 0 }}
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
  font-size: 24px;
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
