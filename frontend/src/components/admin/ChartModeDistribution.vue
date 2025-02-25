<script setup lang="ts">
import { computed, onMounted, ref } from 'vue'
import { API } from '@/lib/api'
import { MODE_LABEL_MAPPING } from '@/lib/referentials'
import type { ChartOptions } from 'chart.js'
import ChartWrapper from '@/components/ChartWrapper.vue'

const modeDistribution = ref({} as Record<string, number>)

const data = computed(() => {
  return {
    labels: Object.keys(modeDistribution.value).map(e => MODE_LABEL_MAPPING[e]),
    datasets: [
      {
        label: 'Modes',
        data: Object.values(modeDistribution.value),
      },
    ],
  }
})

const options: ChartOptions = {
  responsive: true,
  maintainAspectRatio: false,
  plugins: {
    legend: {
      display: true,
      position: 'top',
    },
    title: {
      display: true,
      text: 'Mode distribution',
    },
  },
}

onMounted(async () => await fetchDurationStats())

async function fetchDurationStats() {
  const [data, error] = await API.runInfos.modeDistribution()

  if (!error) {
    modeDistribution.value = data
  }
}
</script>

<template>
  <strong>
    Total: <ins>{{ Object.keys(MODE_LABEL_MAPPING).length }}</ins>
  </strong>

  <div class="chart-wrapper">
    <ChartWrapper type="doughnut" :data="data" :options="options" />
  </div>
</template>

<style scoped>
.chart-wrapper {
  height: 30vh;
}
</style>
