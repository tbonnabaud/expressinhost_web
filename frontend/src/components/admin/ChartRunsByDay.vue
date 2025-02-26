<script setup lang="ts">
import { API } from '@/lib/api'
import type { ChartOptions } from 'chart.js'
import 'chartjs-adapter-date-fns'
import { computed, onMounted, ref } from 'vue'
import ChartWrapper from '@/components/ChartWrapper.vue'

const runsPerDay = ref({} as Record<string, number>)

const totalRuns = computed(() =>
  Object.values(runsPerDay.value).reduce((x, y) => x + y, 0),
)

const data = computed(() => {
  return {
    labels: Object.keys(runsPerDay.value),
    datasets: [
      {
        label: 'Runs',
        data: Object.values(runsPerDay.value),
      },
    ],
  }
})

const options: ChartOptions = {
  responsive: true,
  maintainAspectRatio: false,
  scales: {
    x: {
      type: 'time',
      time: {
        parser: 'yyyy-MM-dd',
        unit: 'month',
        tooltipFormat: 'dd MMMM yyyy',
      },
      max: formatCurrentDay(),
    },
    y: {
      ticks: {
        stepSize: 1,
        precision: 0,
      },
    },
  },
  plugins: {
    legend: {
      display: true,
      position: 'top',
    },
    title: {
      display: true,
      text: 'Runs per day',
    },
  },
}

onMounted(async () => await fetchCountPerDay())

function formatCurrentDay() {
  const currentDate = new Date()
  return currentDate.toISOString().split('T')[0]
}

async function fetchCountPerDay() {
  const [data, error] = await API.runInfos.countPerDay()

  if (!error && data !== null) {
    runsPerDay.value = data
  }
}
</script>

<template>
  <strong>
    Total: <ins>{{ totalRuns }}</ins>
  </strong>

  <div class="chart-container">
    <ChartWrapper type="line" :data="data" :options="options" />
  </div>
</template>

<style scoped>
.chart-container {
  height: 30vh;
}
</style>
