<script setup lang="ts">
import { ref } from 'vue'
import ChartWrapper from '../ChartWrapper.vue'
import type { ChartOptions, ChartType } from 'chart.js'

const props = defineProps<{
  labels: Array<string>
  values: Array<number>
}>()

const chartType = ref('bar' as ChartType)
const chartData = ref({
  labels: props.labels,
  datasets: [
    {
      label: 'Similarity percentages',
      data: props.values,
      borderWidth: 2,
    },
  ],
})

const chartOptions = ref({
  responsive: true,
  plugins: {
    legend: {
      display: true,
      position: 'top',
    },
  },
} as ChartOptions)
</script>

<template>
  <ChartWrapper :type="chartType" :data="chartData" :options="chartOptions" />
</template>
