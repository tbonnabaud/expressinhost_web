<script setup lang="ts">
import { computed } from 'vue'
import ChartWrapper from '../ChartWrapper.vue'
import type { ChartOptions } from 'chart.js'

const props = defineProps<{
  labels: Array<string>
  values: Array<number>
}>()

const data = computed(() => {
  return {
    labels: props.labels,
    datasets: [
      {
        label: 'Similarity percentage',
        data: props.values,
        borderWidth: 2,
      },
    ],
  }
})

const options: ChartOptions = {
  responsive: true,
  scales: {
    y: {
      min: 0,
      max: 100,
    },
  },
  plugins: {
    title: {
      display: true,
      text: 'Similarity percentages between the input and the tuned sequences',
    },
    legend: {
      display: true,
      position: 'top',
    },
    tooltip: {
      callbacks: {
        label: context => {
          let label = context.dataset.label || ''

          if (label) {
            label += ': '
          }

          if (context.parsed.y !== null) {
            label += `${context.parsed.y.toFixed()}%`
          }

          return label
        },
      },
    },
  },
}
</script>

<template>
  <ChartWrapper type="bar" :data="data" :options="options" />
</template>
