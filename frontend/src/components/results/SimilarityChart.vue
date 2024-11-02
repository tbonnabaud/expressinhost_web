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
        label: 'Similarity percentages',
        data: props.values,
        borderWidth: 2,
      },
    ],
  }
})

const options: ChartOptions = {
  responsive: true,
  plugins: {
    legend: {
      display: true,
      position: 'top',
    },
  },
}
</script>

<template>
  <ChartWrapper type="bar" :data="data" :options="options" />
</template>
