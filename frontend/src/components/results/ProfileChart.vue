<script setup lang="ts">
import { computed } from 'vue'
import ChartWrapper from '../ChartWrapper.vue'
import type { ChartOptions } from 'chart.js'

const props = defineProps<{
  title: string
  inputValues: Array<number>
  outputValues: Array<number>
}>()

const data = computed(() => {
  return {
    labels: [...props.inputValues.keys()],
    datasets: [
      {
        label: 'Input sequence',
        data: props.inputValues,
        pointStyle: false,
      },
      {
        label: 'Output sequence',
        data: props.outputValues,
        pointStyle: false,
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
    title: {
      display: true,
      text: props.title,
    },
  },
}
</script>

<template>
  <ChartWrapper type="line" :data="data" :options="options" />
</template>
