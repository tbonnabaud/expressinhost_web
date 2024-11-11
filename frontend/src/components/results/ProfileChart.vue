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
        borderColor: '#a52a2a',
        backgroundColor: '#a52a2a',
      },
      {
        label: 'Output sequence',
        data: props.outputValues,
        pointStyle: false,
        borderColor: '#008000',
        backgroundColor: '#008000',
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
      text: props.title,
    },
  },
}
</script>

<template>
  <div class="chart-wrapper">
    <ChartWrapper type="line" :data="data" :options="options" />
  </div>
</template>

<style scoped>
.chart-wrapper {
  height: 20vh;
}
</style>
