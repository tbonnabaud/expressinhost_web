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
        borderColor: '#D55E00',
        backgroundColor: '#D55E00',
        borderWidth: 1,
      },
      {
        label: 'Output sequence',
        data: props.outputValues,
        pointStyle: false,
        borderColor: '#0072B2',
        backgroundColor: '#0072B2',
        borderWidth: 1,
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
  <div class="chart-container">
    <ChartWrapper type="line" :data="data" :options="options" />
  </div>
</template>

<style scoped>
.chart-container {
  height: 30vh;
}
</style>
