<script setup lang="ts">
import { computed } from 'vue'
import ChartWrapper from '../ChartWrapper.vue'
import type { ChartOptions } from 'chart.js'

const props = defineProps<{
  title: string
  inputTitle: string
  inputValues: Array<number>
  outputTitle: string
  outputValues: Array<number>
}>()

const chartWidth = computed(() => {
  const calculatedWidth = props.outputValues.length * 10
  const minWidth = 600 // Minimum for readability
  const maxWidth = 3000 // Prevent excessive width
  return Math.max(minWidth, Math.min(calculatedWidth, maxWidth))
})

const data = computed(() => {
  return {
    labels: [...props.outputValues.keys()],
    datasets: [
      {
        label: props.inputTitle,
        data: props.inputValues,
        pointStyle: false,
        borderColor: '#D55E00',
        backgroundColor: '#D55E00',
        borderWidth: 2,
      },
      {
        label: props.outputTitle,
        data: props.outputValues,
        pointStyle: false,
        borderColor: '#0072B2',
        backgroundColor: '#0072B2',
        borderWidth: 2,
      },
    ],
  }
})

const options: ChartOptions = {
  responsive: true,
  maintainAspectRatio: false,
  plugins: {
    legend: {
      display: props.inputValues.length != 0,
      position: 'left',
    },
    title: {
      display: false,
      text: props.title,
    },
  },
}
</script>

<template>
  <div style="overflow-x: scroll">
    <div class="chart-container" :style="{ width: chartWidth + 'px' }">
      <ChartWrapper
        class="wrapper"
        type="line"
        :data="data"
        :options="options"
      />
    </div>
  </div>
</template>

<style scoped>
.chart-container {
  height: clamp(200px, 20vh, 400px);
  min-width: 100%;
  max-width: 100%;
}

@media (max-width: 768px) {
  .chart-container {
    height: clamp(150px, 25vh, 300px);
  }
}
</style>
