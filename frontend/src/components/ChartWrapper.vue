<script lang="ts" setup>
import { ref, onMounted, onBeforeUnmount } from 'vue'
import {
  Chart,
  registerables,
  type ChartOptions,
  type ChartData,
  type ChartType,
} from 'chart.js'

// Register Chart.js components
Chart.register(...registerables)

// Define props types
interface ChartProps {
  type: ChartType
  data: ChartData
  options?: ChartOptions
}

// Props for the chart data and options
const props = defineProps<ChartProps>()

// Reference to the canvas element
const chartCanvas = ref<HTMLCanvasElement | null>(null)
let chartInstance: Chart | null = null

// Create the chart when the component is mounted
onMounted(() => {
  if (chartCanvas.value) {
    chartInstance = new Chart(chartCanvas.value, {
      type: props.type,
      data: props.data,
      options: props.options,
    })
  }
})

// Destroy the chart instance when the component is unmounted
onBeforeUnmount(() => {
  if (chartInstance) {
    chartInstance.destroy()
  }
})
</script>

<template>
  <canvas ref="chartCanvas"></canvas>
</template>
