<script setup lang="ts">
import { computed } from 'vue'
import ChartWrapper from '../ChartWrapper.vue'
import type { ChartOptions } from 'chart.js'

const props = defineProps<{
  title: string
  inputValues: Array<number>
  outputValues: Array<number>
}>()

const chartWidth = computed(() => props.inputValues.length * 10)

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
        borderWidth: 2,
      },
      {
        label: 'Output sequence',
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
  height: 20vh;
  min-width: 100%;
}
</style>
