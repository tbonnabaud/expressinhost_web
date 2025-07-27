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

const chartWidth = computed(() => props.outputValues.length * 10)

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
  height: 20vh;
  min-width: 100%;
}
</style>
