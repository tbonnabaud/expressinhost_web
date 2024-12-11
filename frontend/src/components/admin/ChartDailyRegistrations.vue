<script setup lang="ts">
import type { User } from '@/lib/interfaces'
import type { ChartOptions } from 'chart.js'
import 'chartjs-adapter-date-fns'
import { computed } from 'vue'
import ChartWrapper from '@/components/ChartWrapper.vue'

const props = defineProps<{ userList: Array<User> }>()

const dailyRegistrations = computed(() => {
  return props.userList
    .map(e => {
      const date = new Date(e.creation_date)
      return date.toISOString().split('T')[0]
    })
    .reduce((acc: Record<string, number>, e) => {
      acc[e] = acc[e] ? (acc[e] += 1) : (acc[e] = 1)
      return acc
    }, {})
})

const data = computed(() => {
  return {
    labels: Object.keys(dailyRegistrations.value),
    datasets: [
      {
        label: 'Registrations',
        data: Object.values(dailyRegistrations.value),
      },
    ],
  }
})

const options: ChartOptions = {
  responsive: true,
  maintainAspectRatio: false,
  scales: {
    x: {
      type: 'time',
      time: {
        parser: 'yyyy-MM-dd',
      },
    },
    y: {
      ticks: {
        stepSize: 1,
      },
    },
  },
  plugins: {
    legend: {
      display: true,
      position: 'top',
    },
    title: {
      display: true,
      text: 'Registrations per day',
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
  height: 30vh;
}
</style>
