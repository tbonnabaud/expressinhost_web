<script setup lang="ts">
import type { User } from '@/lib/interfaces'
import { API } from '@/lib/api'
import type { ChartOptions } from 'chart.js'
import 'chartjs-adapter-date-fns'
import { computed, onMounted, ref } from 'vue'
import ChartWrapper from '@/components/ChartWrapper.vue'

const userList = ref([] as Array<User>)

const registrationsByMonth = computed(() => {
  return userList.value
    .map(e => getYearMonthFromDateString(e.creation_date))
    .reduce((acc: Record<string, number>, e) => {
      acc[e] = acc[e] ? acc[e] + 1 : 1
      return acc
    }, {})
})

const data = computed(() => {
  return {
    labels: Object.keys(registrationsByMonth.value),
    datasets: [
      {
        label: 'Registrations',
        data: Object.values(registrationsByMonth.value),
        borderWidth: 2,
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
        parser: 'yyyy-MM',
        unit: 'month',
        tooltipFormat: 'MMMM yyyy',
      },
      max: getYearMonthFromDateString(new Date().toISOString()),
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
      text: 'Registrations per month',
    },
  },
}

onMounted(async () => await fetchUserList())

function getYearMonthFromDateString(dateISOString: string): string {
  const [year, month] = dateISOString.split('-')
  return `${year}-${month}`
}

async function fetchUserList() {
  const [data, error] = await API.users.list()

  if (!error && data !== null) {
    userList.value = data
  }
}
</script>

<template>
  <strong>
    Total: <ins>{{ userList.length }}</ins>
  </strong>

  <div class="chart-container">
    <ChartWrapper type="bar" :data="data" :options="options" />
  </div>
</template>

<style scoped>
.chart-container {
  height: 30vh;
}
</style>
