<script setup lang="ts">
import { API } from '@/lib/api'
import type { User } from '@/lib/interfaces'
import { onMounted, ref } from 'vue'
import ChartDailyRegistrations from './ChartDailyRegistrations.vue'

const userList = ref([] as Array<User>)

onMounted(async () => await fetchUserList())

async function fetchUserList() {
  const [data, error] = await API.users.list()

  if (!error && data !== null) {
    userList.value = data
  }
}
</script>

<template>
  <h3>Number of registered users: {{ userList.length }}</h3>
  <ChartDailyRegistrations :user-list="userList" />
</template>
