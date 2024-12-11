<script setup lang="ts">
import { API } from '@/lib/api'
import type { User } from '@/lib/interfaces'
import { onMounted, ref } from 'vue'
import ChartRegistrationsByDate from './ChartRegistrationsByDate.vue'

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
  <h3>
    Number of registered users: <ins>{{ userList.length }}</ins>
  </h3>
  <ChartRegistrationsByDate :user-list="userList" />
</template>
