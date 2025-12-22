<script setup lang="ts">
import ResultList from '@/components/results/ResultList.vue'
import RequiredAuth from '@/components/RequiredAuth.vue'
import { store } from '@/lib/store'

const user = store.currentUser
</script>

<template>
  <main class="container">
    <template v-if="user">
      <h1>History</h1>
      <hr />

      <ResultList />
    </template>

    <RequiredAuth v-else />
  </main>
</template>
