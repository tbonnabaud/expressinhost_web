<script setup lang="ts">
import { store } from '@/lib/store'

const user = store.currentUser
</script>

<template>
  <main class="container">
    <h1 v-if="user">Welcome {{ user.full_name }}</h1>

    <p>
      <strong>ExpressInHost</strong> tunes codon sequences for recombinant
      protein expression in a host organism
    </p>

    <p>⚠ This web application is under development ⚠</p>
  </main>
</template>

<style scoped>
h1,
p {
  margin-top: 3em;
  text-align: center;
}

p {
  font-size: 2.5em;
}
</style>
