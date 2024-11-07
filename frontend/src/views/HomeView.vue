<script setup lang="ts">
import { store } from '@/lib/store'

const user = store.currentUser
</script>

<template>
  <main class="container">
    <h1 v-if="user">Welcome {{ user.full_name }}</h1>

    <p>
      ExpressInHost is an application for tuning the codon sequence of an mRNA
      for recombinant protein expression in a host microorganism.
    </p>

    <p>⚠ This web application is under development ⚠</p>
  </main>
</template>

<style scoped>
h1,
p {
  margin-top: 5em;
  text-align: center;
}

p {
  font-size: 1.5em;
}
</style>
