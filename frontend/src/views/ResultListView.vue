<script setup lang="ts">
import ResultList from '@/components/results/ResultList.vue'
</script>

<template>
  <main class="container">
    <h1>History</h1>

    <ResultList />
  </main>
</template>

<style scoped>
@media (min-width: 1024px) {
  .container {
    padding-left: 15em;
    padding-right: 15em;
  }
}
</style>
