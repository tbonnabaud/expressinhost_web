<script setup lang="ts">
import TuningForm from '@/components/TuningForm.vue'
</script>

<template>
  <main class="container">
    <h1>Sequence tuning</h1>

    <hr />

    <TuningForm />
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
