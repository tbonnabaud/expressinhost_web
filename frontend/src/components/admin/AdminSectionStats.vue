<script setup lang="ts">
import ChartRegistrationsByDate from './ChartRegistrationsByDate.vue'
import ChartRunsByDay from './ChartRunsByDay.vue'
</script>

<template>
  <div>
    <h2>Statistics</h2>

    <article>
      <h3>Registered users</h3>
      <ChartRegistrationsByDate />
    </article>

    <article>
      <h3>Runs per day</h3>
      <ChartRunsByDay />
    </article>
  </div>
</template>
