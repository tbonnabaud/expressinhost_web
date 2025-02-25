<script setup lang="ts">
import ChartRegistrationsByDate from './ChartRegistrationsByDate.vue'
import ChartRunsByDay from './ChartRunsByDay.vue'
import ChartDurationStats from './ChartDurationStats.vue'
import ChartModeDistribution from './ChartModeDistribution.vue'
</script>

<template>
  <div>
    <h2>Statistics</h2>

    <article>
      <h3>Registered users</h3>
      <ChartRegistrationsByDate />
    </article>

    <article>
      <h3>Runs</h3>
      <ChartRunsByDay />
    </article>

    <div class="stats-row">
      <article>
        <h3>Duration statistics</h3>
        <ChartDurationStats />
      </article>

      <article>
        <h3>Tuning modes</h3>
        <ChartModeDistribution />
      </article>
    </div>
  </div>
</template>

<style scoped>
.stats-row {
  display: flex;
  flex-direction: row;
  column-gap: 1em;
}

.stats-row article {
  flex: 1;
}

@media (max-width: 768px) {
  .stats-row {
    flex-direction: column;
  }
}
</style>
