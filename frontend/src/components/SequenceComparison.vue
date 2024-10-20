<script setup lang="ts">
// import { ref } from 'vue';
import type { TunedSequence } from '@/lib/interfaces'

defineProps<TunedSequence>()
</script>

<template>
  <h3>{{ name }}</h3>

  <p><strong>-> Similarity:</strong> {{ identity_percentage.toFixed(2) }}%</p>

  <div class="sequence-group">
    <div class="sequence" v-html="input"></div>
    <div class="sequence" v-html="output"></div>
  </div>
</template>

<style scoped>
.sequence {
  height: 40px;
  white-space: nowrap;
}

.sequence-group {
  overflow-x: scroll;
  margin: 2em 0;
}

.in {
  background-color: brown;
}

.out {
  background-color: green;
}
</style>
