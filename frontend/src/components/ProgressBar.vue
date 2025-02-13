<script setup lang="ts">
import { computed } from 'vue'

const props = defineProps<{ value: number; max: number }>()

const percentage = computed(() => {
  if (props.max) {
    return (props.value / props.max) * 100
  } else {
    return 0
  }
})
</script>

<template>
  <div class="progress-container">
    <progress class="progress-bar" :value="value" :max="max"></progress>
    <span class="progress-text">{{ percentage.toFixed(0) }}%</span>
  </div>
</template>

<style scoped>
.progress-container {
  position: relative;
  width: 100%;
}

.progress-bar {
  width: 100%;
  height: 100%;
  background-color: #727a8d;
  border-radius: 0.25rem;
}

.progress-text {
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  color: white;
  font-weight: bold;
  pointer-events: none;
}
</style>
