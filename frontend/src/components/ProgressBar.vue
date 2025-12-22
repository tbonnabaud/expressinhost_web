<script setup lang="ts">
import { computed } from 'vue'

const props = defineProps<{ value: number; max: number; color?: string }>()

const percentage = computed(() => {
  if (props.max) {
    return Math.round((props.value / props.max) * 100)
  } else {
    return 0
  }
})
</script>

<template>
  <div class="progress-container">
    <div class="progress-bar" :style="{ width: `${percentage}%` }"></div>
    <div class="progress-text">{{ percentage }}%</div>
  </div>
</template>

<style scoped>
.progress-container {
  background-color: #727a8d;
  border-radius: 0.25rem;
  position: relative;
  overflow: hidden;
}

.progress-bar {
  height: 100%;
  background-color: v-bind('props.color || "#0172ad"');
  border-radius: 0.25rem;
  transition: width 0.3s ease;
}

.progress-text {
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  color: #fff;
  font-weight: bold;
}
</style>
