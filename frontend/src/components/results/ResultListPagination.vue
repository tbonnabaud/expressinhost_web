<script setup lang="ts">
defineProps<{
  total: number
}>()

const current = defineModel()
</script>

<template>
  <div class="paginator">
    <button
      v-for="pageNumber in total"
      :key="pageNumber"
      :class="{
        outline: true,
        secondary: true,
        contrast: current == pageNumber,
      }"
      @click="current = pageNumber"
    >
      {{ pageNumber }}
    </button>
  </div>
</template>

<style scoped>
.paginator {
  margin: 1em auto 1em;
  text-align: center;
}

.paginator > button {
  margin: 5px;
}
</style>
