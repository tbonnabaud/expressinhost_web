<script setup lang="ts">
import type { RestrictionSite } from '@/lib/interfaces'
defineProps<{ site: RestrictionSite; deletable: boolean }>()
defineEmits(['delete'])
</script>

<template>
  <span class="restriction-site-tag">
    <i>{{ site.enzyme }}</i> - {{ site.sequence }}
    <button
      v-if="deletable"
      type="button"
      class="close-button"
      @click="$emit('delete', site.sequence)"
    >
      X
    </button>
  </span>
</template>

<style scoped>
.restriction-site-tag {
  display: inline-block;
  margin: 5px 5px;
  background-color: #3ca33c;
  border-radius: 5px;
  padding: 10px 15px;
  font-weight: bold;
  color: white;
}

.close-button {
  background: none;
  padding: 0 5px;
  margin: 0 0 0 1em;
  border: 1px white solid;
}

@media (prefers-color-scheme: dark) {
  .restriction-site-tag {
    background-color: #297029;
  }
}
</style>
