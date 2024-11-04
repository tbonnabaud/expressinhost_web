<script setup lang="ts">
import type { Result } from '@/lib/interfaces'

defineProps<{ result: Result }>()
</script>

<template>
  <RouterLink :to="'/results/' + result.id" class="contrast">
    <article class="zoomable">
      <header>
        <h4>{{ result.id }}</h4>
      </header>

      <p>
        <strong>Host codon table:</strong> {{ result.host_codon_table_name }}
        <br />
        <strong>Mode:</strong> {{ result.mode }}
      </p>
    </article>
  </RouterLink>
</template>

<style scoped>
article {
  margin: 0.5em;
}

a {
  text-decoration: none;
}

.zoomable {
  transition: transform 0.2s;
}

.zoomable:hover {
  transform: scale(1.03);
}
</style>
