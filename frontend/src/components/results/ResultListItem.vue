<script setup lang="ts">
import type { Result } from '@/lib/interfaces'

defineProps<{ result: Result }>()
</script>

<template>
  <article class="zoomable">
    <header>
      <RouterLink :to="'/results/' + result.id" class="contrast">
        <h4>{{ result.id }}</h4>
      </RouterLink>
    </header>

    <p>{{ result.host_codon_table_name }}</p>
  </article>
</template>
