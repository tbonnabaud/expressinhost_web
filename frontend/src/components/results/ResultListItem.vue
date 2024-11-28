<script setup lang="ts">
import type { Result } from '@/lib/interfaces'
import { MODE_LABEL_MAPPING } from '@/lib/referentials'
import { formatToLocaleDateString } from '@/lib/helpers'

defineProps<{ result: Result }>()
</script>

<template>
  <RouterLink :to="'/results/' + result.id" class="contrast">
    <article class="zoomable">
      <header>
        <h5>
          <span class="creation-date">
            {{ formatToLocaleDateString(result.creation_date) }}
          </span>
          &nbsp; <i>{{ result.host_codon_table.organism }}</i> -
          {{ result.host_codon_table.name }}
        </h5>
      </header>

      <p><strong>Mode:</strong> {{ MODE_LABEL_MAPPING[result.mode] }}</p>
      <div class="grid">
        <p>
          <strong>Slow speed thresold:</strong>
          {{ result.slow_speed_threshold }}
        </p>
        <p>
          <strong>Conservation thresold:</strong>
          {{ result.conservation_threshold }}
        </p>
      </div>
    </article>
  </RouterLink>
</template>

<style scoped>
article {
  margin: 0.5em;
  padding-bottom: 0.3em;
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

h5 {
  margin: 0.5em auto;
}

.creation-date {
  border: 1px solid var(--pico-color);
  border-radius: 5px;
  padding: 3px 5px;
}

p {
  text-align: center;
}
</style>
