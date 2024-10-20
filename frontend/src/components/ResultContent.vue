<script setup lang="ts">
import type { TuningOutput } from '@/lib/interfaces'
import { computed } from 'vue'

const props = defineProps<TuningOutput>()

const modeLabelMapping: Record<string, string> = {
  direct_mapping: 'Direct mapping',
  optimisation_and_conservation_1: 'Optimisation and conservation 1',
  optimisation_and_conservation_2: 'Optimisation and conservation 2',
}

const hostCondonTable = props.result.host_codon_table_name
const mode = computed(
  () => modeLabelMapping[props.result.mode] || 'Unknown mode',
)
</script>

<template>
  <h2>Expression in {{ hostCondonTable }}</h2>

  <p><strong>Mode:</strong> {{ mode }}</p>

  <pre>{{ JSON.stringify(tuned_sequences, null, 2) }}</pre>
</template>
