<script setup lang="ts">
import { computed } from 'vue'
import type { TuningOutput } from '@/lib/interfaces'
import SequenceComparison from '@/components/SequenceComparison.vue'

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
  <br />
  <h2>Expression in {{ hostCondonTable }}</h2>

  <p><strong>Mode:</strong> {{ mode }}</p>

  <h3>Sequence comparisons</h3>
  <hr />

  <SequenceComparison
    v-for="(item, index) in tuned_sequences"
    :tuned-sequence="item"
    :open="index == 0"
    :key="index"
  />
</template>
