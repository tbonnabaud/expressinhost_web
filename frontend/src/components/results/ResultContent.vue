<script setup lang="ts">
import { computed } from 'vue'
import type { TuningOutput } from '@/lib/interfaces'
import SequenceComparison from '@/components/results/SequenceComparison.vue'
import DownloadResult from './DownloadResult.vue'
import SimilarityChart from './SimilarityChart.vue'

const props = defineProps<TuningOutput>()

const modeLabelMapping: Record<string, string> = {
  direct_mapping: 'Direct mapping',
  optimisation_and_conservation_1: 'Optimisation and conservation 1',
  optimisation_and_conservation_2: 'Optimisation and conservation 2',
}

const mode = computed(
  () => modeLabelMapping[props.result.mode] || 'Unknown mode',
)
const percentageLabels = computed(() => props.tuned_sequences.map(e => e.name))
const percentageValues = computed(() =>
  props.tuned_sequences.map(e => e.identity_percentage),
)
</script>

<template>
  <br />
  <h2>Expression in {{ result.host_codon_table_name }}</h2>

  <p><strong>Mode:</strong> {{ mode }}</p>

  <h3>Sequence comparisons</h3>
  <hr />

  <SequenceComparison
    v-for="(item, index) in tuned_sequences"
    :tuned-sequence="item"
    :open="index == 0"
    :key="index"
  />

  <SimilarityChart :labels="percentageLabels" :values="percentageValues" />

  <DownloadResult :result="result" :tuned_sequences="tuned_sequences" />
</template>
