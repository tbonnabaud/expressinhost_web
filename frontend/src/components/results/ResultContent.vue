<script setup lang="ts">
import { computed } from 'vue'
import type { TuningOutput } from '@/lib/interfaces'
import { formatToLocaleDateString } from '@/lib/helpers'
import { MODE_LABEL_MAPPING } from '@/lib/referentials'
import SequenceComparison from '@/components/results/SequenceComparison.vue'
import DownloadResult from './DownloadResult.vue'
import SimilarityChart from './SimilarityChart.vue'

const props = defineProps<TuningOutput>()

const mode = computed(
  () => MODE_LABEL_MAPPING[props.result.mode] || 'Unknown mode',
)
const percentageLabels = computed(() => props.tuned_sequences.map(e => e.name))
const percentageValues = computed(() =>
  props.tuned_sequences.map(e => e.identity_percentage),
)
</script>

<template>
  <br />
  <h2>Expression in {{ result.host_codon_table_name }}</h2>

  <i>Created on {{ formatToLocaleDateString(result.creation_date) }}.</i>

  <hr />

  <h3>Parameters</h3>
  <hr />

  <table>
    <thead>
      <tr>
        <th scope="col">Mode</th>
        <th scope="col">Slow speed thresold</th>
        <th scope="col">Conservation thresold</th>
      </tr>
    </thead>

    <tbody>
      <tr>
        <td>{{ mode }}</td>
        <td>{{ result.slow_speed_threshold }}</td>
        <td>{{ result.conservation_threshold }}</td>
      </tr>
    </tbody>
  </table>

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
