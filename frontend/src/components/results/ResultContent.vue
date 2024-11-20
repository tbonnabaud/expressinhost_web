<script setup lang="ts">
import { computed, onMounted, ref } from 'vue'
import type { TuningOutput, CodonTable } from '@/lib/interfaces'
import { formatToLocaleDateString } from '@/lib/helpers'
import { MODE_LABEL_MAPPING } from '@/lib/referentials'
import { API } from '@/lib/api'
import SequenceComparison from '@/components/results/SequenceComparison.vue'
import DownloadResult from './DownloadResult.vue'
import SimilarityChart from './SimilarityChart.vue'

const props = defineProps<TuningOutput>()

const hostCodonTable = ref(null as CodonTable | null)

const mode = computed(
  () => MODE_LABEL_MAPPING[props.result.mode] || 'Unknown mode',
)
const percentageLabels = computed(() => props.tuned_sequences.map(e => e.name))
const percentageValues = computed(() =>
  props.tuned_sequences.map(e => e.identity_percentage),
)

onMounted(async () => await fetchHostCodonTable())

async function fetchHostCodonTable() {
  const [data, error] = await API.codonTables.get(
    props.result.host_codon_table_id,
  )

  if (!error) {
    hostCodonTable.value = data
  }
}
</script>

<template>
  <br />
  <h2>
    Expression in <i>{{ hostCodonTable?.organism }}</i> -
    {{ hostCodonTable?.name }}
  </h2>

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
