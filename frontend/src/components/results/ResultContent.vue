<script setup lang="ts">
import { computed, ref } from 'vue'
import { useRouter } from 'vue-router'
import type { TuningOutput } from '@/lib/interfaces'
import { formatToLocaleDateString } from '@/lib/helpers'
import { MODE_LABEL_MAPPING } from '@/lib/referentials'
import { API } from '@/lib/api'
import SequenceComparison from '@/components/results/SequenceComparison.vue'
import DownloadResult from './DownloadResult.vue'
import SimilarityChart from './SimilarityChart.vue'
import BaseModal from '../BaseModal.vue'

const router = useRouter()
const props = defineProps<TuningOutput>()

const openDeleteModal = ref(false)

const mode = computed(
  () => MODE_LABEL_MAPPING[props.result.mode] || 'Unknown mode',
)
const percentageLabels = computed(() => props.tuned_sequences.map(e => e.name))
const percentageValues = computed(() =>
  props.tuned_sequences.map(e => e.identity_percentage),
)

async function deleteResult() {
  if (props.result.id) {
    const [, error] = await API.results.delete(props.result.id)

    if (!error) {
      console.log('Result deleted successfully.')
      // Redirect to results
      router.push('/results')
    }
  }

  openDeleteModal.value = false
}
</script>

<template>
  <BaseModal
    :open="openDeleteModal"
    title="Confirm the deletion"
    @close="openDeleteModal = false"
  >
    <p>Do you really want to delete this result?</p>

    <footer>
      <button class="secondary" @click="openDeleteModal = false">Cancel</button>
      <button class="danger" @click="deleteResult">Delete</button>
    </footer>
  </BaseModal>

  <div class="flex-container">
    <h2>
      {{ result.name }}
    </h2>

    <button
      v-if="result.id"
      id="deleteButton"
      class="danger"
      @click="openDeleteModal = true"
    >
      Delete
    </button>
  </div>

  <p id="createdOn">
    Created on {{ formatToLocaleDateString(result.creation_date) }}.
  </p>

  <hr />

  <h3>Host codon table</h3>
  <hr />
  <p v-if="result.host_codon_table">
    <i>{{ result.host_codon_table.organism }}</i> -
    {{ result.host_codon_table.name }}
  </p>

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
        <td>{{ result.slow_speed_threshold * 100 }}%</td>
        <td>
          {{
            result.conservation_threshold &&
            result.conservation_threshold * 100
          }}%
        </td>
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

<style scoped>
#deleteButton {
  margin-left: auto;
}

#createdOn {
  font-style: italic;
}
</style>
