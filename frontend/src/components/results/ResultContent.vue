<script setup lang="ts">
import { computed, ref } from 'vue'
import { useRouter } from 'vue-router'
import { TuningModeName, type TuningOutput } from '@/lib/interfaces'
import { formatToLocaleDateString } from '@/lib/helpers'
import {
  MODE_LABEL_MAPPING,
  FIVE_PRIME_REGION_TUNING_LABEL_MAPPING,
} from '@/lib/referentials'
import { API } from '@/lib/api'
import SequenceComparison from '@/components/results/SequenceComparison.vue'
import SequenceFromProtein from '@/components/results/SequenceFromProtein.vue'
import DownloadResult from './DownloadResult.vue'
import SimilarityChart from './SimilarityChart.vue'
import BaseModal from '../BaseModal.vue'
import RestrictionSiteTag from '../RestrictionSiteTag.vue'

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
        <th v-if="result.slow_speed_threshold" scope="col">
          Slow speed thresold
        </th>
        <th v-if="result.conservation_threshold" scope="col">
          Conservation thresold
        </th>
      </tr>
    </thead>

    <tbody>
      <tr>
        <td>{{ mode }}</td>
        <td v-if="result.slow_speed_threshold">
          {{ result.slow_speed_threshold * 100 }}%
        </td>
        <td v-if="result.conservation_threshold">
          {{ result.conservation_threshold * 100 }}%
        </td>
      </tr>
    </tbody>
  </table>

  <h3>Five prime region tuning</h3>
  <hr />

  <div v-if="result.five_prime_region_tuning">
    <table>
      <thead>
        <tr>
          <th scope="col">Mode</th>
          <th
            scope="col"
            v-if="result.five_prime_region_tuning.mode == 'partial_untuning'"
          >
            Untuned codon number
          </th>
          <th
            scope="col"
            v-else-if="result.five_prime_region_tuning.mode == 'fine_tuning'"
          >
            Codon window size
          </th>
        </tr>
      </thead>

      <tbody>
        <tr>
          <td>
            {{
              FIVE_PRIME_REGION_TUNING_LABEL_MAPPING[
                result.five_prime_region_tuning.mode
              ]
            }}
          </td>
          <td v-if="result.five_prime_region_tuning.mode == 'partial_untuning'">
            {{ result.five_prime_region_tuning.untuned_codon_number }}
          </td>
          <td v-else-if="result.five_prime_region_tuning.mode == 'fine_tuning'">
            {{ result.five_prime_region_tuning.codon_window_size }}
          </td>
        </tr>
      </tbody>
    </table>

    <article
      id="utrSequence"
      v-if="result.five_prime_region_tuning.mode == 'fine_tuning'"
    >
      <details name="utrAccordion">
        <summary>UTR sequence</summary>
        {{ result.five_prime_region_tuning.utr }}
      </details>
    </article>
  </div>

  <div v-else>None</div>

  <h3>Restriction enzyme recognition sites to avoid</h3>
  <hr />

  <div v-if="result.restriction_sites">
    <RestrictionSiteTag
      v-for="site in result.restriction_sites"
      :site
      :key="site.enzyme"
      :deletable="false"
    />
  </div>

  <div v-else>None</div>

  <h3>Sequence profiles</h3>
  <hr />

  <div v-if="result.mode == TuningModeName.PROTEIN_STRUCTURE_ANALYSIS">
    <SequenceFromProtein
      v-for="(item, index) in tuned_sequences"
      :tuned-sequence="item"
      :key="index"
    />
  </div>

  <div v-else class="sequence-comparisons">
    <SequenceComparison
      v-for="(item, index) in tuned_sequences"
      :tuned-sequence="item"
      :native-codon-table-id="result.sequences_native_codon_tables[item.name]"
      :open="index == 0"
      :key="index"
    />
  </div>

  <SimilarityChart
    v-if="result.mode != TuningModeName.PROTEIN_STRUCTURE_ANALYSIS"
    :labels="percentageLabels"
    :values="percentageValues"
  />

  <DownloadResult :result="result" :tuned_sequences="tuned_sequences" />
</template>

<style scoped>
table th,
table td {
  text-align: center;
}

#utrSequence {
  padding-bottom: 0.3em;
}

#deleteButton {
  margin-left: auto;
}

#createdOn {
  font-style: italic;
}
</style>
