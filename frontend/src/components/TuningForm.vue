<script setup lang="ts">
import { ref, onMounted, reactive, watch, computed } from 'vue'
import {
  TuningModeName,
  type CodonTable,
  type RunTrainingForm,
} from '@/lib/interfaces'
import { API } from '@/lib/api'
import { store } from '@/lib/store'
import { Status, useStreamState } from '@/lib/streamedState'
import CodonTableSearchSelect from '@/components/codon-tables/CodonTableSearchSelect.vue'
import ToolTip from '@/components/ToolTip.vue'
import AlertError from '@/components/AlertError.vue'
import ProgressBar from '@/components/ProgressBar.vue'
import TuningModeSelector from '@/components/tuning-form/TuningModeSelector.vue'
import SlowSpeedThresholdSelector from '@/components/tuning-form/SlowSpeedThresholdSelector.vue'
import ConservationThresholdSelector from '@/components/tuning-form/ConservationThresholdSelector.vue'
import FivePrimeRegionTuning from '@/components/tuning-form/FivePrimeRegionTuning.vue'
import FastaInput from '@/components/tuning-form/FastaInput.vue'
import ClustalInput from '@/components/tuning-form/ClustalInput.vue'
import PdbInput from '@/components/tuning-form/PdbInput.vue'
import RestrictionSiteSelector from '@/components/tuning-form/RestrictionSiteSelector.vue'
import RsaThresholdSelector from '@/components/tuning-form/RsaThresholdSelector.vue'

const emit = defineEmits(['submit'])

const baseForm: RunTrainingForm = reactive({
  name: '',
  nucleotide_file_content: '',
  pdb_file_content: '',
  clustal_file_content: '',
  host_codon_table_id: '',
  sequences_native_codon_tables: {},
  mode: TuningModeName.DIRECT_MAPPING,
  slow_speed_threshold: 0.5,
  conservation_threshold: null,
  rsa_threshold: null,
  five_prime_region_tuning: null,
  restriction_sites: [],
  send_email: false,
})

const baseFormErrors = reactive({
  nucleotide_file_content: [] as string[],
  clustal_file_content: [] as string[],
})

const selectedHostCodonTable = ref(null as CodonTable | null)
const selectedSequencesNativeCodonTables = ref(
  {} as Record<string, CodonTable | null>,
)

const codonTableList = ref([] as Array<CodonTable>)
const tuningLoading = ref(false)
const currentJobId = ref(null as string | null)

const {
  state: tuningState,
  startStream: startTuningStream,
  stopStream: stopTuningStream,
} = useStreamState(localStorage.getItem('accessToken') || undefined)

const clustalIsRequired = computed(() =>
  [
    'optimisation_and_conservation_1',
    'optimisation_and_conservation_2',
  ].includes(baseForm.mode),
)

const isGuest = computed(() => store.currentUser.value === null)

onMounted(async () => await fetchCodonTables())
onMounted(async () => {
  currentJobId.value = localStorage.getItem('tuningJobId')

  if (currentJobId.value) {
    await startTuningStream(`/api/tuning/state/${currentJobId.value}`)
  }
})

watch(
  () => baseForm.nucleotide_file_content,
  content => {
    // Reset object
    selectedSequencesNativeCodonTables.value = {}

    if (content && baseFormErrors.nucleotide_file_content.length == 0) {
      const sequenceNames = parseFastaSequenceNames(content)

      for (const seq of sequenceNames) {
        const table = findCorrespondingTable(seq)
        selectedSequencesNativeCodonTables.value[seq] = table
      }
    }
  },
)

watch(tuningState, state => {
  if (state) {
    // After the job has been ended remove tuningJobId from local storage
    if (![Status.STARTED, Status.QUEUED].includes(state.status)) {
      localStorage.removeItem('tuningJobId')
      currentJobId.value = null
    }

    if (state.status == Status.FINISHED) {
      emit('submit', state.result)
    }
  }
})

watch(
  () => baseForm.mode,
  value => {
    // Set default values
    baseForm.conservation_threshold =
      value === TuningModeName.OPTIMISATION_AND_CONSERVATION_2 ? 0.75 : null
    baseForm.rsa_threshold =
      value === TuningModeName.PROTEIN_STRUCTURE_ANALYSIS ? 0.25 : null
  },
)

async function fetchCodonTables() {
  const [data, error] = await API.codonTables.list()

  if (!error) {
    codonTableList.value = data
  }
}

function parseFastaSequenceNames(content: string) {
  // Match the group after ">" symbol
  const fastaSeqIdentifierRegex = /^\>(\S+ ?.*)/gm

  return Array.from(content.matchAll(fastaSeqIdentifierRegex), m => m[1])
}

/**
 * Try to guess the first corresponding codon table name using organism name.
 * @param {string} sequenceName - Name of the sequence
 * @returns First corresponding codon table name or empty string if not found.
 */
function findCorrespondingTable(sequenceName: string) {
  const lowerCaseSequenceName = sequenceName
    .toLowerCase()
    .replace(/[\s_-]+/g, ' ')

  for (const tableName of codonTableList.value) {
    if (lowerCaseSequenceName.includes(tableName.organism.toLowerCase())) {
      return tableName
    }
  }

  return null
}

function formIsValid() {
  for (const [key, value] of Object.entries(
    selectedSequencesNativeCodonTables.value,
  )) {
    if (value === null) {
      alert(`Missing table for sequence ${key}`)
      return false
    }
  }

  if (selectedHostCodonTable.value === null) {
    alert('Missing table for host organism')
    return false
  }

  return true
}

async function runTuning() {
  tuningLoading.value = true

  if (formIsValid()) {
    const form = { ...baseForm }

    if (selectedHostCodonTable.value) {
      form.host_codon_table_id = selectedHostCodonTable.value.id
    }

    Object.entries(selectedSequencesNativeCodonTables.value).forEach(
      ([key, value]) => {
        if (value) {
          form.sequences_native_codon_tables[key] = value.id
        }
      },
    )

    const [jobId, error] = await API.runTraining(form)

    if (!error) {
      localStorage.setItem('tuningJobId', jobId)
      currentJobId.value = jobId
      await startTuningStream(`/api/tuning/state/${jobId}`)
    }
  } else {
    tuningLoading.value = false
  }
}

async function cancelTuning() {
  if (currentJobId.value) {
    const [data, error] = await API.cancelTuning(currentJobId.value)

    if (!error) {
      console.log(data)
    } else {
      console.error(error)
    }
  } else {
    console.log('No job ID.')
  }
  stopTuningStream()
  localStorage.removeItem('tuningJobId')
  currentJobId.value = null
  tuningState.value = null
  tuningLoading.value = false
}
</script>

<template>
  <form id="runTuningForm" @submit.prevent="runTuning">
    <section>
      <input
        type="text"
        placeholder="Name of the run (optional)"
        v-model="baseForm.name"
      />
    </section>

    <section>
      <ToolTip>
        <h2>
          Host organism
          <span class="material-icons question-marks">question_mark</span>
        </h2>
        <template #tooltip>
          The organism in which the mRNA is to be expressed. The tRNA-GCN codon
          table of the host organism is utilised for codon tuning. For custom
          codon tables, head to "Codon tables" section in the navigation bar.
        </template>
      </ToolTip>

      <CodonTableSearchSelect
        v-model="selectedHostCodonTable"
        :options="codonTableList"
      />
    </section>

    <section>
      <h2>Mode</h2>

      <TuningModeSelector v-model="baseForm.mode" />
    </section>

    <section v-if="baseForm.mode == TuningModeName.PROTEIN_STRUCTURE_ANALYSIS">
      <h2>Structure of the native organism</h2>
      <PdbInput id="pdbInput" v-model="baseForm.pdb_file_content" />
    </section>

    <template v-else>
      <section>
        <h2>Sequences of native organisms</h2>

        <FastaInput
          id="fastaInput"
          v-model="baseForm.nucleotide_file_content"
        />
      </section>

      <section>
        <h2>Codon tables for the native organisms</h2>

        <table v-if="Object.keys(selectedSequencesNativeCodonTables).length">
          <thead>
            <tr>
              <th>Sequence name</th>
              <th>Codon table</th>
            </tr>
          </thead>

          <tbody>
            <tr
              v-for="seq in Object.keys(selectedSequencesNativeCodonTables)"
              :key="seq"
            >
              <td>{{ seq }}</td>
              <td class="select-cell">
                <CodonTableSearchSelect
                  v-model="selectedSequencesNativeCodonTables[seq]"
                  :options="codonTableList"
                />
              </td>
            </tr>
          </tbody>
        </table>

        <p v-else>No sequence. Please provide a valid FASTA file.</p>
      </section>
    </template>

    <section v-if="clustalIsRequired">
      <h2>Alignments</h2>

      <ClustalInput
        id="clustalInput"
        v-model="baseForm.clustal_file_content"
        :fasta-content="baseForm.nucleotide_file_content"
      />
    </section>

    <hr />

    <section>
      <h2>Thresholds</h2>

      <div class="flex-container">
        <SlowSpeedThresholdSelector v-model="baseForm.slow_speed_threshold" />
        <ConservationThresholdSelector
          v-if="baseForm.mode == TuningModeName.OPTIMISATION_AND_CONSERVATION_2"
          v-model="baseForm.conservation_threshold"
        />
        <RsaThresholdSelector
          v-else-if="baseForm.mode == TuningModeName.PROTEIN_STRUCTURE_ANALYSIS"
          v-model="baseForm.rsa_threshold"
        />
      </div>
    </section>

    <hr />

    <section>
      <h2>Specific tuning of mRNA's 5’ region</h2>

      <FivePrimeRegionTuning
        v-model="baseForm.five_prime_region_tuning"
        :mode="baseForm.mode"
      />
    </section>

    <hr />

    <section>
      <h2>Restriction enzyme recognition sites to avoid</h2>

      <div id="restrictionSiteSelector">
        <RestrictionSiteSelector v-model="baseForm.restriction_sites" />
      </div>
    </section>

    <template v-if="!isGuest">
      <hr />

      <div id="sendEmail">
        <input
          id="sendEmailCheckbox"
          type="checkbox"
          v-model="baseForm.send_email"
        />
        <label for="sendEmailCheckbox">
          Send an e-mail when the job is completed.
        </label>
      </div>
    </template>

    <hr />

    <div id="dataConsent">
      <input id="dataConsentCheckbox" type="checkbox" checked required />
      <label for="dataConsentCheckbox">
        Uploaded data are sent to our server for computing but are not shared
        with any third-party. Results are not saved unless you are logged-in. We
        take reasonable measures to protect the data you upload and store on our
        servers. While we strive to maintain a secure environment, we cannot
        guarantee absolute security. We are not liable for any unauthorized
        access, loss, or disclosure of data. By checking this box, you
        acknowledge that you have read and understand this Data Privacy Policy
        and agree to the collection, processing, and storage of your data as
        described herein.
      </label>
    </div>

    <AlertError
      :show="tuningState?.status == Status.FAILED"
      @close="tuningState = null"
    >
      {{ tuningState?.message }}
    </AlertError>

    <div id="tuningProgressWrapper">
      <p
        v-if="
          tuningState &&
          [Status.STARTED, Status.QUEUED].includes(tuningState.status)
        "
        id="tuningProgressText"
      >
        {{ tuningState.message }}
      </p>

      <template v-if="tuningState || tuningLoading">
        <ProgressBar
          id="tuningProgress"
          :value="tuningState?.step || 0"
          :max="tuningState?.total || 0"
        />
        <button
          v-if="!isGuest"
          id="cancelTuningButton"
          type="button"
          class="danger"
          @click="cancelTuning"
        >
          Cancel
        </button>
      </template>

      <button v-else id="runTuningButton" type="submit">Run tuning</button>
    </div>
  </form>
</template>

<style scoped>
#runTuningForm {
  padding-bottom: 3em;
}

section {
  margin-top: 2em;
}

.input-file {
  width: 100%;
}

.select-cell details {
  width: 100%;
  margin: 0;
}

td {
  width: 50%;
}

#codonWindow {
  margin-top: 1em;
}

#dataConsent {
  margin-bottom: 1em;
  display: flex;
  flex-direction: row;
  align-items: center;
}

#dataConsent label {
  font-size: small;
  flex: 1;
  text-align: justify;
}

#tuningProgressWrapper {
  height: 3em;
  margin-bottom: 3em;
  text-align: center;
}

#tuningProgressText {
  text-align: center;
  font-style: italic;
}

#runTuningButton {
  height: 100%;
  margin-bottom: 0;
}

#cancelTuningButton {
  height: 100%;
  margin-top: 1em;
  width: 20%;
}

#tuningProgress {
  height: 100%;
}

#tuningError {
  border: 1px solid #721c24;
  color: #721c24;
  background-color: #f8d7da;
  border-radius: 0.25rem;
  padding: 0.75rem 1.25rem;
  margin: 0.75rem 0;
}

#restrictionSiteSelector {
  margin: auto 15rem;
}

@media (max-width: 786px) {
  #cancelTuningButton {
    width: 50%;
  }
}
</style>
