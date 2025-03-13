<script setup lang="ts">
import { ref, onMounted, reactive, watch, computed } from 'vue'
import { readTextFile } from '@/lib/helpers'
import type { CodonTable, RunTrainingForm } from '@/lib/interfaces'
import { API } from '@/lib/api'
import { Status, useStreamState } from '@/lib/streamedState'
import CodonTableSearchSelect from '@/components/codon-tables/CodonTableSearchSelect.vue'
import ToolTip from '@/components/ToolTip.vue'
import WithAlertError from '@/components/WithAlertError.vue'
import AlertError from '@/components/AlertError.vue'
import ProgressBar from '@/components/ProgressBar.vue'
import TuningModeSelector from '@/components/tuning-form/TuningModeSelector.vue'
import SlowSpeedThresholdSelector from '@/components/tuning-form/SlowSpeedThresholdSelector.vue'
import ConservationThresholdSelector from '@/components/tuning-form/ConservationThresholdSelector.vue'
import FivePrimeRegionTuning from '@/components/tuning-form/FivePrimeRegionTuning.vue'
import {
  checkClustal,
  checkClustalMatchingFasta,
  checkFasta,
} from '@/lib/checkers'

const emit = defineEmits(['submit'])

const baseForm: RunTrainingForm = reactive({
  name: '',
  nucleotide_file_content: '',
  clustal_file_content: '',
  host_codon_table_id: '',
  sequences_native_codon_tables: {},
  mode: 'direct_mapping',
  slow_speed_threshold: 0.5,
  conservation_threshold: null,
  five_prime_region_tuning: null,
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

const { state: tuningState, startStream: startTuningStream } = useStreamState(
  localStorage.getItem('accessToken') || undefined,
)

const clustalIsRequired = computed(() =>
  [
    'optimisation_and_conservation_1',
    'optimisation_and_conservation_2',
  ].includes(baseForm.mode),
)

onMounted(async () => await fetchCodonTables())

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

watch(
  [() => baseForm.nucleotide_file_content, () => baseForm.clustal_file_content],
  ([fastaContent, clustalContent]) => {
    baseFormErrors.clustal_file_content = []

    if (fastaContent && clustalContent) {
      const errors = checkClustalMatchingFasta(clustalContent, fastaContent)
      baseFormErrors.clustal_file_content =
        baseFormErrors.clustal_file_content.concat(errors)
    }
  },
)

watch(tuningState, state => {
  if (state && state.status == Status.FINISHED) {
    emit('submit', state.result)
  }
})

watch(
  () => baseForm.mode,
  value => {
    if (value === 'optimisation_and_conservation_2') {
      // Set a default value
      baseForm.conservation_threshold = 0.75
    } else {
      baseForm.conservation_threshold = null
    }
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

function checkFastaContent() {
  const errors = checkFasta(baseForm.nucleotide_file_content)
  baseFormErrors.nucleotide_file_content = errors
}

async function setFastaContent(event: Event) {
  baseForm.nucleotide_file_content = await readTextFile(event)
}

function checkClustalContent() {
  if (baseForm.clustal_file_content) {
    const errors = checkClustal(baseForm.clustal_file_content)
    baseFormErrors.clustal_file_content = errors
  } else {
    baseFormErrors.clustal_file_content = []
  }
}

async function setClustalContent(event: Event) {
  baseForm.clustal_file_content = await readTextFile(event)
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

    // // To remove
    // console.log(JSON.stringify(form))

    const [jobId, error] = await API.runTraining(form)

    if (!error) {
      await startTuningStream(`/api/tuning/state/${jobId}`)
    }
  }

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
      <h2>Sequences</h2>

      <div id="fastaInput">
        <ToolTip>
          <label for="fastaContent">
            Fasta format
            <span class="material-icons question-marks">question_mark</span>
          </label>
          <template #tooltip>
            A text containing the mRNA coding sequence(s), with each sequence
            preceded by a carat (">"), followed by an unique sequence
            identifier.
          </template>
        </ToolTip>
        <WithAlertError :errors="baseFormErrors.nucleotide_file_content">
          <textarea
            id="fastaContent"
            rows="10"
            spellcheck="false"
            v-model="baseForm.nucleotide_file_content"
            @input="checkFastaContent"
          ></textarea>
        </WithAlertError>

        <div class="input-file">
          <input type="file" id="fasta" @change="setFastaContent" required />

          <i>
            You can download an example sequence file
            <a href="/examples/Rad51_nucleotide.txt" download>here</a>.
          </i>
        </div>
      </div>

      <div id="clustalInput" v-if="clustalIsRequired">
        <ToolTip>
          <label for="clustalContent">
            Alignments (CLUSTAL, optional)
            <span class="material-icons question-marks">question_mark</span>
          </label>
          <template #tooltip>
            A text containing multiple sequence alignment data of orthologous
            proteins from different organisms. The alignment must include the
            amino acid sequence corresponding to the uploaded FASTA sequence,
            and must strictly follow the same order.
          </template>
        </ToolTip>
        <WithAlertError :errors="baseFormErrors.clustal_file_content">
          <textarea
            id="clustalContent"
            rows="10"
            spellcheck="false"
            v-model="baseForm.clustal_file_content"
            @input="checkClustalContent"
          ></textarea>
        </WithAlertError>

        <div class="input-file">
          <input
            type="file"
            id="clustal"
            @change="setClustalContent"
            :required="clustalIsRequired"
          />
          <i>
            You can download an example alignment file
            <a href="/examples/Rad51_CLUSTAL.txt" download>here</a>.
          </i>
        </div>
      </div>
    </section>

    <section>
      <h2>Native organisms</h2>

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

    <section>
      <h2>Mode</h2>

      <TuningModeSelector v-model="baseForm.mode" />
    </section>

    <section>
      <h2>Thresholds</h2>

      <div class="flex-container">
        <SlowSpeedThresholdSelector v-model="baseForm.slow_speed_threshold" />
        <ConservationThresholdSelector
          v-if="
            baseForm.mode == 'optimisation_and_conservation_2' &&
            baseForm.conservation_threshold !== null
          "
          v-model="baseForm.conservation_threshold"
        />
      </div>
    </section>

    <section>
      <h2>Specific tuning of mRNA's 5’ region</h2>

      <FivePrimeRegionTuning v-model="baseForm.five_prime_region_tuning" />
    </section>

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
      <ProgressBar
        v-if="tuningLoading"
        id="tuningProgress"
        :value="tuningState?.step || 0"
        :max="tuningState?.total || 0"
      />

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
}

#tuningProgressText {
  text-align: center;
  font-style: italic;
}

#runTuningButton {
  height: 100%;
  margin-bottom: 0;
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
</style>
