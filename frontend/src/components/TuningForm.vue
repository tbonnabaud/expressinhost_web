<script setup lang="ts">
import { ref, onMounted, reactive, watch, computed } from 'vue'
import { readTextFile, toFixedFloat } from '@/lib/helpers'
import type { CodonTable } from '@/lib/interfaces'
import { API } from '@/lib/api'
import { Status, useStreamState } from '@/lib/streamedState'
import CodonTableSearchSelect from '@/components/codon-tables/CodonTableSearchSelect.vue'
import ToolTip from '@/components/ToolTip.vue'
import WithAlertError from './WithAlertError.vue'
import ProgressBar from './ProgressBar.vue'
import {
  checkClustal,
  checkClustalMatchingFasta,
  checkFasta,
} from '@/lib/checkers'

const emit = defineEmits(['submit'])

const baseForm = reactive({
  name: '',
  nucleotide_file_content: '',
  clustal_file_content: '',
  host_codon_table_id: '',
  sequences_native_codon_tables: {} as Record<string, string>,
  mode: 'direct_mapping',
  slow_speed_threshold: 0.5,
  conservation_threshold: 0.75,
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
  '/api/run-tuning',
  'POST',
  localStorage.getItem('accessToken') || undefined,
)

const clustalIsRequired = computed(() => baseForm.mode != 'direct_mapping')

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
    if (fastaContent && clustalContent) {
      const errors = checkClustalMatchingFasta(clustalContent, fastaContent)

      if (errors.length) {
        baseFormErrors.clustal_file_content =
          baseFormErrors.clustal_file_content.concat(errors)
      }
    }
  },
)

watch(tuningState, state => {
  if (state && state.status == Status.SUCCESS) {
    emit('submit', state.result)
  }
})

async function fetchCodonTables() {
  const [data, error] = await API.codonTables.list()

  if (!error) {
    codonTableList.value = data
  }
}

function parseFastaSequenceNames(content: string) {
  // Match the group after ">" symbol
  const fastaSeqIdentifierRegex = /^\>([^\s]+)/gm

  return Array.from(content.matchAll(fastaSeqIdentifierRegex), m => m[1])
}

async function setFastaContent(event: Event) {
  // Reset error list
  baseFormErrors.nucleotide_file_content = []
  const content = await readTextFile(event)
  const errors = checkFasta(content)

  if (errors.length) {
    baseFormErrors.nucleotide_file_content = errors
  } else {
    baseForm.nucleotide_file_content = content
  }
}

async function setClustalContent(event: Event) {
  // Reset error list
  baseFormErrors.clustal_file_content = []
  const content = await readTextFile(event)
  const errors = checkClustal(content)

  if (errors.length) {
    baseFormErrors.clustal_file_content = errors
  } else {
    baseForm.clustal_file_content = content
  }
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

    // const [data, error] = await API.runTraining(form)
    await startTuningStream(form)

    // if (!error) {
    //   emit('submit', data)
    // }
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
        <h2>Host organism 🯄</h2>
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
      <h2>Data files</h2>

      <div class="flex-container">
        <div class="input-file">
          <ToolTip>
            <label for="fasta">Sequence file (FASTA) 🯄</label>
            <template #tooltip>
              A text-based file containing the mRNA coding sequence(s), with
              each sequence preceded by a carat (">"), followed by an unique
              sequence identifier.
            </template>
          </ToolTip>

          <WithAlertError :errors="baseFormErrors.nucleotide_file_content">
            <input type="file" id="fasta" @change="setFastaContent" required />
          </WithAlertError>

          <i>
            You can download an example sequence file
            <a href="/examples/Rad51_nucleotide.txt" download>here</a>.
          </i>
        </div>

        <div class="input-file">
          <ToolTip>
            <label for="clustal">Alignment file (CLUSTAL, optional) 🯄</label>
            <template #tooltip>
              A text-based file containing multiple sequence alignment data of
              orthologous proteins from different organisms. The alignment must
              include the amino acid sequence corresponding to the uploaded
              FASTA sequence, and must strictly follow the same order.
            </template>
          </ToolTip>

          <WithAlertError :errors="baseFormErrors.clustal_file_content">
            <input
              type="file"
              id="clustal"
              @change="setClustalContent"
              :required="clustalIsRequired"
            />
          </WithAlertError>
          <i>
            You can download an example alignment file
            <a href="/examples/Rad51_CLUSTAL.txt" download>here</a>.
          </i>
        </div>
      </div>
    </section>

    <section>
      <h2>Sequences</h2>

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

      <div id="mode-selector">
        <input
          id="direct_mapping"
          type="radio"
          value="direct_mapping"
          v-model="baseForm.mode"
          required
        />
        <label htmlFor="direct_mapping">
          <ToolTip>
            Direct mapping 🯄
            <template #tooltip>
              Tuning mode, mimics the translation speed profile from the native
              organism into the host organism.
            </template>
          </ToolTip>
        </label>

        <input
          id="optimisation_and_conservation_1"
          type="radio"
          value="optimisation_and_conservation_1"
          v-model="baseForm.mode"
        />
        <label htmlFor="optimisation_and_conservation_1">
          <ToolTip>
            Optimisation and conservation 1 🯄
            <template #tooltip>
              Tuning mode, performs a protein sequence similarity analysis to
              identify conserved amino acids across a set of orthologous
              proteins from different organisms. The translation speed profile
              is maximised excepted at the conserved positions. Requirement:
              CLUSTAL alignment file.
            </template>
          </ToolTip>
        </label>

        <input
          id="optimisation_and_conservation_2"
          type="radio"
          value="optimisation_and_conservation_2"
          v-model="baseForm.mode"
        />
        <label htmlFor="optimisation_and_conservation_2">
          <ToolTip>
            Optimisation and conservation 2 🯄
            <template #tooltip>
              Tuning mode, individually analyses the translation speed profile
              of each sequence in the set of orthologous proteins from different
              organisms, and identifies conserved slow translation codons. The
              translation speed profile is maximised excepted at the conserved
              positions. Requirement: CLUSTAL alignment file.
            </template>
          </ToolTip>
        </label>
      </div>
    </section>

    <section>
      <h2>Thresholds</h2>

      <div class="flex-container">
        <div class="input-range">
          <ToolTip>
            <label>
              Slow speed threshold =
              {{ toFixedFloat(baseForm.slow_speed_threshold * 100, 1) }}% 🯄
            </label>
            <template #tooltip>
              Applies to the mode Optimisation and Conservation 2. The speed
              index range spans from the slowest codon to the median translation
              speed for the native mRNA sequence. A threshold of 0.5 selects
              codons within the slowest 50% of this range.
            </template>
          </ToolTip>
          <input
            type="range"
            min="0"
            max="1"
            step="0.01"
            v-model="baseForm.slow_speed_threshold"
          />
        </div>

        <div class="input-range">
          <ToolTip>
            <label>
              Conservation threshold =
              {{ toFixedFloat(baseForm.conservation_threshold * 100, 1) }}% 🯄
            </label>
            <template #tooltip>
              Applies to the mode Optimisation and Conservation 2. A threshold
              of 0.75 identifies a codon as conserved if 75% of the sequences in
              the set of orthologous proteins share a slow codon at the same
              position.
            </template>
          </ToolTip>
          <input
            type="range"
            min="0"
            max="1"
            step="0.01"
            v-model="baseForm.conservation_threshold"
          />
        </div>
      </div>
    </section>

    <hr />

    <div id="tuningProgressWrapper">
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
section {
  margin-top: 2em;
}

.input-file {
  width: 100%;
}

.input-range input {
  width: 90%;
  margin-bottom: 1em;
}

.input-range {
  flex: auto;
  text-align: center;
  width: 100%;
}

.select-cell details {
  width: 100%;
  margin: 0;
}

td {
  width: 50%;
}

#mode-selector {
  text-align: left;
}

#runTuningButton {
  width: 100%;
  height: 100%;
}

#tuningProgressWrapper,
#tuningProgress {
  height: 100%;
}

#tuningProgressWrapper {
  height: 3em;
}
</style>
