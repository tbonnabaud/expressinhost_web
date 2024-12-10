<script setup lang="ts">
import { ref, onMounted, reactive, watch, computed } from 'vue'
import { readTextFile, toFixedFloat } from '@/lib/helpers'
import type { CodonTable } from '@/lib/interfaces'
import { API } from '@/lib/api'
import CodonTableSearchSelect from '@/components/codon-tables/CodonTableSearchSelect.vue'

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

const selectedHostCodonTable = ref(null as CodonTable | null)
const selectedSequencesNativeCodonTables = ref(
  {} as Record<string, CodonTable | null>,
)

const codonTableList = ref([] as Array<CodonTable>)
const tuningLoading = ref(false)

const clustalIsRequired = computed(() => baseForm.mode != 'direct_mapping')

onMounted(async () => await fetchCodonTables())

watch(
  () => baseForm.nucleotide_file_content,
  content => {
    // Reset object
    selectedSequencesNativeCodonTables.value = {}
    const sequenceNames = parseFastaSequenceNames(content)

    for (const seq of sequenceNames) {
      const table = findCorrespondingTable(seq)
      selectedSequencesNativeCodonTables.value[seq] = table
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
  const fastaSeqNameRegex = /^\>\s*(.*\w)/gm
  return Array.from(content.matchAll(fastaSeqNameRegex), m => m[1])
}

async function setFastaContent(event: Event) {
  baseForm.nucleotide_file_content = await readTextFile(event)
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

    // To remove
    console.log(JSON.stringify(form))

    const [data, error] = await API.runTraining(form)

    if (!error) {
      emit('submit', data)
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
      <h2>Host organism</h2>
      <CodonTableSearchSelect
        v-model="selectedHostCodonTable"
        :options="codonTableList"
      />
    </section>

    <section>
      <h2>Data files</h2>

      <div class="flex-container">
        <div class="input-file">
          <label for="fasta">Sequence file (FASTA)</label>
          <input type="file" id="fasta" @change="setFastaContent" required />
          <i>
            You can download an example sequence file
            <a href="/examples/Rad51_nucleotide.txt" download>here</a>.
          </i>
        </div>

        <div class="input-file">
          <label for="clustal">Alignment file (CLUSTAL, optional)</label>
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
          type="radio"
          value="direct_mapping"
          v-model="baseForm.mode"
          required
        />
        <label>Direct mapping </label>
        <input
          type="radio"
          value="optimisation_and_conservation_1"
          v-model="baseForm.mode"
        />
        <label>Optimisation and conservation 1 </label>
        <input
          type="radio"
          value="optimisation_and_conservation_2"
          v-model="baseForm.mode"
        />
        <label>Optimisation and conservation 2 </label>
      </div>
    </section>

    <section>
      <h2>Thresholds</h2>

      <div class="flex-container">
        <div class="input-range">
          <label
            >Slow speed threshold =
            {{ toFixedFloat(baseForm.slow_speed_threshold * 100, 1) }}%</label
          >
          <input
            type="range"
            min="0"
            max="1"
            step="0.01"
            v-model="baseForm.slow_speed_threshold"
          />
        </div>

        <div class="input-range">
          <label
            >Conservation threshold =
            {{ toFixedFloat(baseForm.conservation_threshold * 100, 1) }}%</label
          >
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

    <button :aria-busy="tuningLoading" id="runTuningButton" type="submit">
      {{ tuningLoading ? 'Please wait...' : 'Run tuning' }}
    </button>
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

#runTuningButton {
  width: 100%;
}

#mode-selector {
  text-align: left;
}
</style>
