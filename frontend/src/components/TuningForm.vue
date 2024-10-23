<script setup lang="ts">
import { ref, onMounted, reactive, watch, computed } from 'vue'
import { readTextFile } from '@/lib/helpers'
import { CODON_TABLE_LIST } from '@/lib/constants'
import type { CodonTable, RunTrainingForm } from '@/lib/interfaces'
import { API } from '@/lib/api'
import SearchSelect from '@/components/SearchSelect.vue'

const emit = defineEmits(['submit'])

const form = reactive({
  nucleotide_file_content: '',
  clustal_file_content: '',
  host_codon_table_name: '',
  sequences_native_codon_tables: {},
  mode: 'direct_mapping',
  slow_speed_threshold: 0.5,
  conservation_threshold: 0.75,
} as RunTrainingForm)

const codonTableList = ref([] as Array<CodonTable>)

const clustalIsRequired = computed(() => form.mode != 'direct_mapping')
const codonTableNameList = computed(() => codonTableList.value.map(e => e.name))

onMounted(() => (codonTableList.value = CODON_TABLE_LIST))

watch(
  () => form.nucleotide_file_content,
  content => {
    // Reset object
    form.sequences_native_codon_tables = {}
    const sequenceNames = parseFastaSequenceNames(content)

    for (const seq of sequenceNames) {
      form.sequences_native_codon_tables[seq] = selectTableName(seq)
    }
  },
)

function parseFastaSequenceNames(content: string) {
  // Match the group after ">" symbol
  const fastaSeqNameRegex = /^\>\s*(.*\w)/gm
  return Array.from(content.matchAll(fastaSeqNameRegex), m => m[1])
}

async function setFastaContent(event: Event) {
  form.nucleotide_file_content = await readTextFile(event)
}

async function setClustalContent(event: Event) {
  form.clustal_file_content = await readTextFile(event)
}

/**
 * Try to guess the first corresponding codon table name.
 * @param {string} sequenceName - Name of the sequence
 * @returns First corresponding codon table name or empty string if not found.
 */
function selectTableName(sequenceName: string) {
  const lowerCaseSequenceName = sequenceName.toLowerCase()

  for (const tableName of codonTableList.value) {
    if (lowerCaseSequenceName.includes(tableName.name.toLowerCase())) {
      return tableName.name
    }
  }

  return ''
}

async function runTuning() {
  // console.log(JSON.stringify(form))
  const output = await API.runTraining(form)
  // console.log(JSON.stringify(output))
  emit('submit', output)
}
</script>

<template>
  <form id="runTuningForm" @submit.prevent="runTuning">
    <section>
      <h2>Data files</h2>

      <div class="flex-container">
        <div class="input-file">
          <label for="fasta">Sequence file (FASTA)</label>
          <input type="file" id="fasta" @change="setFastaContent" required />
        </div>

        <div class="input-file">
          <label for="clustal">Alignment file (CLUSTAL, optional)</label>
          <input
            type="file"
            id="clustal"
            @change="setClustalContent"
            :required="clustalIsRequired"
          />
        </div>
      </div>
    </section>

    <section>
      <h2>Sequences</h2>

      <table v-if="Object.keys(form.sequences_native_codon_tables).length">
        <thead>
          <tr>
            <th>Sequence name</th>
            <th>Codon table</th>
          </tr>
        </thead>

        <tbody>
          <tr
            v-for="seq in Object.keys(form.sequences_native_codon_tables)"
            :key="seq"
          >
            <td>{{ seq }}</td>
            <td class="select-cell">
              <SearchSelect
                v-model="form.sequences_native_codon_tables[seq]"
                :options="codonTableNameList"
              />
              <!-- <select
                required
                v-model="form.sequences_native_codon_tables[seq]"
              >
                <option value=""></option>
                <option
                  v-for="codonTableName in codonTableList"
                  :key="codonTableName.name"
                  :value="codonTableName.name"
                >
                  {{ codonTableName.name }}
                </option>
              </select> -->
            </td>
          </tr>
        </tbody>
      </table>

      <p v-else>No sequence. Please provide a valid FASTA file.</p>
    </section>

    <section>
      <h2>Host organism</h2>
      <SearchSelect
        v-model="form.host_codon_table_name"
        :options="codonTableNameList"
      />
      <!-- <div id="hostCodonTable">
        <label>Table name</label>
        <select
          id="hostCodonTableName"
          v-model="form.host_codon_table_name"
          required
        >
          <option value=""><input type="text" placeholder="coucou" /></option>
          <option
            v-for="codonTableName in codonTableList"
            :key="codonTableName.name"
            :value="codonTableName.name"
          >
            {{ codonTableName.name }}
          </option>
        </select>
      </div> -->
    </section>

    <section>
      <h2>Mode</h2>

      <div id="mode-selector">
        <input
          type="radio"
          value="direct_mapping"
          v-model="form.mode"
          required
        />
        <label>Direct mapping </label>
        <input
          type="radio"
          value="optimisation_and_conservation_1"
          v-model="form.mode"
        />
        <label>Optimisation and conservation 1 </label>
        <input
          type="radio"
          value="optimisation_and_conservation_2"
          v-model="form.mode"
        />
        <label>Optimisation and conservation 2 </label>
      </div>
    </section>

    <section>
      <h2>Thresholds</h2>

      <div class="flex-container">
        <div class="input-range">
          <label>Slow speed threshold = {{ form.slow_speed_threshold }}</label>
          <input
            type="range"
            min="0"
            max="1"
            step="0.01"
            v-model="form.slow_speed_threshold"
          />
        </div>

        <div class="input-range">
          <label
            >Conservation threshold = {{ form.conservation_threshold }}</label
          >
          <input
            type="range"
            min="0"
            max="1"
            step="0.01"
            v-model="form.conservation_threshold"
          />
        </div>
      </div>
    </section>

    <hr />

    <button id="runTuningButton" type="submit">Run tuning</button>
  </form>
</template>

<style scoped>
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
}

.select-cell select {
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
