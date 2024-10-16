<script setup lang="ts">
import { computed, ref, onMounted } from 'vue'
import { readTextFile } from '@/utils/helpers'
import { CODON_TABLE_LIST } from '@/utils/constants'

interface CodonTable {
  name: string
  organism: string
  custom: boolean
}

const fastaContent = ref('')
const clustalContent = ref('')

const codonTableList = ref([] as Array<CodonTable>)

const sequenceNameList = computed(() => {
  // Match the group after ">" symbol
  const fastaSeqNameRegex = /^\>\s*(.*\w)/gm
  return Array.from(fastaContent.value.matchAll(fastaSeqNameRegex), m => m[1])
})

onMounted(() => (codonTableList.value = CODON_TABLE_LIST))

async function setFastaContent(event: Event) {
  fastaContent.value = await readTextFile(event)
}

async function setClustalContent(event: Event) {
  clustalContent.value = await readTextFile(event)
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
</script>

<template>
  <form>
    <section>
      <h2>Data files</h2>

      <div class="flex-container">
        <div class="input-file">
          <label for="fasta">Sequence file (FASTA)</label>
          <input type="file" name="" id="fasta" @change="setFastaContent" />
        </div>

        <div class="input-file">
          <label for="clustal">Alignment file (CLUSTAL, optional)</label>
          <input type="file" name="" id="clustal" @change="setClustalContent" />
        </div>
      </div>
    </section>

    <section>
      <h2>Sequences</h2>

      <table v-if="sequenceNameList.length">
        <thead>
          <tr>
            <th>Sequence name</th>
            <th>Codon table</th>
          </tr>
        </thead>

        <tbody>
          <tr v-for="seq in sequenceNameList" :key="seq">
            <td>{{ seq }}</td>
            <td class="select-cell">
              <select name="" id="" required :value="selectTableName(seq)">
                <option value=""></option>
                <option
                  v-for="codonTableName in codonTableList"
                  :key="codonTableName.name"
                  :value="codonTableName.name"
                >
                  {{ codonTableName.name }}
                </option>
              </select>
            </td>
          </tr>
        </tbody>
      </table>

      <p v-else>No sequence. Please provide a valid FASTA file.</p>
    </section>

    <section>
      <h2>Host organism</h2>

      <div class="flex-container">
        <!-- <div id="organism">
          <label for="">Host organism</label>
          <select name="">
            <option value="">Escherichia_coli</option>
            <option value="">Danio_rerio</option>
          </select>
        </div> -->

        <div id="hostCodonTable">
          <!-- <label for="">Table name</label> -->
          <select name="" id="" required>
            <option value=""></option>
            <option
              v-for="codonTableName in codonTableList"
              :key="codonTableName.name"
              :value="codonTableName.name"
            >
              {{ codonTableName.name }}
            </option>
          </select>
        </div>
      </div>
    </section>

    <section>
      <h2>Mode</h2>

      <div id="mode-selector">
        <input type="radio" name="mode" value="direct_mapping" checked />
        <label>Direct mapping </label>
        <input
          type="radio"
          name="mode"
          value="optimisation_and_conservation_1"
        />
        <label>Optimisation and conservation 1 </label>
        <input
          type="radio"
          name="mode"
          value="optimisation_and_conservation_2"
        />
        <label>Optimisation and conservation 1 </label>
      </div>
    </section>

    <section>
      <h2>Thresholds</h2>

      <div class="flex-container">
        <div class="input-range">
          <label for="">Slow speed threshold</label>
          <input type="range" name="" id="" min="0" max="1" step="0.01" />
        </div>

        <div class="input-range">
          <label for="">Conservation threshold</label>
          <input type="range" name="" id="" min="0" max="1" step="0.01" />
        </div>
      </div>
    </section>

    <hr />

    <button id="run">Run tuning</button>
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

#organism,
#hostCodonTable {
  flex: auto;
  margin: 10px;
  width: 100%;
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

#run {
  width: 100%;
  /* background-color: #294f29; */
}

#mode-selector {
  text-align: left;
}
</style>
