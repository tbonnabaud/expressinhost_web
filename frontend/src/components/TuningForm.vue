<script setup lang="ts">
import { computed, ref } from 'vue'

const fastaContent = ref('')
const clustalContent = ref('')

const codonTableNameList = ref([
  'Arabidopsis_thaliana',
  'Bacillus_subtilis',
  'Caenorhabditis_elegans',
  'Danio_rerio',
  'Drosophila_melanogaster',
  'Escherichia_coli',
  'Gallus_gallus',
  'Homo_sapiens',
  'Komagataella_pastoris',
  'Methanocaldococcus_jannaschii',
  'Saccharomyces_cerevisiae',
  'Staphylococcus_aureus',
  'Xenopus_laevis',
])

const sequenceNames = computed(() => {
  const fastaSeqRegex = /^\> ?([\w ]*\w)/gm
  return Array.from(fastaContent.value.matchAll(fastaSeqRegex), m => m[1])
})

async function readText(event: Event) {
  const target = event.target as HTMLInputElement

  if (target.files && target.files.length > 0) {
    const file = target.files[0]
    const text = await file.text()
    return text
  }

  return ''
}

async function setFastaContent(event: Event) {
  fastaContent.value = await readText(event)
}

async function setClustalContent(event: Event) {
  clustalContent.value = await readText(event)
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

      <table v-if="sequenceNames.length">
        <tbody>
          <tr v-for="seq in sequenceNames" :key="seq">
            <td>{{ seq }}</td>
            <td class="select-cell">
              <select name="" id="" required>
                <option value=""></option>
                <option
                  v-for="codonTableName in codonTableNameList"
                  :key="codonTableName"
                  :value="codonTableName"
                >
                  {{ codonTableName }}
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

        <div id="tableName">
          <label for="">Table name</label>
          <select name="" id="" required>
            <option value=""></option>
            <option
              v-for="codonTableName in codonTableNameList"
              :key="codonTableName"
              :value="codonTableName"
            >
              {{ codonTableName }}
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
#tableName {
  flex: auto;
  margin: 10px;
  width: 100%;
}

.input-range {
  flex: auto;
  text-align: center;
}

td select {
  width: 100%;
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
