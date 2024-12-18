<script setup lang="ts">
import { computed, onMounted, ref } from 'vue'
import type { CodonTable, TunedSequence } from '@/lib/interfaces'
import { API } from '@/lib/api'
import ProfileChart from './ProfileChart.vue'

const props = defineProps<{
  tunedSequence: TunedSequence
  nativeCodonTableId: string
  open?: boolean
}>()

const openDetails = ref(false)
const nativeCodonTable = ref<CodonTable | null>(null)

const seqComparison = computed(() =>
  colorSequences(props.tunedSequence.input, props.tunedSequence.output),
)

onMounted(async () => {
  openDetails.value = Boolean(props.open)
  await fetchNativeCodonTable()
})

function toggleDetails() {
  openDetails.value = !openDetails.value
}

async function fetchNativeCodonTable() {
  const [data, error] = await API.codonTables.get(props.nativeCodonTableId)

  if (!error && data) {
    nativeCodonTable.value = data
  }
}

/**
 * Compare two arrays and return a new one with the elements of the first array surrounded
 * by span tags and "modified-codon" class if value different from the other table.
 * @param arr1
 * @param arr2
 */
function compareAndMapToArrayOfSpan(arr1: Array<string>, arr2: Array<string>) {
  return arr1.map((item, i) => {
    if (item != arr2[i]) {
      return `<span data-tooltip="${i}" data-placement="bottom" class="codon modified-codon">${item}</span>`
    }
    return `<span class="codon" data-tooltip="${i}" data-placement="bottom">${item}</span>`
  })
}

function colorSequences(inputSequence: string, outputSequence: string) {
  const inputCodonArray = inputSequence.match(/.{3}/g) || []
  const outputCodonArray = outputSequence.match(/.{3}/g) || []

  return {
    input: compareAndMapToArrayOfSpan(inputCodonArray, outputCodonArray).join(
      '',
    ),
    output: compareAndMapToArrayOfSpan(outputCodonArray, inputCodonArray).join(
      '',
    ),
  }
}
</script>

<template>
  <details :open="openDetails" @click.prevent>
    <summary @click="toggleDetails">
      <div>
        {{ tunedSequence.name }}
      <span v-if="nativeCodonTable">
        (<i>{{ nativeCodonTable.organism }}</i> - {{ nativeCodonTable.name }})
      </span>
      </div>
    </summary>

    <p>
      &rarr; Similarity: {{ tunedSequence.identity_percentage.toFixed(2) }}%
    </p>

    <div class="flex-container sequence-comparison">
      <div class="sequence-group-label">
        <label>Input: </label>
        <label>Output: </label>
      </div>

      <div class="sequence-group">
        <div class="sequence input-sequence" v-html="seqComparison.input"></div>
        <div
          class="sequence output-sequence"
          v-html="seqComparison.output"
        ></div>
      </div>
    </div>

    <div v-if="openDetails">
      <ProfileChart
        title="Speed profiles"
        :input-values="tunedSequence.input_profiles.speed"
        :output-values="tunedSequence.output_profiles.speed"
      />

      <ProfileChart
        title="Rank profiles"
        :input-values="tunedSequence.input_profiles.rank"
        :output-values="tunedSequence.output_profiles.rank"
      />
    </div>
  </details>

  <hr />
</template>

<style scoped>
.sequence {
  height: 2.5em;
  white-space: nowrap;
  font-family: monospace;
  font-size: 1.5em;
}

.sequence-group {
  overflow-x: scroll;
  /* overflow-y: visible; */
  margin: 1em 0 0 1em;
}

.sequence-group-label {
  margin: 1em 0;
}

.sequence-group-label label {
  height: 2.5em;
  font-weight: bold;
  text-decoration: underline;
}

.sequence-comparison {
  height: 145px;
}

details summary {
  font-size: 1.2em;
}
</style>
