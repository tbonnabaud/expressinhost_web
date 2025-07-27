<script setup lang="ts">
import type { TunedSequence } from '@/lib/interfaces'
import ProfileChart from './results/ProfileChart.vue'
import SequenceComparatorDownload from './SequenceComparatorDownload.vue'
import { computed } from 'vue'

const props = defineProps<{
  comparedSequences: TunedSequence
}>()

const seqComparison = computed(() =>
  colorSequences(props.comparedSequences.input, props.comparedSequences.output),
)

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
  <div>
    <h2>Result</h2>

    <h3>{{ comparedSequences.name }}</h3>

    <p>
      <strong>&rarr; Similarity: </strong>
      {{ comparedSequences.identity_percentage.toFixed() }}%
    </p>

    <div class="flex-container sequence-comparison">
      <div class="sequence-group-label">
        <label>Sequence 1:</label>
        <label>Sequence 2:</label>
      </div>

      <div class="sequence-group">
        <div class="sequence input-sequence" v-html="seqComparison.input"></div>
        <div
          class="sequence output-sequence"
          v-html="seqComparison.output"
        ></div>
      </div>
    </div>

    <div>
      <div class="profiles">
        <h5 class="profile-title">Speed profiles</h5>
        <ProfileChart
          v-if="comparedSequences.input_profiles"
          title="Speed profiles"
          input-title="Sequence 1"
          :input-values="comparedSequences.input_profiles.speed"
          output-title="Sequence 2"
          :output-values="comparedSequences.output_profiles.speed"
        />
      </div>

      <div class="profiles">
        <h5 class="profile-title">Rank profiles</h5>
        <ProfileChart
          v-if="comparedSequences.input_profiles"
          title="Rank profiles"
          input-title="Sequence 1"
          :input-values="comparedSequences.input_profiles.rank"
          output-title="Sequence 2"
          :output-values="comparedSequences.output_profiles.rank"
        />
      </div>
    </div>

    <SequenceComparatorDownload :compared-sequences="comparedSequences" />
  </div>
</template>

<style scoped>
.profiles {
  margin-top: 1em;
}

.profiles .profile-title {
  margin-bottom: 0;
}

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

.warning-icons {
  font-size: 0.75rem;
}
</style>
