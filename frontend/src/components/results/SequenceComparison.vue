<script setup lang="ts">
import { computed, onMounted, ref } from 'vue'
import type { TunedSequence } from '@/lib/interfaces'
import ProfileChart from './ProfileChart.vue'

const props = defineProps<{
  tunedSequence: TunedSequence
  open?: boolean
}>()

const openDetails = ref(false)

const seqComparison = computed(() =>
  colorSequences(props.tunedSequence.input, props.tunedSequence.output),
)

onMounted(() => {
  openDetails.value = Boolean(props.open)
})

function toggleDetails(event: Event) {
  // Prevent the default toggle behavior of HTML <details> element
  event.preventDefault()
  openDetails.value = !openDetails.value
}

function colorSequences(inputSequence: string, outputSequence: string) {
  const inputCodonArray = inputSequence.match(/.{3}/g) || []
  const outputCodonArray = outputSequence.match(/.{3}/g) || []

  return [
    inputCodonArray
      .map((item, i) => {
        if (item != outputCodonArray[i]) {
          return `<span class="in-seq">${item}</span>`
        }
        return item
      })
      .join(''),
    outputCodonArray
      .map((item, i) => {
        if (item != inputCodonArray[i]) {
          return `<span class="out-seq">${item}</span>`
        }
        return item
      })
      .join(''),
  ]
}
</script>

<template>
  <details :open="openDetails" @click.prevent>
    <summary @click="toggleDetails">{{ tunedSequence.name }}</summary>

    <p>
      &rarr; Similarity: {{ tunedSequence.identity_percentage.toFixed(2) }}%
    </p>

    <div class="flex-container">
      <div class="sequence-group-label">
        <label>Input: </label>
        <label>Output: </label>
      </div>

      <div class="sequence-group">
        <div class="sequence" v-html="seqComparison[0]"></div>
        <div class="sequence" v-html="seqComparison[1]"></div>
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

details summary {
  font-size: 1.2em;
}
</style>
