<script setup lang="ts">
import type { TunedSequence } from '@/lib/interfaces'
import ProfileChart from './ProfileChart.vue'
import { computed } from 'vue'

const props = defineProps<{
  tunedSequence: TunedSequence
}>()

const outputCodonList = computed(() => {
  return props.tunedSequence.output.match(/.{3}/g) || []
})
</script>

<template>
  <h4>{{ tunedSequence.name }}</h4>

  <div class="flex-container sequence-comparison">
    <div class="sequence-group-label">
      <label>Input: </label>
      <label>Output: </label>
    </div>

    <div class="sequence-group">
      <div class="sequence input-sequence">
        <span
          v-for="(item, index) in tunedSequence.input"
          :key="index"
          class="codon"
          :data-tooltip="index"
          data-placement="bottom"
        >
          &nbsp;{{ item }}
        </span>
      </div>
      <div class="sequence output-sequence">
        <span
          v-for="(item, index) in outputCodonList"
          :key="index"
          class="codon"
          :data-tooltip="index"
          data-placement="bottom"
        >
          {{ item }}
        </span>
      </div>
    </div>
  </div>

  <div class="profiles">
    <h5 class="profile-title">Speed profiles</h5>
    <ProfileChart
      title="Speed profiles"
      :input-values="[]"
      :output-values="tunedSequence.output_profiles.speed"
    />
  </div>

  <div class="profiles">
    <h5 class="profile-title">Rank profiles</h5>
    <ProfileChart
      title="Rank profiles"
      :input-values="[]"
      :output-values="tunedSequence.output_profiles.rank"
    />
  </div>
</template>

<style scoped>
.profiles {
  margin-top: 1em;
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
</style>
