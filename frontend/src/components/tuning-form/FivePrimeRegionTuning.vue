<script setup lang="ts">
import { ref, watch } from 'vue'
import type { PartialUntuningMode, FineTuningMode } from '@/lib/interfaces'
import ToolTip from '@/components/ToolTip.vue'

const mode = ref<string | null>(null)
const model = defineModel<PartialUntuningMode | FineTuningMode | null>()

watch(mode, value => {
  if (value === null) {
    model.value = null
  } else if (value == 'partial_untuning') {
    model.value = {
      mode: value,
      untuned_codons: 5,
    } as PartialUntuningMode
  } else if (value == 'fine_tuning') {
    model.value = {
      mode: value,
      codon_window_size: 5,
      utr: '',
    } as FineTuningMode
  } else {
    console.error(`Invalid ${value} mode.`)
  }
})
</script>

<template>
  <div id="fivePrimeRegionTuningSelector">
    <div>
      <label>
        <input type="radio" :value="null" v-model="mode" />
        <ToolTip>
          Like the whole sequence
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            Here nothing to be done, the output sequence is the one from the
            modes as they are.
          </template>
        </ToolTip>
      </label>
    </div>

    <div>
      <label>
        <input type="radio" value="partial_untuning" v-model="mode" />
        <ToolTip>
          Untuned
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            User selects a number of codons (something like a range slider –
            like for thresholds) between 0 and 50 codons. This part of the
            sequence is exactly like the one that was initially input by the
            user (no tuning has taken place).
          </template>
        </ToolTip>
      </label>
    </div>

    <div>
      <label>
        <input type="radio" value="fine_tuning" v-model="mode" />
        <ToolTip>
          Fine-tuned
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            User selects a number of codons (something like a range slider)
            between 0 and 10 codons, that defines the codon window.
          </template>
        </ToolTip>
      </label>
    </div>
  </div>

  <div
    id="partialUntuningModeOptions"
    v-if="mode == 'partial_untuning' && model"
  >
    <div class="input-range" v-if="'untuned_codons' in model">
      <ToolTip>
        <label>
          Number of untuned codons = {{ model.untuned_codons }} codons
          <span class="material-icons question-marks">question_mark</span>
        </label>
        <template #tooltip> Lorem ipsum </template>
      </ToolTip>
      <input
        type="range"
        min="1"
        max="50"
        step="1"
        v-model="model.untuned_codons"
      />
    </div>
  </div>

  <div id="fineTuningModeOptions" v-else-if="mode == 'fine_tuning' && model">
    <div class="input-range" v-if="'codon_window_size' in model">
      <ToolTip>
        <label>
          Codon window size = {{ model.codon_window_size }} codons
          <span class="material-icons question-marks">question_mark</span>
        </label>
        <template #tooltip> Lorem ipsum </template>
      </ToolTip>
      <input
        type="range"
        min="1"
        max="10"
        step="1"
        v-model="model.codon_window_size"
      />
    </div>

    <div v-if="'utr' in model">
      <label id="utrSequenceLabel" for="utrSequence">UTR sequence</label>
      <textarea
        id="utrSequence"
        v-model="model.utr"
        placeholder="Put UTR sequence"
        rows="3"
      ></textarea>
    </div>
  </div>
</template>

<style scoped>
#fivePrimeRegionTuningSelector {
  display: flex;
  column-gap: 2em;
  justify-content: center;
}

#partialUntuningModeOptions,
#fineTuningModeOptions {
  margin-top: 2em;
}

#fineTuningModeOptions {
  text-align: center;
}
</style>
