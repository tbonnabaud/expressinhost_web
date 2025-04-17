<script setup lang="ts">
import { computed, ref, watch } from 'vue'
import { store } from '@/lib/store'
import type {
  PartialUntuningMode,
  FineTuningMode,
  SlowedDownMode,
  FivePrimeRegionTuningMode,
} from '@/lib/interfaces'
import { TuningModeName } from '@/lib/interfaces'
import ToolTip from '@/components/ToolTip.vue'
import UtrSequenceInput from './five-prime-region/UtrSequenceInput.vue'

defineProps<{ mode: TuningModeName }>()

const fivePrimeMode = ref<string | null>(null)
const model = defineModel<FivePrimeRegionTuningMode | null>()

const isGuest = computed(() => store.currentUser.value === null)

watch(fivePrimeMode, value => {
  if (value === null) {
    model.value = null
  } else if (value == 'partial_untuning') {
    model.value = {
      mode: value,
      untuned_codon_number: 5,
    } as PartialUntuningMode
  } else if (value == 'fine_tuning') {
    model.value = {
      mode: value,
      codon_window_size: 5,
      utr: '',
    } as FineTuningMode
  } else if (value == 'slowed_down') {
    model.value = { mode: value, slowed_down_codon_number: 5 } as SlowedDownMode
  } else {
    console.error(`Invalid ${value} mode.`)
  }
})
</script>

<template>
  <div id="fivePrimeRegionTuningSelector">
    <div>
      <label>
        <input type="radio" :value="null" v-model="fivePrimeMode" />
        <ToolTip>
          Like the whole sequence
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            Here nothing to be done, the output sequences are the ones from the
            modes as they are.
          </template>
        </ToolTip>
      </label>
    </div>

    <div v-if="mode != TuningModeName.PROTEIN_STRUCTURE_ANALYSIS">
      <label>
        <input type="radio" value="partial_untuning" v-model="fivePrimeMode" />
        <ToolTip>
          Untuned
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            In this mode the N first selected codons of the sequences are
            exactly like the ones that was initially input by the user (no
            tuning has taken place).
          </template>
        </ToolTip>
      </label>
    </div>

    <div>
      <label>
        <input type="radio" value="slowed_down" v-model="fivePrimeMode" />
        Slowed down
        <!-- <ToolTip>
          Slowed down
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            TODO
          </template>
        </ToolTip> -->
      </label>
    </div>

    <div>
      <label :aria-disabled="isGuest">
        <input
          type="radio"
          value="fine_tuning"
          v-model="fivePrimeMode"
          :disabled="isGuest"
        />
        <span v-if="isGuest">Fine-tuned (for logged user only)</span>
        <ToolTip v-else>
          Fine-tuned
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            In this mode the sequences are optimised with OSTIR software. User
            selects a number of codons between 1 and 10 codons, that defines the
            codon window. The larger the codon window, the longer the
            calculation time, due to the greater number of combinations (of the
            order of the factorial). For example, per sequence, it may take
            about 3 seconds for a window of 5 codons, but more than half an hour
            for 10.
          </template>
        </ToolTip>
      </label>
    </div>
  </div>

  <div
    id="partialUntuningModeOptions"
    v-if="fivePrimeMode == 'partial_untuning' && model"
  >
    <div class="input-range" v-if="'untuned_codon_number' in model">
      <label>
        Number of untuned codons = {{ model.untuned_codon_number }} codons
      </label>
      <input
        type="range"
        min="1"
        max="50"
        step="1"
        v-model="model.untuned_codon_number"
      />
    </div>
  </div>

  <div
    id="partialUntuningModeOptions"
    v-else-if="fivePrimeMode == 'slowed_down' && model"
  >
    <div class="input-range" v-if="'slowed_down_codon_number' in model">
      <label>
        Number of slowed down codons =
        {{ model.slowed_down_codon_number }} codons
      </label>
      <input
        type="range"
        min="1"
        max="50"
        step="1"
        v-model="model.slowed_down_codon_number"
      />
    </div>
  </div>

  <div
    id="fineTuningModeOptions"
    v-else-if="!isGuest && fivePrimeMode == 'fine_tuning' && model"
  >
    <div class="input-range" v-if="'codon_window_size' in model">
      <label> Codon window size = {{ model.codon_window_size }} codons </label>
      <input
        type="range"
        min="1"
        max="10"
        step="1"
        v-model="model.codon_window_size"
      />
    </div>

    <div v-if="'utr' in model">
      <UtrSequenceInput v-model="model.utr" />
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
