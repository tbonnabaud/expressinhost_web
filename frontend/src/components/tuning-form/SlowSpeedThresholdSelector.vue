<script setup lang="ts">
import { toFixedFloat } from '@/lib/helpers'
import ToolTip from '@/components/ToolTip.vue'

const model = defineModel<number>()
</script>

<template>
  <div class="input-range">
    <ToolTip>
      <label v-if="model">
        Conservation threshold =
        {{ toFixedFloat(model * 100, 1) }}%
        <span class="material-icons question-marks">question_mark</span>
      </label>
      <template #tooltip>
        Applies to the mode Optimisation and Conservation 2. A threshold of 0.75
        identifies a codon as conserved if 75% of the sequences in the set of
        orthologous proteins share a slow codon at the same position.
      </template>
    </ToolTip>
    <input type="range" min="0" max="1" step="0.01" v-model="model" />
  </div>
</template>
