<script setup lang="ts">
import { toFixedFloat } from '@/lib/helpers'
import ToolTip from '@/components/ToolTip.vue'

const model = defineModel<number | null>()
</script>

<template>
  <div class="input-range">
    <ToolTip>
      <label v-if="model">
        Slow speed threshold =
        {{ toFixedFloat(model * 100, 1) }}%
        <span class="material-icons question-marks">question_mark</span>
      </label>
      <template #tooltip>
        The speed index range spans from the slowest codon to the median
        translation speed for the native mRNA sequence. A threshold of 0.5
        selects codons within the slowest 50% of this range.
      </template>
    </ToolTip>
    <input type="range" min="0" max="1" step="0.01" v-model="model" />
  </div>
</template>
