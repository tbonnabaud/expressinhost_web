<script setup lang="ts">
import { toFixedFloat } from '@/lib/helpers'
import ToolTip from '@/components/ToolTip.vue'

const model = defineModel<number | null>()
</script>

<template>
  <div class="input-range">
    <ToolTip>
      <label v-if="model">
        RSA threshold =
        {{ toFixedFloat(model * 100, 1) }}%
        <span class="material-icons question-marks">question_mark</span>
      </label>
      <template #tooltip>
        RSA is the Relative Solvent Accessibility. If the RSA is greater than
        the specified threshold, then the amino-acid is considered to be outside
        the structure and therefore slow.
      </template>
    </ToolTip>
    <input type="range" min="0" max="0.5" step="0.01" v-model="model" />
  </div>
</template>
