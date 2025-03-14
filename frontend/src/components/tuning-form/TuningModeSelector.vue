<script setup lang="ts">
import ToolTip from '@/components/ToolTip.vue'

const model = defineModel<string>()
</script>

<template>
  <div id="tuningModeSelector">
    <div>
      <input
        id="direct_mapping"
        type="radio"
        value="direct_mapping"
        v-model="model"
        required
      />
      <label for="direct_mapping">
        <ToolTip>
          Direct mapping
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            Tuning mode, mimics the translation speed profile from the native
            organism into the host organism.
          </template>
        </ToolTip>
      </label>
    </div>

    <div>
      <input
        id="optimisation_and_conservation_1"
        type="radio"
        value="optimisation_and_conservation_1"
        v-model="model"
      />
      <label for="optimisation_and_conservation_1">
        <ToolTip>
          Optimisation and conservation 1
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            Tuning mode, performs a protein sequence similarity analysis to
            identify conserved amino acids across a set of orthologous proteins
            from different organisms. The translation speed profile is maximised
            excepted at the conserved positions. Requirement: CLUSTAL alignment
            file.
          </template>
        </ToolTip>
      </label>
    </div>

    <div>
      <input
        id="optimisation_and_conservation_2"
        type="radio"
        value="optimisation_and_conservation_2"
        v-model="model"
      />
      <label for="optimisation_and_conservation_2">
        <ToolTip>
          Optimisation and conservation 2
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            Tuning mode, individually analyses the translation speed profile of
            each sequence in the set of orthologous proteins from different
            organisms, and identifies conserved slow translation codons. The
            translation speed profile is maximised excepted at the conserved
            positions. Requirement: CLUSTAL alignment file.
          </template>
        </ToolTip>
      </label>
    </div>
  </div>
</template>

<style scoped>
#tuningModeSelector {
  display: flex;
  column-gap: 2em;
  justify-content: center;
}

@media (max-width: 1024px) {
  #tuningModeSelector {
    flex-direction: column;
    row-gap: 1em;
  }
}
</style>
