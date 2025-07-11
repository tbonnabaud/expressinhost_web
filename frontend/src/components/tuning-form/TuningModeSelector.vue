<script setup lang="ts">
import ToolTip from '@/components/ToolTip.vue'
import { TuningModeName } from '@/lib/interfaces'

const model = defineModel<TuningModeName>()
</script>

<template>
  <div id="tuningModeSelector">
    <div>
      <input
        id="protein_structure_analysis"
        type="radio"
        :value="TuningModeName.PROTEIN_STRUCTURE_ANALYSIS"
        v-model="model"
        required
      />
      <label for="protein_structure_analysis">
        <ToolTip>
          Protein structure analysis
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            Analyses the 3D structure of the protein to determine the degree of
            burial or exposure of each amino acid residue and selects codons
            based on these values (Relative Solvent Accessibility, RSA).
          </template>
        </ToolTip>
      </label>
    </div>

    <div>
      <input
        id="direct_mapping"
        type="radio"
        :value="TuningModeName.DIRECT_MAPPING"
        v-model="model"
        required
      />
      <label for="direct_mapping">
        <ToolTip>
          Direct mapping
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            Mimics the translation speed profile from the native organism into
            the host organism.
          </template>
        </ToolTip>
      </label>
    </div>

    <div>
      <input
        id="optimisation_and_conservation_1"
        type="radio"
        :value="TuningModeName.OPTIMISATION_AND_CONSERVATION_1"
        v-model="model"
      />
      <label for="optimisation_and_conservation_1">
        <ToolTip>
          Optimisation and conservation 1
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            Performs a protein sequence similarity analysis to identify
            conserved amino acids across a set of orthologous proteins from
            different organisms. The translation speed profile is maximised
            excepted at the conserved positions. Requirement: CLUSTAL alignment
            file + mRNA sequences of orthologous proteins.
          </template>
        </ToolTip>
      </label>
    </div>

    <div>
      <input
        id="optimisation_and_conservation_2"
        type="radio"
        :value="TuningModeName.OPTIMISATION_AND_CONSERVATION_2"
        v-model="model"
      />
      <label for="optimisation_and_conservation_2">
        <ToolTip>
          Optimisation and conservation 2
          <span class="material-icons question-marks">question_mark</span>
          <template #tooltip>
            Individually analyses the translation speed profile of each sequence
            in the set of orthologous proteins from different organisms, and
            identifies conserved slow translation codons. The translation speed
            profile is maximised excepted at the conserved positions.
            Requirement: CLUSTAL alignment file + mRNA sequences of orthologous
            proteins.
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
  text-align: center;
}

@media (max-width: 1024px) {
  #tuningModeSelector {
    flex-direction: column;
    row-gap: 1em;
  }
}
</style>
