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
            based on these values (Relative Solvent Accessibility, RSA). Only
            one sequence at a time can be tuned.
            <br /><strong>Requirement:</strong> .pdb file (Protein Data Bank).
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
            the host one. Any number of sequences can be tuned simultaneously.
            <br /><strong>Requirement:</strong> mRNA sequences (FASTA format)
            from the native organisms.
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
            Utilises a protein sequence similarity analysis to identify
            conserved amino acids across a set of orthologous proteins from
            different organisms. The translation speed profile of the tuned mRNA
            is maximised except at the conserved amino acid positions, where it
            mimics its native speed. Any number of sequences can be in the
            CLUSTAL amino acid alignment, that must include the one(s) for the
            mRNA(s) to be tuned.
            <br /><strong>Requirement:</strong> mRNA sequences (FASTA format)
            from the native organism(s) + CLUSTAL amino acid alignment.
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
            Analyses the translation speed profiles of the mRNA sequences
            aligned according to the corresponding CLUSTAL amino acid alignment,
            and identifies conserved slow translation codons. The translation
            speed profiles of all mRNAs are maximised except at the conserved
            slow positions, where they mimic their specific native speeds. Any
            number of sequences can be tuned simultaneously.
            <br /><strong>Requirement:</strong> mRNA sequences (FASTA format)
            from the native organisms + CLUSTAL alignment file of the
            corresponding amino acid sequences.
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
