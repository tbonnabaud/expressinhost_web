<script setup lang="ts">
import { watch } from 'vue'
import { CODON_LIST } from '@/lib/referentials'
import type { CodonTranslation } from '@/lib/interfaces'

const model = defineModel<CodonTranslation>({
  default: {
    codon: '',
    anticodon: '',
    amino_acid: '',
    trna_gcn: 1,
    corresp_codon: 'OOO',
    wobble_rate: 0,
  },
})
const trnaGcn = defineModel('trnaGcn', { default: 1 })
const correspCodon = defineModel('correspCodon', { default: 'OOO' })
const wobbleRate = defineModel('wobbleRate', { default: 0 })

watch(trnaGcn, value => {
  if (value == 0) {
    correspCodon.value = 'GCU'
    wobbleRate.value = 0.35
  }
})

watch(correspCodon, value => {
  if (value == 'OOO') {
    wobbleRate.value = 0
  } else {
    trnaGcn.value = 0
  }
})

watch(wobbleRate, value => {
  if (value == 0) {
    correspCodon.value = 'OOO'
  } else {
    trnaGcn.value = 0
  }
})
</script>

<template>
  <tr>
    <td>{{ model.anticodon }}</td>
    <td>{{ model.codon }}</td>
    <td>
      <input type="number" v-model.number="model.trna_gcn" min="0" required />
    </td>
    <td>
      <select v-model="model.corresp_codon" required>
        <option value="OOO">OOO</option>
        <option v-for="codon in CODON_LIST" :value="codon" :key="codon">
          {{ codon }}
        </option>
      </select>
    </td>
    <td>
      <input
        type="number"
        v-model.number="model.wobble_rate"
        min="0"
        max="1"
        step="0.01"
        required
      />
    </td>
  </tr>
</template>

<style scoped>
td {
  margin: 0;
}

input,
select {
  margin: 0;
}
</style>
