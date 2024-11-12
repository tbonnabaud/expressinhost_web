<script setup lang="ts">
import { CODON_LIST } from '@/lib/constants'

defineProps<{
  anticodon: string
  codon: string
}>()

const trnaGcn = defineModel('trnaGcn', { default: 0 })
const correspCodon = defineModel('correspCodon', { default: 'OOO' })
const wobbleRate = defineModel('wobbleRate', { default: 0.35 })
</script>

<template>
  <tr>
    <td>{{ anticodon }}</td>
    <td>{{ codon }}</td>
    <td>
      <input
        type="number"
        v-model.number="trnaGcn"
        min="0"
        step="0.01"
        required
      />
    </td>
    <td>
      <select v-model="correspCodon" required>
        <option value="OOO">OOO</option>
        <option v-for="codon in CODON_LIST" :value="codon" :key="codon">
          {{ codon }}
        </option>
      </select>
    </td>
    <td>
      <input
        type="number"
        v-model.number="wobbleRate"
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
