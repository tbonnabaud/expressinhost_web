<script setup lang="ts">
import { watch } from 'vue'
import { CODON_LIST } from '@/lib/referentials'

defineProps<{
  codon: string
  anticodon: string
}>()

const trnaGcn = defineModel('trnaGcn', { default: 1 })
const correspCodon = defineModel('correspCodon', { default: '---' })
const wobbleRate = defineModel('wobbleRate', { default: 0 })

watch(trnaGcn, value => {
  if (value > 0) {
    correspCodon.value = '---'
    wobbleRate.value = 0
  }
})

watch(correspCodon, value => {
  if (value == '---') {
    wobbleRate.value = 0
  } else {
    trnaGcn.value = 0
  }
})

watch(wobbleRate, value => {
  if (value == 0) {
    correspCodon.value = '---'
  } else {
    trnaGcn.value = 0
  }
})
</script>

<template>
  <tr>
    <td>{{ anticodon }}</td>
    <td>{{ codon }}</td>
    <td>
      <input type="number" v-model.number="trnaGcn" min="0" required />
    </td>
    <td>
      <select v-model="correspCodon" required>
        <option value="---">---</option>
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
