<script setup lang="ts">
import { ref, watch } from 'vue'
import { CODON_LIST } from '@/lib/referentials'

defineProps<{
  codon: string
  anticodon: string
}>()

const trnaGcn = defineModel('trnaGcn', { default: 1 })
const wobbleCodon = defineModel('wobbleCodon', { default: '---' })
const wobbleRate = defineModel('wobbleRate', { default: 0 })

const activity = ref(100)

watch(trnaGcn, value => {
  if (value > 0) {
    wobbleCodon.value = '---'
    wobbleRate.value = 0
  }
})

watch(wobbleCodon, value => {
  if (value == '---') {
    wobbleRate.value = 0
  } else if (value != '---' && wobbleRate.value == 0) {
    wobbleRate.value = 0.35
  } else {
    trnaGcn.value = 0
  }
})

watch(wobbleRate, value => {
  activity.value = (1 - value) * 100

  if (value == 0) {
    wobbleCodon.value = '---'
  } else {
    trnaGcn.value = 0
  }
})

function updateWobbleRate() {
  // Round the number with a maximum of three decimals
  wobbleRate.value = parseFloat((1 - activity.value / 100).toFixed(3))
}
</script>

<template>
  <tr>
    <td>{{ anticodon }}</td>
    <td>{{ codon }}</td>
    <td>
      <input type="number" v-model.number="trnaGcn" min="0" required />
    </td>
    <td>
      <select v-model="wobbleCodon" required>
        <option value="---">---</option>
        <option v-for="codon in CODON_LIST" :value="codon" :key="codon">
          {{ codon }}
        </option>
      </select>
    </td>
    <td>
      <input
        type="number"
        v-model="activity"
        min="0"
        max="100"
        step="1"
        @change="updateWobbleRate"
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
