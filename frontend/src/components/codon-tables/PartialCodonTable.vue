<script setup lang="ts">
import PartialCodonTableRow from './PartialCodonTableRow.vue'
import type { CodonTranslation } from '@/lib/interfaces'

defineProps<{
  title: string
}>()

const model = defineModel<Array<CodonTranslation>>()
</script>

<template>
  <section>
    <h2>{{ title }}</h2>

    <table>
      <thead>
        <tr>
          <th>Anticodon</th>
          <th>Codon</th>
          <th>tRNA GCN</th>
          <th>Wobble codon</th>
          <th>Wobble rate</th>
        </tr>
      </thead>

      <tbody v-if="model">
        <PartialCodonTableRow
          v-for="(row, index) in model"
          :key="row.codon"
          :codon="row.codon"
          :anticodon="row.anticodon"
          :value="row"
          v-model:trna-gcn="model[index].trna_gcn"
          v-model:wobble-codon="model[index].wobble_codon"
          v-model:wobble-rate="model[index].wobble_rate"
        />
      </tbody>
    </table>
  </section>
</template>

<style scoped>
section {
  padding: 15px;
  border: solid 2px #c2c7d0;
  border-radius: 10px;
}

section:hover {
  border: solid 2px #5d6b89;
  border-radius: 10px;
  padding: 15px;
}

@media (prefers-color-scheme: dark) {
  section {
    border: solid 2px #5d6b89;
  }

  section:hover {
    border: solid 2px #c2c7d0;
  }
}

th {
  white-space: nowrap;
}
</style>
