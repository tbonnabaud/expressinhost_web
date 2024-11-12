<script setup lang="ts">
import { computed } from 'vue'
import PartialCodonTableRow from './PartialCodonTableRow.vue'
import { BASE_CODON_TABLE } from '@/lib/constants'

const props = defineProps<{
  title: string
  aminoAcid: string
}>()

const rowList = computed(() =>
  BASE_CODON_TABLE.filter(e => e.aminoAcid === props.aminoAcid),
)
</script>

<template>
  <section>
    <h2>{{ title }}</h2>

    <table>
      <thead>
        <tr>
          <th>anticodon</th>
          <th>codon</th>
          <th>trna_gcn</th>
          <th>corresp_codon</th>
          <th>wobble_rate</th>
        </tr>
      </thead>

      <tbody>
        <PartialCodonTableRow
          v-for="row in rowList"
          :anticodon="row.anticodon"
          :codon="row.codon"
          :key="row.codon"
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
</style>
