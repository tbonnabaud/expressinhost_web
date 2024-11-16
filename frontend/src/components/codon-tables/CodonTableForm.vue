<script setup lang="ts">
import { onMounted, reactive, ref } from 'vue'
import { API } from '@/lib/api'
import type { CodonTable } from '@/lib/interfaces'
import PartialCodonTable from './PartialCodonTable.vue'
import CodonTableSearchSelect from './CodonTableSearchSelect.vue'

const codonTableList = ref([] as Array<CodonTable>)

// eslint-disable-next-line @typescript-eslint/no-unused-vars
const table = reactive({
  Ala: [],
  Arg: [],
  Asn: [],
  Asp: [],
  Cys: [],
  Gln: [],
  Glu: [],
  Gly: [],
  His: [],
  Ile: [],
  Leu: [],
  Lys: [],
  Met: [],
  Phe: [],
  Pro: [],
  Ser: [],
  Thr: [],
  Trp: [],
  Tyr: [],
  Val: [],
})

onMounted(async () => await fetchCodonTableList())

async function fetchCodonTableList() {
  const [data, error] = await API.codonTables.list()

  if (!error) {
    codonTableList.value = data
  }
}
</script>

<template>
  <CodonTableSearchSelect id="codonTableSelect" :options="codonTableList" />

  <div class="row">
    <div class="column">
      <PartialCodonTable title="Alanine (Ala)" amino-acid="Ala" />
      <PartialCodonTable title="Aspartic acid (Asp)" amino-acid="Asp" />
      <PartialCodonTable title="Glutamic acid (Glu)" amino-acid="Glu" />
      <PartialCodonTable title="Isoleucine (Ile)" amino-acid="Ile" />
      <PartialCodonTable title="Methionine (Met)" amino-acid="Met" />
      <PartialCodonTable title="Serine (Ser)" amino-acid="Ser" />
      <PartialCodonTable title="Tyrosine (Tyr)" amino-acid="Tyr" />
    </div>

    <div class="column">
      <PartialCodonTable title="Arginine (Arg)" amino-acid="Arg" />
      <PartialCodonTable title="Cysteine (Cys)" amino-acid="Cys" />
      <PartialCodonTable title="Glycine (Gly)" amino-acid="Gly" />
      <PartialCodonTable title="Leucine (Leu)" amino-acid="Leu" />
      <PartialCodonTable title="Phenylalanine (Phe)" amino-acid="Phe" />
      <PartialCodonTable title="Threonine (Thr)" amino-acid="Thr" />
    </div>

    <div class="column">
      <PartialCodonTable title="Asparagine (Asn)" amino-acid="Asn" />
      <PartialCodonTable title="Glutamine (Gln)" amino-acid="Gln" />
      <PartialCodonTable title="Histidine (His)" amino-acid="His" />
      <PartialCodonTable title="Lysine (Lys)" amino-acid="Lys" />
      <PartialCodonTable title="Proline (Pro)" amino-acid="Pro" />
      <PartialCodonTable title="Tryptophan (Trp)" amino-acid="Trp" />
      <PartialCodonTable title="Valine (Val)" amino-acid="Val" />
    </div>
  </div>
</template>

<style scoped>
.row {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  grid-gap: 20px;
}

#codonTableSelect {
  width: 50%;
  margin-bottom: 2em;
}
</style>
