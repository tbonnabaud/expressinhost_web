<script setup lang="ts">
import PartialCodonTable from './PartialCodonTable.vue'
import { reactive } from 'vue'

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
</script>

<template>
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
</style>
