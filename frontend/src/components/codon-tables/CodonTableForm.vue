<script setup lang="ts">
import { onMounted, reactive, ref, watch } from 'vue'
import { API } from '@/lib/api'
import type { CodonTable } from '@/lib/interfaces'
import { AMINO_ACID_MAPPING } from '@/lib/referentials'
import { groupByAminoAcid } from '@/lib/helpers'
import PartialCodonTable from './PartialCodonTable.vue'
import CodonTableSearchSelect from './CodonTableSearchSelect.vue'

const codonTableList = ref([] as Array<CodonTable>)
const selectedCodonTable = ref(null as CodonTable | null)

const translationMapping = reactive(AMINO_ACID_MAPPING)

onMounted(async () => await fetchCodonTableList())

watch(selectedCodonTable, async value => {
  if (value) {
    await fetchCodonTableTranslations(value.id)
  } else {
    for (const [key, value] of Object.entries(AMINO_ACID_MAPPING)) {
      translationMapping[key] = value
    }
  }
})

async function fetchCodonTableList() {
  const [data, error] = await API.codonTables.list()

  if (!error) {
    codonTableList.value = data
  }
}

async function fetchCodonTableTranslations(codonTableId: string) {
  const [data, error] = await API.codonTables.getTranslations(codonTableId)

  if (!error) {
    const newMapping = groupByAminoAcid(data)

    for (const [key, value] of Object.entries(newMapping)) {
      translationMapping[key] = value
    }
  }
}
</script>

<template>
  <CodonTableSearchSelect
    id="codonTableSelect"
    :options="codonTableList"
    v-model="selectedCodonTable"
  />

  <div class="row">
    <div class="column">
      <PartialCodonTable
        title="Alanine (Ala)"
        v-model="translationMapping.Ala"
      />
      <PartialCodonTable
        title="Asparagine (Asn)"
        v-model="translationMapping.Asn"
      />
      <PartialCodonTable
        title="Arginine (Arg)"
        v-model="translationMapping.Arg"
      />
      <PartialCodonTable
        title="Aspartic acid (Asp)"
        v-model="translationMapping.Asp"
      />
      <PartialCodonTable
        title="Cysteine (Cys)"
        v-model="translationMapping.Cys"
      />
      <PartialCodonTable
        title="Glutamic acid (Glu)"
        v-model="translationMapping.Glu"
      />
      <PartialCodonTable
        title="Glutamine (Gln)"
        v-model="translationMapping.Gln"
      />
    </div>

    <div class="column">
      <PartialCodonTable
        title="Glycine (Gly)"
        v-model="translationMapping.Gly"
      />
      <PartialCodonTable
        title="Histidine (His)"
        v-model="translationMapping.His"
      />
      <PartialCodonTable
        title="Isoleucine (Ile)"
        v-model="translationMapping.Ile"
      />
      <PartialCodonTable
        title="Leucine (Leu)"
        v-model="translationMapping.Leu"
      />
      <PartialCodonTable
        title="Lysine (Lys)"
        v-model="translationMapping.Lys"
      />
      <PartialCodonTable
        title="Methionine (Met)"
        v-model="translationMapping.Met"
      />
      <PartialCodonTable
        title="Phenylalanine (Phe)"
        v-model="translationMapping.Phe"
      />
    </div>

    <div class="column">
      <PartialCodonTable
        title="Proline (Pro)"
        v-model="translationMapping.Pro"
      />
      <PartialCodonTable
        title="Serine (Ser)"
        v-model="translationMapping.Ser"
      />
      <PartialCodonTable
        title="Threonine (Thr)"
        v-model="translationMapping.Thr"
      />
      <PartialCodonTable
        title="Tryptophan (Trp)"
        v-model="translationMapping.Trp"
      />
      <PartialCodonTable
        title="Tyrosine (Tyr)"
        v-model="translationMapping.Tyr"
      />
      <PartialCodonTable
        title="Valine (Val)"
        v-model="translationMapping.Val"
      />
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
