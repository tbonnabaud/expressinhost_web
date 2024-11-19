<script setup lang="ts">
import { onMounted, reactive, ref, watch } from 'vue'
import { API } from '@/lib/api'
import type { CodonTable, CodonTranslation } from '@/lib/interfaces'
import { AMINO_ACID_MAPPING } from '@/lib/referentials'
import { groupByAminoAcid } from '@/lib/helpers'
import PartialCodonTable from './PartialCodonTable.vue'
import CodonTableSearchSelect from './CodonTableSearchSelect.vue'

const codonTableList = ref([] as Array<CodonTable>)
const selectedCodonTable = ref(null as CodonTable | null)

// Make a deep copy to avoid referential modification
const translationMapping = reactive(
  JSON.parse(JSON.stringify(AMINO_ACID_MAPPING)),
)

onMounted(async () => fetchCodonTableList())

watch(selectedCodonTable, async value => {
  if (value) {
    await fetchCodonTableTranslations(value.id)
  } else {
    fillTranslationMapping(AMINO_ACID_MAPPING)
  }
})

function fillTranslationMapping(mapping: Record<string, CodonTranslation[]>) {
  for (const [key, value] of Object.entries(mapping)) {
    translationMapping[key] = value
  }
}

async function fetchCodonTableList() {
  const [data, error] = await API.codonTables.list()

  if (!error && data) {
    codonTableList.value = data
  }
}

async function fetchCodonTableTranslations(codonTableId: string) {
  const [data, error] = await API.codonTables.getTranslations(codonTableId)

  if (!error) {
    const newMapping = groupByAminoAcid(data)
    fillTranslationMapping(newMapping)
  } else {
    // Reset with default values
    fillTranslationMapping(AMINO_ACID_MAPPING)
  }
}
</script>

<template>
  <div class="flex-container actions">
    <CodonTableSearchSelect
      id="codonTableSelect"
      :options="codonTableList"
      v-model="selectedCodonTable"
    />

    <div class="action-button-group">
      <button>New</button>
      <button>Save</button>
      <button>Save as...</button>
      <button class="secondary">Delete</button>
    </div>
  </div>

  <div class="grid">
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
#codonTableSelect {
  width: 50%;
}

.action-button-group {
  width: 50%;
  display: inline-flex;
  margin-bottom: 2em;
}

.action-button-group > button {
  margin-left: 1em;
  flex: 1 1 auto;
}

.actions {
  margin-bottom: 1em;
}
</style>
