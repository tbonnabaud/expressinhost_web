<script setup lang="ts">
import { computed, onMounted, reactive, ref, watch, useTemplateRef } from 'vue'
import { API } from '@/lib/api'
import type {
  CodonTable,
  CodonTableForm,
  CodonTranslation,
} from '@/lib/interfaces'
import { AMINO_ACID_MAPPING } from '@/lib/referentials'
import { groupByAminoAcid } from '@/lib/helpers'
import PartialCodonTable from './PartialCodonTable.vue'
import CodonTableSearchSelect from './CodonTableSearchSelect.vue'
import BaseModal from '../BaseModal.vue'

const codonTableList = ref([] as Array<CodonTable>)
const selectedCodonTable = ref(null as CodonTable | null)

const codonTableOrganism = ref('')
const codonTableName = ref('')

const isEditable = ref(false)
const openDeleteModal = ref(false)

const formRef = useTemplateRef('table-form')

// Make a deep copy to avoid referential modification
const translationMapping = reactive(
  JSON.parse(JSON.stringify(AMINO_ACID_MAPPING)) as Record<
    string,
    CodonTranslation[]
  >,
)

onMounted(async () => fetchCodonTableList())

const existingOrganisms = computed(() => [
  ...new Set(codonTableList.value.map(e => e.organism)),
])

watch(selectedCodonTable, async value => {
  if (value) {
    isEditable.value = value.user_id ? true : false
    codonTableOrganism.value = value.organism
    codonTableName.value = value.name
    await fetchCodonTableTranslations(value.id)
  } else {
    isEditable.value = false
    resetToDefault()
  }
})

function fillTranslationMapping(mapping: Record<string, CodonTranslation[]>) {
  for (const [key, value] of Object.entries(mapping)) {
    translationMapping[key] = value
  }
}

function resetToDefault() {
  codonTableOrganism.value = ''
  codonTableName.value = ''
  fillTranslationMapping(AMINO_ACID_MAPPING)
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
    resetToDefault()
  }
}

function prepareForm(): CodonTableForm | null {
  if (formRef.value?.checkValidity()) {
    return {
      organism: codonTableOrganism.value,
      name: codonTableName.value,
      translations: Object.values(translationMapping).flat(),
    }
  }

  return null
}

async function addNewCodonTable() {
  const form = prepareForm()

  if (form) {
    const [data, error] = await API.codonTables.add(form)

    if (!error && data) {
      console.log('Codon table added succesfully.')
      alert('Codon table added succesfully.')
      // Refresh list
      await fetchCodonTableList()
    }
  }
}

async function updateExistingCodonTable() {
  const form = prepareForm()

  if (form && selectedCodonTable.value) {
    const [, error] = await API.codonTables.update(
      selectedCodonTable.value.id,
      form,
    )

    if (!error) {
      console.log('Codon table updated succesfully.')
      alert('Codon table updated succesfully.')
      await fetchCodonTableList()
    }
  }
}

async function deleteCodonTable() {
  if (selectedCodonTable.value) {
    const [, error] = await API.codonTables.delete(selectedCodonTable.value.id)

    if (!error) {
      console.log('Codon table deleted succesfully.')
      // Refresh list
      selectedCodonTable.value = null
      await fetchCodonTableList()
    }
  }

  openDeleteModal.value = false
}
</script>

<template>
  <BaseModal
    :open="openDeleteModal"
    title="Confirm the deletion"
    @close="openDeleteModal = false"
  >
    <p v-if="selectedCodonTable">
      Do you really want to delete the codon table
      <strong>
        <i>{{ selectedCodonTable.organism }}</i> -
        {{ selectedCodonTable.name }} </strong
      >?
    </p>

    <p v-else>No codon table selected.</p>

    <footer>
      <button class="secondary" @click="openDeleteModal = false">Cancel</button>
      <button class="danger" @click="deleteCodonTable">Delete</button>
    </footer>
  </BaseModal>

  <form @submit.prevent @keydown.enter.prevent ref="table-form">
    <div id="actions" class="flex-container">
      <label id="codonTableSelect">
        Existing codon table
        <CodonTableSearchSelect
          :options="codonTableList"
          v-model="selectedCodonTable"
        />
      </label>

      <fieldset class="grid table-infos">
        <label>
          Organism
          <input
            type="text"
            placeholder="Organism"
            v-model="codonTableOrganism"
            list="existingOrganisms"
            required
          />
        </label>
        <datalist id="existingOrganisms">
          <option
            v-for="org in existingOrganisms"
            :value="org"
            :key="org"
          ></option>
        </datalist>

        <label>
          Table name
          <input
            type="text"
            placeholder="Table name"
            v-model="codonTableName"
            required
          />
        </label>

        <div class="action-button-group">
          <button @click="resetToDefault">New</button>
          <button :disabled="!isEditable" @click="updateExistingCodonTable">
            Update
          </button>
          <button @click="addNewCodonTable">Save as new</button>
          <button
            class="danger"
            :disabled="!isEditable"
            @click="openDeleteModal = true"
          >
            Delete
          </button>
        </div>
      </fieldset>
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
  </form>
</template>

<style scoped>
#codonTableSelect {
  width: 50%;
}

.table-infos {
  margin-left: 1em;
}

.action-button-group {
  display: inline-flex;
  /* Same value as for input tag */
  min-height: calc(
    1rem * var(--pico-line-height) + var(--pico-form-element-spacing-vertical) *
      2 + var(--pico-border-width) * 2
  );
  margin-top: auto;
}

.action-button-group > button:not(:first-child) {
  margin-left: 1em;
  flex: 1 1 auto;
}

#actions {
  margin-bottom: 1em;
}

#actions input {
  margin-bottom: 0;
}

#actions label {
  margin-bottom: 0;
}
</style>
