<script setup lang="ts">
import CodonTableSearchSelect from './codon-tables/CodonTableSearchSelect.vue'
import SequenceInput from '@/components/SequenceInput.vue'
import { onMounted, reactive, ref, watch } from 'vue'
import type { CodonTable, SequenceComparatorForm } from '@/lib/interfaces'
import { checkSequence } from '@/lib/checkers'
import { API } from '@/lib/api'

const emit = defineEmits(['done'])

const form = reactive<SequenceComparatorForm>({
  sequence1: '',
  sequence2: '',
  host_codon_table_id: '',
})

const formErrors = reactive({
  sequence1: [] as string[],
  sequence2: [] as string[],
})

const isFormValid = ref(false)
const selectedCodonTable = ref<CodonTable | null>(null)
const codonTableList = ref<CodonTable[]>([])

onMounted(fetchCodonTables)

watch(selectedCodonTable, value => {
  if (value !== null) {
    form.host_codon_table_id = value.id
  }
})

watch(form, value => {
  formErrors.sequence1 = checkSequence(value.sequence1)
  formErrors.sequence2 = checkSequence(value.sequence2)

  if (value.sequence1.length === 0) {
    formErrors.sequence1.push('Is empty.')
  }

  if (value.sequence2.length === 0) {
    formErrors.sequence2.push('Is empty.')
  }

  if (value.sequence1.length !== value.sequence2.length) {
    formErrors.sequence2.push('The two sequences have not the same length.')
  }
})

watch(formErrors, value => {
  let isValid = true

  for (const errors of Object.values(value)) {
    if (errors.length > 0) {
      isValid = false
      break
    }
  }

  isFormValid.value = isValid
})

async function fetchCodonTables() {
  const [data, error] = await API.codonTables.list()

  if (!error) {
    codonTableList.value = data
  }
}

async function runComparison() {
  if (form.host_codon_table_id) {
    const [data, error] = await API.compareSequences(form)

    if (!error) {
      emit('done', data)
    }
  } else {
    alert('Please select a codon table.')
  }
}

function resetForm() {
  Object.assign(form, {
    sequence1: '',
    sequence2: '',
    host_codon_table_id: '',
  })
  emit('done', null)
}
</script>

<template>
  <form @submit.prevent="runComparison">
    <div>
      <label>Host organism</label>
      <CodonTableSearchSelect
        :options="codonTableList"
        v-model="selectedCodonTable"
      />
    </div>

    <div>
      <label for="sequence1Input">mRNA sequence 1</label>
      <SequenceInput
        id="sequence1Input"
        v-model="form.sequence1"
        :errors="formErrors.sequence1"
        required
      />
    </div>

    <div>
      <label for="sequence2Input">mRNA sequence 2</label>
      <SequenceInput
        id="sequence2Input"
        v-model="form.sequence2"
        :errors="formErrors.sequence2"
        required
      />
    </div>

    <button id="submitButton" type="submit" :disabled="!isFormValid">
      Compare
    </button>
    <button type="button" class="secondary" @click="resetForm">Reset</button>
  </form>
</template>

<style scoped>
#submitButton {
  margin-top: 1em;
}

.secondary {
  width: 100%;
}
</style>
