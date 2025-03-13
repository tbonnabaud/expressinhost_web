<script setup lang="ts">
import { ref, watch } from 'vue'
import { checkClustal, checkClustalMatchingFasta } from '@/lib/checkers'
import { readTextFile } from '@/lib/helpers'
import ToolTip from '@/components/ToolTip.vue'
import WithAlertError from '@/components/WithAlertError.vue'

const props = defineProps<{ fastaContent: string }>()
const model = defineModel<string | null>()
const clustalContent = ref('')
const errors = ref<string[]>([])
const matchingErrors = ref<string[]>([])

watch(
  [() => props.fastaContent, () => clustalContent.value],
  ([fasta, clustal]) => {
    if (fasta && clustal) {
      const errors = checkClustalMatchingFasta(clustal, fasta)
      matchingErrors.value = errors
    }
  },
)

watch(clustalContent, value => {
  model.value = value
})

function checkClustalContent() {
  if (clustalContent.value) {
    errors.value = checkClustal(clustalContent.value)
  }
}

async function setClustalContent(event: Event) {
  clustalContent.value = await readTextFile(event)
}
</script>

<template>
  <div id="fastaInput">
    <ToolTip>
      <label for="fastaContent">
        Fasta format
        <span class="material-icons question-marks">question_mark</span>
      </label>
      <template #tooltip>
        A text containing the mRNA coding sequence(s), with each sequence
        preceded by a carat (">"), followed by an unique sequence identifier.
      </template>
    </ToolTip>
    <WithAlertError :errors="errors.concat(matchingErrors)">
      <textarea
        id="fastaContent"
        rows="10"
        spellcheck="false"
        v-model="clustalContent"
        @input="checkClustalContent"
      ></textarea>
    </WithAlertError>

    <div class="input-file">
      <input type="file" id="fasta" @change="setClustalContent" required />

      <i>
        You can download an example sequence file
        <a href="/examples/Rad51_nucleotide.txt" download>here</a>.
      </i>
    </div>
  </div>
</template>
