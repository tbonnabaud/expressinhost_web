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

watch(clustalContent, content => {
  model.value = content

  if (content) {
    errors.value = checkClustal(content)
  }
})

watch(
  [() => props.fastaContent, () => clustalContent.value],
  ([fasta, clustal]) => {
    matchingErrors.value = checkClustalMatchingFasta(clustal, fasta)
  },
)

async function setClustalContent(event: Event) {
  clustalContent.value = await readTextFile(event)
}
</script>

<template>
  <div id="fastaInput">
    <ToolTip>
      <label for="fastaContent">
        CLUSTAL format
        <span class="material-icons question-marks">question_mark</span>
      </label>
      <template #tooltip>
        A text containing multiple sequence alignment data of orthologous
        proteins from different organisms. The alignment must include the amino
        acid sequence corresponding to the uploaded FASTA sequence, and must
        strictly follow the same order.
      </template>
    </ToolTip>
    <WithAlertError :errors="errors.concat(matchingErrors)">
      <textarea
        id="fastaContent"
        rows="10"
        spellcheck="false"
        v-model="clustalContent"
        placeholder="Put your CLUSTAL here or import a file using the button below."
        required
      ></textarea>
    </WithAlertError>

    <div class="input-file">
      <input type="file" id="clustal" @change="setClustalContent" />
    </div>
  </div>
</template>

<style scoped>
textarea {
  font-family: 'Courier New', Courier, monospace;
}

@media (max-width: 1024px) {
  textarea {
    font-size: 65%;
  }
}
</style>
