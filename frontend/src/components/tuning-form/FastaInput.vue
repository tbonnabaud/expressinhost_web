<script setup lang="ts">
import { ref, watch } from 'vue'
import { checkFasta } from '@/lib/checkers'
import { readTextFile } from '@/lib/helpers'
import ToolTip from '@/components/ToolTip.vue'
import WithAlertError from '@/components/WithAlertError.vue'

const model = defineModel<string>()
const fastaContent = ref('')
const errors = ref<string[]>([])

watch(fastaContent, content => {
  model.value = content

  if (content) {
    errors.value = checkFasta(content)
  }
})

async function setFastaContent(event: Event) {
  fastaContent.value = await readTextFile(event)
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
    <WithAlertError :errors="errors">
      <textarea
        id="fastaContent"
        rows="10"
        spellcheck="false"
        v-model="fastaContent"
        placeholder="Put your FASTA sequences here or import a file using the button below."
        required
      ></textarea>
    </WithAlertError>

    <div class="input-file">
      <input type="file" id="fasta" @change="setFastaContent" />
    </div>
  </div>
</template>

<style scoped>
textarea {
  font-family: 'Courier New', Courier, monospace;
}

@media (max-width: 1024px) {
  textarea {
    font-size: 90%;
  }
}
</style>
