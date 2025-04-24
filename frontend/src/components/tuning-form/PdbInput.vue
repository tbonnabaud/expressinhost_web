<script setup lang="ts">
import { ref, watch } from 'vue'
import { readTextFile } from '@/lib/helpers'
import { checkPdb } from '@/lib/checkers'
import ToolTip from '@/components/ToolTip.vue'
import WithAlertError from '@/components/WithAlertError.vue'

const model = defineModel<string | null>()
const pdbContent = ref('')
const errors = ref<string[]>([])

watch(pdbContent, content => {
  model.value = content

  if (content) {
    errors.value = checkPdb(content)
  }
})

async function setPdbContent(event: Event) {
  pdbContent.value = await readTextFile(event)
}
</script>

<template>
  <div id="pdbInput">
    <ToolTip>
      <label for="pdbContent">
        PDB format
        <span class="material-icons question-marks">question_mark</span>
      </label>
      <template #tooltip>
        A text containing the structure of a protein, in PDB format.
      </template>
    </ToolTip>
    <WithAlertError :errors="errors">
      <textarea
        id="pdbContent"
        rows="10"
        spellcheck="false"
        v-model="pdbContent"
        placeholder="Put your PDB content here or import a file using the button below."
        required
      ></textarea>
    </WithAlertError>

    <div class="input-file">
      <input type="file" id="pdb" @change="setPdbContent" />
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
