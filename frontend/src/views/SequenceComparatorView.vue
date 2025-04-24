<script setup lang="ts">
import SequenceComparatorForm from '@/components/SequenceComparatorForm.vue'
import SequenceComparatorResult from '@/components/SequenceComparatorResult.vue'
import type { TunedSequence } from '@/lib/interfaces'
import { ref } from 'vue'

const comparedSequences = ref<TunedSequence | null>(null)

function setComparedSequences(value: TunedSequence | null) {
  comparedSequences.value = value
}
</script>

<template>
  <main class="container-fluid">
    <h1>Sequence comparator</h1>

    <div id="sequenceComparator">
      <SequenceComparatorForm
        id="sequenceComparatorForm"
        @done="setComparedSequences"
      />

      <div id="sequenceComparatorResult">
        <SequenceComparatorResult
          v-if="comparedSequences"
          :compared-sequences="comparedSequences"
        />
      </div>
    </div>
  </main>
</template>

<style scoped>
#sequenceComparator {
  display: flex;
  gap: 2em;
  height: 100%;
}

#sequenceComparatorForm {
  flex: 2;
  padding: 1em 2em;
  border: 2px solid grey;
  border-radius: 10px;
}

#sequenceComparatorResult {
  flex: 3;
  width: 500px;
  margin: auto 1em;
}

@media (max-width: 1024px) {
  #sequenceComparator {
    flex-direction: column;
  }

  #sequenceComparatorResult {
    width: 100%;
  }
}
</style>
