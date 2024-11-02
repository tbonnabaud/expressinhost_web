<script setup lang="ts">
import { API } from '@/lib/api'
import { onMounted, ref, watch } from 'vue'
import ResultContent from '@/components/results/ResultContent.vue'
import type { Result, TunedSequence } from '@/lib/interfaces'

const props = defineProps<{ id: string }>()

const result = ref({} as Result)
const tunedSequences = ref([] as Array<TunedSequence>)

onMounted(async () => {
  result.value = await fetchResult(props.id)
})

watch(result, async value => {
  if (value.id) {
    tunedSequences.value = await fetchTunedSequences(value.id)
  }
})

async function fetchResult(id: string) {
  const [data, error] = await API.results.get(id)

  return !error ? data : {}
}

async function fetchTunedSequences(resultId: string) {
  const [data, error] = await API.tunedSequences.list(resultId)

  return !error ? data : []
}
</script>

<template>
  <main class="container">
    <ResultContent :result="result" :tuned_sequences="tunedSequences" />
  </main>
</template>

<style scoped>
@media (min-width: 1024px) {
  .container {
    padding-left: 15em;
    padding-right: 15em;
  }
}
</style>
