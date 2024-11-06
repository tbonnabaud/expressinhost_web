<script setup lang="ts">
import { onMounted, ref, watch } from 'vue'
import { API } from '@/lib/api'
import type { Result, TunedSequence } from '@/lib/interfaces'
import { store } from '@/lib/store'
import ResultContent from '@/components/results/ResultContent.vue'
import RequiredAuth from '@/components/RequiredAuth.vue'

const props = defineProps<{ id: string }>()

const user = store.currentUser

const result = ref({} as Result)
const tunedSequences = ref([] as Array<TunedSequence>)

onMounted(async () => {
  if (user) {
    result.value = await fetchResult(props.id)
  }
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
    <ResultContent
      v-if="user"
      :result="result"
      :tuned_sequences="tunedSequences"
    />
    <RequiredAuth v-else />
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
