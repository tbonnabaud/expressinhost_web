<script setup lang="ts">
import { onMounted, ref, watch } from 'vue'
import type { Result } from '@/lib/interfaces'
import { API } from '@/lib/api'
import { store } from '@/lib/store'
import ResultListItem from './ResultListItem.vue'

const currentUser = store.currentUser
const resultList = ref([] as Array<Result>)

onMounted(async () => await fetchResultList())

// Watch needed because query depends on a value in the store
watch(currentUser, async () => await fetchResultList())

async function fetchResultList() {
  if (currentUser.value) {
    const [data, error] = await API.results.list(currentUser.value.id)

    if (!error) {
      resultList.value = data
    }
  }
}
</script>

<template>
  <div class="grid-of-three">
    <ResultListItem
      v-for="(result, index) in resultList"
      :result="result"
      :key="index"
    />
  </div>
</template>

<style scoped>
@media (min-width: 1024px) {
  .grid-of-three {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
  }
}
</style>
