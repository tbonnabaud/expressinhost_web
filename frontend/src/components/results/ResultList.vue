<script setup lang="ts">
import { computed, onMounted, ref, watch } from 'vue'
import type { Result } from '@/lib/interfaces'
import { API } from '@/lib/api'
import { store } from '@/lib/store'
import ResultListItem from './ResultListItem.vue'
import ResultListPagination from './ResultListPagination.vue'

const currentUser = store.currentUser
const resultList = ref([] as Array<Result>)

// Pagination
const totalResultCount = ref(0)
const limitPerPage = ref(8)
const currentPageNumber = ref(1)

const numberOfPages = computed(() =>
  Math.ceil(totalResultCount.value / limitPerPage.value),
)

onMounted(async () => {
  await fetchTotalResultCount()
  await fetchResultList()
})

// Watch needed because query depends on a value in the store
watch(currentUser, async () => await fetchResultList())
watch(currentPageNumber, async () => await fetchResultList())

async function fetchResultList() {
  if (currentUser.value) {
    const [data, error] = await API.results.list(
      limitPerPage.value,
      (currentPageNumber.value - 1) * limitPerPage.value,
    )

    if (!error) {
      resultList.value = data
    }
  }
}

async function fetchTotalResultCount() {
  if (currentUser.value) {
    const [data, error] = await API.results.count()

    if (!error) {
      totalResultCount.value = data
    }
  }
}
</script>

<template>
  <h3 id="resultTotal">
    Total: <ins>{{ totalResultCount }}</ins>
  </h3>

  <ResultListPagination
    id="resultPaginationTop"
    v-model="currentPageNumber"
    :total="numberOfPages"
  />

  <div class="grid-of-three">
    <ResultListItem
      v-for="(result, index) in resultList"
      :result="result"
      :key="index"
    />
  </div>

  <ResultListPagination
    id="resultPaginationBottom"
    v-model="currentPageNumber"
    :total="numberOfPages"
  />
</template>

<style scoped>
@media (min-width: 1024px) {
  .grid-of-three {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
  }
}

#resultTotal {
  text-align: center;
}
</style>
