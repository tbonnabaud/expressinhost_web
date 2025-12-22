<script setup lang="ts">
import { onMounted, ref, watch } from 'vue'
import type { ResultWithId } from '@/lib/interfaces'
import { API } from '@/lib/api'
import { store } from '@/lib/store'
import ResultListItem from './ResultListItem.vue'
import BaseModal from '../BaseModal.vue'
import PaginationDualColWrapper from '../PaginationDualColWrapper.vue'

const currentUser = store.currentUser
const resultList = ref([] as Array<ResultWithId>)
const openDeleteAllModal = ref(false)

onMounted(async () => await fetchResultList())

// Watch needed because query depends on a value in the store
watch(currentUser, async () => await fetchResultList())

async function fetchResultList() {
  if (currentUser.value) {
    const [data, error] = await API.results.list()

    if (!error) {
      resultList.value = data
    }
  }
}

async function deleteAllResults() {
  if (currentUser.value) {
    const [, error] = await API.results.deleteAll()

    if (!error) {
      resultList.value = []
    }
  }

  openDeleteAllModal.value = false
}
</script>

<template>
  <BaseModal
    :open="openDeleteAllModal"
    title="Confirm the deletion"
    @close="openDeleteAllModal = false"
  >
    <p>Do you really want to delete all your results?</p>

    <footer>
      <button class="secondary" @click="openDeleteAllModal = false">
        Cancel
      </button>
      <button class="danger" @click="deleteAllResults">Delete all</button>
    </footer>
  </BaseModal>

  <div class="grid">
    <h3 id="resultTotal">
      Total: <ins>{{ resultList.length }}</ins>
    </h3>

    <button
      id="deleteAllResults"
      type="button"
      class="danger"
      @click="openDeleteAllModal = true"
    >
      Delete all results
    </button>
  </div>

  <PaginationDualColWrapper
    id="resultPaginationTop"
    :items="resultList"
    :per-page="20"
  >
    <template v-slot="{ item }">
      <ResultListItem :result="item" />
    </template>
  </PaginationDualColWrapper>
</template>

<style scoped>
#resultPaginationTop {
  margin-top: auto;
}

#deleteAllResults {
  margin-left: auto;
  height: fit-content;
}
</style>
