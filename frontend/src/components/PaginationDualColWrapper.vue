<script setup lang="ts" generic="T extends { id: string }">
import { ref, computed } from 'vue'

const props = defineProps<{
  items: T[]
  perPage?: number
}>()

const perPage = props.perPage ?? 5
const currentPage = ref(1)
const totalPages = computed(() => Math.ceil(props.items.length / perPage))

const paginatedItems = computed(() => {
  const start = (currentPage.value - 1) * perPage
  return props.items.slice(start, start + perPage)
})

const nextPage = () => {
  if (currentPage.value < totalPages.value) {
    currentPage.value++
  }
}

const prevPage = () => {
  if (currentPage.value > 1) {
    currentPage.value--
  }
}
</script>

<template>
  <div class="pagination">
    <button class="secondary" @click="prevPage" :disabled="currentPage === 1">
      Prev
    </button>
    <span>Page {{ currentPage }} of {{ totalPages }}</span>
    <button
      class="secondary"
      @click="nextPage"
      :disabled="currentPage === totalPages"
    >
      Next
    </button>
  </div>

  <div id="content">
    <div v-for="item in paginatedItems" :key="item.id">
      <slot :item="item"></slot>
    </div>
  </div>

  <div class="pagination">
    <button class="secondary" @click="prevPage" :disabled="currentPage === 1">
      Prev
    </button>
    <span>Page {{ currentPage }} of {{ totalPages }}</span>
    <button
      class="secondary"
      @click="nextPage"
      :disabled="currentPage === totalPages"
    >
      Next
    </button>
  </div>
</template>

<style scoped>
@media (min-width: 768px) {
  #content {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 15px;
  }
}

#content {
  margin-bottom: 20px;
}

.pagination {
  display: flex;
  justify-content: center;
  align-items: center;
  gap: 10px;
  margin-bottom: 20px;
}
</style>
