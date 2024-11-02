<script setup lang="ts">
import type { Result } from '@/lib/interfaces'
import { onMounted, ref } from 'vue'
import { RouterLink } from 'vue-router'
import { API } from '@/lib/api'
import { store } from '@/lib/store'

const currentUser = store.currentUser
const resultList = ref([] as Array<Result>)

onMounted(async () => await fetchResultList())

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
  <main class="container">
    <h1>Result list</h1>

    <div class="grid-of-three">
      <article
        v-for="(result, index) in resultList"
        :key="index"
        class="zoomable"
      >
        <header>
          <RouterLink :to="'/results/' + result.id" class="contrast">
            <h4>{{ result.id }}</h4>
          </RouterLink>
        </header>

        <p>{{ result.host_codon_table_name }}</p>
      </article>
    </div>
  </main>
</template>

<style scoped>
@media (min-width: 1024px) {
  .container {
    padding-left: 15em;
    padding-right: 15em;
  }
}

.grid-of-three {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
}

.grid-of-three > article {
  margin: 0.5em;
}

.zoomable {
  transition: transform 0.2s;
}

.zoomable:hover {
  transform: scale(1.03);
}
</style>
