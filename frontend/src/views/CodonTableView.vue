<script setup lang="ts">
import CodonTableForm from '@/components/codon-tables/CodonTableForm.vue'
import RequiredAuth from '@/components/RequiredAuth.vue'
import { store } from '@/lib/store'

const user = store.currentUser
</script>

<template>
  <main class="container-fluid">
    <template v-if="user">
      <h1>tRNA GCN tables</h1>
      <hr />

      <CodonTableForm />
    </template>

    <RequiredAuth v-else />
  </main>
</template>
