<script setup lang="ts">
import RequiredAdmin from '@/components/RequiredAdmin.vue'
import RequiredAuth from '@/components/RequiredAuth.vue'
import AdminWebScraping from '@/components/admin/AdminWebScraping.vue'
import { store } from '@/lib/store'

const user = store.currentUser
</script>

<template>
  <main class="container">
    <template v-if="user && user.role == 'admin'">
      <h1>Administration</h1>
      <hr />

      <AdminWebScraping />
    </template>

    <RequiredAdmin v-else-if="user && user.role != 'admin'" />

    <RequiredAuth v-else />
  </main>
</template>

<style scoped></style>
