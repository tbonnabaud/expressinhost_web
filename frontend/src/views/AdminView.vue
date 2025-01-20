<script setup lang="ts">
import { store } from '@/lib/store'
import RequiredAdmin from '@/components/RequiredAdmin.vue'
import RequiredAuth from '@/components/RequiredAuth.vue'
import AdminWebScraping from '@/components/admin/AdminWebScraping.vue'
import AdminUsers from '@/components/admin/AdminUsers.vue'
import AdminDowloadLogFile from '@/components/admin/AdminDowloadLogFile.vue'

const user = store.currentUser
</script>

<template>
  <main class="container">
    <template v-if="user && user.role == 'admin'">
      <h1>Administration</h1>
      <hr />

      <section>
        <AdminWebScraping />
      </section>

      <hr />

      <section>
        <h2>Users</h2>
        <AdminUsers />
      </section>

      <hr />

      <section>
        <h2>Log files</h2>
        <AdminDowloadLogFile />
      </section>
    </template>

    <RequiredAdmin v-else-if="user && user.role != 'admin'" />

    <RequiredAuth v-else />
  </main>
</template>

<style scoped></style>
