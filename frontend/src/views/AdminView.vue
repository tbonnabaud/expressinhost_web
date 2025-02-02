<script setup lang="ts">
import { ref, type Component } from 'vue'
import { store } from '@/lib/store'
import RequiredAdmin from '@/components/RequiredAdmin.vue'
import RequiredAuth from '@/components/RequiredAuth.vue'
import AdminSectionWebScraping from '@/components/admin/AdminSectionWebScraping.vue'
import AdminSectionStats from '@/components/admin/AdminSectionStats.vue'
import AdminSectionUsers from '@/components/admin/AdminSectionUsers.vue'
import AdminSectionLogFile from '@/components/admin/AdminSectionLogFile.vue'

const user = store.currentUser
const currentSection = ref(null as string | null)
const sectionMapping: Record<string, Component> = {
  scraping: AdminSectionWebScraping,
  users: AdminSectionUsers,
  stats: AdminSectionStats,
  logfile: AdminSectionLogFile,
}
</script>

<template>
  <main class="container">
    <template v-if="user && user.role == 'admin'">
      <h1>Administration</h1>
      <hr />

      <div class="grid">
        <button
          class="secondary"
          :class="{ contrast: currentSection == 'scraping' }"
          @click="currentSection = 'scraping'"
        >
          Codon table web scraping
        </button>
        <button
          class="secondary"
          :class="{ contrast: currentSection == 'users' }"
          @click="currentSection = 'users'"
        >
          Management of users
        </button>
        <button
          class="secondary"
          :class="{ contrast: currentSection == 'stats' }"
          @click="currentSection = 'stats'"
        >
          Statistics
        </button>
        <button
          class="secondary"
          :class="{ contrast: currentSection == 'logfile' }"
          @click="currentSection = 'logfile'"
        >
          Log file
        </button>
      </div>

      <hr />

      <component v-if="currentSection" :is="sectionMapping[currentSection]" />
      <p v-else>Select a section.</p>
    </template>

    <RequiredAdmin v-else-if="user && user.role != 'admin'" />

    <RequiredAuth v-else />
  </main>
</template>

<style scoped></style>
