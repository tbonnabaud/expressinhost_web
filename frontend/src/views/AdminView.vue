<script setup lang="ts">
import { ref } from 'vue'
import { store } from '@/lib/store'
import type { ComponentMeta } from '@/lib/interfaces'
import RequiredAdmin from '@/components/RequiredAdmin.vue'
import RequiredAuth from '@/components/RequiredAuth.vue'
import AdminSectionWebScraping from '@/components/admin/AdminSectionWebScraping.vue'
import AdminSectionStats from '@/components/admin/AdminSectionStats.vue'
import AdminSectionUsers from '@/components/admin/AdminSectionUsers.vue'
import AdminSectionLogFile from '@/components/admin/AdminSectionLogFile.vue'

const user = store.currentUser
const currentSection = ref(null as string | null)

const sectionMapping: Record<string, ComponentMeta> = {
  scraping: {
    name: 'Codon table web scraping',
    component: AdminSectionWebScraping,
  },
  users: {
    name: 'Management of users',
    component: AdminSectionUsers,
  },
  stats: {
    name: 'Statistics',
    component: AdminSectionStats,
  },
  logfile: {
    name: 'Server log file',
    component: AdminSectionLogFile,
  },
}

function setCurrentSection(key: string) {
  if (key == currentSection.value) {
    currentSection.value = null
  } else {
    currentSection.value = key
  }
}
</script>

<template>
  <main class="container">
    <template v-if="user && user.role == 'admin'">
      <h1>Administration</h1>

      <hr />

      <div class="grid">
        <button
          v-for="(section, key) in sectionMapping"
          class="secondary"
          :class="{ contrast: currentSection == key }"
          @click="setCurrentSection(key)"
          :key="key"
        >
          {{ section.name }}
        </button>
      </div>

      <hr />

      <component
        v-if="currentSection"
        :is="sectionMapping[currentSection].component"
      />
      <p v-else>Select a section.</p>
    </template>

    <RequiredAdmin v-else-if="user && user.role != 'admin'" />

    <RequiredAuth v-else />
  </main>
</template>

<style scoped></style>
