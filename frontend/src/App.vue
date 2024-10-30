<script setup lang="ts">
import { RouterLink, RouterView } from 'vue-router'
import LoginForm from '@/components/LoginForm.vue'
import { ref, onMounted } from 'vue'
import { store } from './lib/store'
import { API } from './lib/api'

const openLoginForm = ref(false)
const user = store.currentUser

onMounted(async () => API.users.isLoggedIn() && (await store.setCurrentUser()))

function logout() {
  API.users.logout()
  store.emptyCurrentUser()
}
</script>

<template>
  <header>
    <nav>
      <ul>
        <li><img id="logo" src="./assets/logo.png" alt="Logo" /></li>
        <li><strong>ExpressInHost</strong></li>
        <li><RouterLink to="/" class="secondary">Home</RouterLink></li>
        <!-- <li>
          <RouterLink to="/codon-tables" class="secondary"
            >Codon tables</RouterLink
          >
        </li> -->
        <li><RouterLink to="/tuning" class="secondary">Tuning</RouterLink></li>
        <!-- <li>
          <RouterLink to="/result-list" class="secondary"
            >Result list</RouterLink
          >
        </li> -->
      </ul>

      <ul>
        <li><RouterLink to="/about" class="secondary">About</RouterLink></li>
        <li v-if="user" @click="logout">
          <button class="outline secondary">Logout</button>
        </li>
        <li v-else>
          <button class="outline" @click="openLoginForm = true">Log in</button>
        </li>
      </ul>
    </nav>
  </header>

  <LoginForm :open="openLoginForm" @close="openLoginForm = false" />

  <RouterView />
</template>

<style scoped>
nav {
  padding: 0 2em;
  margin: 0;
}

nav li {
  font-size: 1.5em;
}

nav a.router-link-exact-active {
  color: var(--color-text);
  text-decoration: underline;
}

nav a.router-link-exact-active:hover {
  background-color: transparent;
}

#logo {
  max-height: 1.5em;
}
</style>
