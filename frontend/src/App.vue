<script setup lang="ts">
import { ref, onMounted } from 'vue'
import { RouterLink, RouterView } from 'vue-router'
import { store } from './lib/store'
import { API, setCurrentUserInStore } from './lib/api'
import LoginForm from '@/components/LoginForm.vue'

const openLoginForm = ref(false)
const user = store.currentUser

onMounted(async () => API.auth.isLoggedIn() && (await setCurrentUserInStore()))

function logout() {
  API.auth.logout()
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
        <li v-if="user">
          <RouterLink to="/codon-tables" class="secondary">
            Codon tables
          </RouterLink>
        </li>
        <li><RouterLink to="/tuning" class="secondary">Tuning</RouterLink></li>
        <li v-if="user">
          <RouterLink to="/results" class="secondary">Results</RouterLink>
        </li>
      </ul>

      <ul>
        <li v-if="user && user.role == 'admin'">
          <RouterLink to="/admin" class="secondary">Administration</RouterLink>
        </li>
        <li><RouterLink to="/about" class="secondary">About</RouterLink></li>
        <span class="user-buttons">
          <li v-if="user">
            <RouterLink to="/user-profile">
              <button class="outline secondary">Your profile</button>
            </RouterLink>
          </li>
          <li v-if="user">
            <button class="outline secondary" @click="logout">Logout</button>
          </li>
          <li v-else>
            <button class="outline" @click="openLoginForm = true">
              Log in
            </button>
          </li>
        </span>
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

nav a.router-link-exact-active button {
  color: var(--color-text);
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

.user-buttons {
  margin-left: 1em;
}
</style>
