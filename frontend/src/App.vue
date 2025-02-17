<script setup lang="ts">
import { ref, onMounted, onBeforeUnmount } from 'vue'
import { RouterLink, RouterView, useRouter } from 'vue-router'
import { store } from './lib/store'
import { API, setCurrentUserInStore } from './lib/api'
import LoginForm from '@/components/LoginForm.vue'

const router = useRouter()
const user = store.currentUser
const openLoginForm = ref(false)
const openMenu = ref(window.innerWidth > 768)

onMounted(async () => {
  API.auth.isLoggedIn()
  await setCurrentUserInStore()
  window.addEventListener('resize', handleResize)
})

onBeforeUnmount(() => {
  window.removeEventListener('resize', handleResize)
})

function handleResize() {
  openMenu.value = window.innerWidth > 768
}

function logout() {
  API.auth.logout()
  store.emptyCurrentUser()
  router.push('/')
}
</script>

<template>
  <header>
    <nav>
      <ul id="navBrand">
        <li><img id="logo" src="./assets/cropped_logo.png" alt="Logo" /></li>
        <li><strong>ExpressInHost</strong></li>
        <li class="burger-menu" @click="openMenu = !openMenu">☰</li>
      </ul>

      <ul v-show="openMenu" class="nav-links">
        <li><RouterLink to="/" class="secondary">Home</RouterLink></li>
        <li v-if="user">
          <RouterLink to="/codon-tables" class="secondary"
            >Codon tables</RouterLink
          >
        </li>
        <li><RouterLink to="/tuning" class="secondary">Tuning</RouterLink></li>
        <li v-if="user">
          <RouterLink to="/results" class="secondary">Results</RouterLink>
        </li>
        <li v-if="user && user.role == 'admin'">
          <RouterLink to="/admin" class="secondary">Administration</RouterLink>
        </li>
        <li><RouterLink to="/about" class="secondary">About</RouterLink></li>
      </ul>

      <ul v-show="openMenu" class="user-buttons">
        <li v-if="user">
          <RouterLink to="/user-profile">
            <button class="outline secondary">Your profile</button>
          </RouterLink>
        </li>
        <li v-if="user">
          <button class="outline secondary" @click="logout">Logout</button>
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
  display: flex;
  flex-wrap: wrap;
  justify-content: space-between;
  align-items: center;
  padding: 0 2em;
}

nav ul {
  display: flex;
  align-items: center;
  gap: 1em;
  list-style: none;
  padding: 0;
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

.nav-links {
  flex-grow: 1;
  display: flex;
  margin-left: 1.5em;
  flex-wrap: wrap;
}

.user-buttons {
  display: flex;
  gap: 0.5em;
}

#logo {
  height: 2.5em;
  border-radius: 10px;
}

.burger-menu {
  display: none;
  cursor: pointer;
  font-size: 1.5em;
}

@media (max-width: 1024px) {
  nav {
    padding: 0 1em;
  }

  nav ul {
    gap: 0;
  }

  nav li {
    font-size: 1em;
  }

  .nav-links {
    margin-left: auto;
  }
}

@media (max-width: 768px) {
  nav {
    flex-direction: column;
    align-items: center;
  }

  .burger-menu {
    display: block;
  }

  #navBrand li {
    font-size: 1.5em;
  }

  .nav-links {
    margin-left: 0;
    justify-content: center;
  }
}
</style>
