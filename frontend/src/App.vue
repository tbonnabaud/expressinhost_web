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
  // Try to load user from cookie-based authentication
  await setCurrentUserInStore()
  window.addEventListener('resize', handleResize)
})

onBeforeUnmount(() => {
  window.removeEventListener('resize', handleResize)
})

function handleResize() {
  openMenu.value = window.innerWidth > 768
}

async function logout() {
  await API.auth.logout()
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
            >tRNA GCN tables</RouterLink
          >
        </li>
        <li><RouterLink to="/tuning" class="secondary">Tuning</RouterLink></li>
        <li v-if="user">
          <RouterLink to="/results" class="secondary">Results</RouterLink>
        </li>
        <li>
          <RouterLink to="/sequence-comparator" class="secondary">
            Sequence comparator
          </RouterLink>
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
  gap: clamp(0.5rem, 1.5vw, 1.5rem);
  list-style: none;
  padding: 0;
  margin: 0;
}

nav li {
  font-size: clamp(0.9rem, 1.5vw, 1.25rem);
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
  font-size: 2em;
  padding: 0.25em;
  user-select: none;
  transition: transform 0.2s ease;
}

.burger-menu:hover {
  transform: scale(1.1);
}

.burger-menu:active {
  transform: scale(0.95);
}

/* Intermediate screens (tablets, small laptops 1024-1366px) */
@media (min-width: 769px) and (max-width: 1366px) {
  nav {
    padding: 0 1em;
  }

  nav ul {
    gap: clamp(0.3rem, 1vw, 0.8rem);
  }

  nav li {
    font-size: clamp(0.85rem, 1.2vw, 1rem);
  }

  .nav-links {
    margin-left: 1em;
  }

  .user-buttons {
    gap: 0.3rem;
  }

  .user-buttons button {
    padding: 0.5rem 0.75rem;
    font-size: 0.9rem;
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
    flex-direction: column;
    width: 100%;
    text-align: center;
  }

  .user-buttons {
    width: 100%;
    justify-content: center;
  }

  .nav-links li,
  .user-buttons li {
    width: 100%;
    padding: 0.5rem 0;
  }

  .nav-links a,
  .user-buttons button {
    width: 100%;
  }
}
</style>
