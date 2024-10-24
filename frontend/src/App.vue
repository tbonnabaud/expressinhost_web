<script setup lang="ts">
import { RouterLink, RouterView } from 'vue-router'
import LoginForm from '@/components/LoginForm.vue'
import { ref } from 'vue'

const openLoginForm = ref(false)
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
        <li><RouterLink to="/about" class="secondary">About</RouterLink></li>
      </ul>

      <ul>
        <li>
          <button class="outline secondary" @click="openLoginForm = true">
            Sign in
          </button>
        </li>
        <!-- <li><RouterLink to="/sign-up" class="secondary">Sign up</RouterLink></li> -->
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
