<script setup lang="ts">
import { RouterLink, RouterView } from 'vue-router'
</script>

<template>
  <header>
    <nav>
      <ul>
        <li><strong>ExpressInHost</strong></li>
      </ul>

      <ul>
        <li><RouterLink to="/" class="secondary">Home</RouterLink></li>
        <li>
          <RouterLink to="/codon-tables" class="secondary"
            >Codon tables</RouterLink
          >
        </li>
        <li><RouterLink to="/tuning" class="secondary">Tuning</RouterLink></li>
        <li>
          <RouterLink to="/result-list" class="secondary"
            >Result list</RouterLink
          >
        </li>
        <li><RouterLink to="/about" class="secondary">About</RouterLink></li>
      </ul>
    </nav>
  </header>

  <RouterView />
</template>

<style scoped>
nav {
  margin: 0 2em;
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
</style>
