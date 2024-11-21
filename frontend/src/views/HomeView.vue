<script setup lang="ts">
import { store } from '@/lib/store'

const user = store.currentUser
</script>

<template>
  <main class="container-fluid">
    <h1 v-if="user" id="welcome">Welcome {{ user.full_name }}</h1>

    <div class="image-container">
      <img id="banner" src="../assets/banner.png" alt="Banner" />

      <p id="description">
        <strong>ExpressInHost</strong> tunes codon sequences for recombinant
        protein expression in a host organism
      </p>
    </div>

    <p id="warning">⚠ This web application is under development ⚠</p>
  </main>
</template>

<style scoped>
.container-fluid {
  padding-left: 0;
  padding-right: 0;
}

#welcome {
  margin-bottom: 1em;
  text-align: center;
}

.image-container {
  margin-top: 2em;
  position: relative;
  text-align: center;
}

#banner {
  width: 100%;
}

#description {
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  font-size: 3em;
  backdrop-filter: blur(15px);
  padding: 1em;
  color: whitesmoke;
  text-shadow: 0 0 5px #303030;
}

#warning {
  margin-top: 1em;
  font-size: 2.5em;
  text-align: center;
}
</style>
