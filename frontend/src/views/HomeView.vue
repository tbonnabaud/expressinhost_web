<script setup lang="ts">
import { store } from '@/lib/store'

const user = store.currentUser
</script>

<template>
  <main class="container-fluid">
    <h1 v-if="user" id="welcome">Welcome {{ user.full_name }}</h1>

    <div class="image-container">
      <img id="banner" src="../assets/banner.png" alt="Banner" />

      <p id="summary">
        <strong>ExpressInHost</strong> tunes codon sequences for recombinant
        protein expression in a host organism
      </p>
    </div>

    <p id="warning">⚠ This web application is under development ⚠</p>

    <div id="description">
      <p>
        Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed euismod
        maximus libero, vel facilisis ipsum malesuada nec. Mauris molestie ipsum
        quis risus dignissim placerat. Phasellus at suscipit erat, sit amet
        congue odio. Vivamus sed lorem vehicula, consectetur urna a, tempor
        erat. Aenean convallis gravida orci, at pharetra elit tincidunt vel.
        Integer vulputate viverra nisl sit amet facilisis. Aenean condimentum
        dapibus metus non viverra. Maecenas sagittis sem ac sollicitudin
        sollicitudin. Praesent iaculis varius sagittis. Etiam a accumsan risus.
        Praesent viverra non eros sed commodo. Quisque arcu dolor, pellentesque
        sed velit sit amet, fermentum blandit sem. Fusce blandit urna dapibus,
        sodales erat ut, euismod nisl. Sed et placerat nulla.
      </p>
    </div>
  </main>

  <footer id="contact">
    Contact: <a href="mailto:Adelaide.Raguin@hhu.de">Adélaïde Raguin</a>
  </footer>
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

#summary {
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%);
  font-size: 2.5em;
  backdrop-filter: blur(15px);
  padding: 0.5em;
  color: whitesmoke;
  text-shadow: 0 0 5px #303030;
}

#warning {
  margin-top: 1em;
  font-size: 2em;
  text-align: center;
}

#description {
  margin: 2em 5em;
  font-size: 1.5em;
  text-align: justify;
}

#contact {
  text-align: center;
  font-style: italic;
}
</style>
