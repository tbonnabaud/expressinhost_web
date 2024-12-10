<script setup lang="ts">
import { store } from '@/lib/store'
import SlideShow from '@/components/SlideShow.vue'

const user = store.currentUser
</script>

<template>
  <main class="container-fluid">
    <h1 v-if="user" id="welcome">Welcome {{ user.full_name }}</h1>


    <div class="image-container">
      <!-- <img id="banner" src="/banners/banner0.png" alt="Banner" /> -->
      <SlideShow id="banner" />

      <p id="summary">
        <strong>ExpressInHost</strong> tunes codon sequences for recombinant
        protein expression in a host organism
      </p>
    </div>

    <p id="warning">⚠ This web application is under development ⚠</p>

    <div id="description">
      <p>
        In the rapidly evolving field of recombinant protein expression,
        ExpressInHost offers novel innovative solutions for codon optimisation.
        ExpressInHost tunes mRNA sequences to enhance their expression while
        considering the challenge of proper protein folding in a host organism.
        For thousands of different organisms, and customisable ones, leveraging
        tRNA gene copy number data, it analyses and compares translation speed
        properties in native and host organisms, possibly also incorporating
        evolutionary insights and optimising mRNA structural properties for
        ribosome binding.
      </p>

      <p>
        ExpressInHost is a free, user-friendly, and fully open source web
        application that builds on our
        <a
          href="https://openresearchsoftware.metajnl.com/articles/10.5334/jors.385"
          target="_blank"
          >previous work</a
        >
        as well as on open source work from other groups for
        <a href="https://gtrnadb.ucsc.edu/" target="_blank">tRNA data</a>
        collection and
        <a
          href="https://joss.theoj.org/papers/10.21105/joss.03362"
          target="_blank"
          >ribosome binding</a
        >
        properties.
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
