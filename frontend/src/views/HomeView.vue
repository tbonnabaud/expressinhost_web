<script setup lang="ts">
import { store } from '@/lib/store'
import BannerSlideshow from '@/components/BannerSlideshow.vue'

const user = store.currentUser
</script>

<template>
  <main class="container-fluid">
    <h1 v-if="user" id="welcome">Welcome {{ user.full_name }}</h1>

    <div class="image-container">
      <BannerSlideshow id="banner" />

      <p id="summary">
        <strong>ExpressInHost</strong> tunes genes<br />for recombinant protein
        expression
      </p>
    </div>

    <div id="description">
      <p>
        In the rapidly evolving field of recombinant protein expression,
        <strong>ExpressInHost</strong> offers novel and innovative solutions for
        codon optimisation. <strong>ExpressInHost</strong> tunes mRNA sequences
        to enhance their expression while considering the challenge of proper
        protein folding in a host organism. Utilising tRNA gene copy number data
        from thousands of organisms, including customisable ones, it analyses
        and compares translational speed profiles between native and host
        systems. <strong>ExpressInHost</strong> also incorporates evolutionary
        perspectives and optimises mRNA structural features for efficient
        ribosome binding.
      </p>

      <p>
        <strong>ExpressInHost</strong> is a free, user-friendly, and fully open
        source web application that builds on our
        <a
          href="https://openresearchsoftware.metajnl.com/articles/10.5334/jors.385"
          target="_blank"
          >previous work</a
        >
        as well as on open source results from other groups for
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

#description {
  margin: 2em 5em;
  font-size: 1.5em;
  text-align: justify;
}

#contact {
  text-align: center;
  font-style: italic;
}

@media (max-width: 1024px) {
  #summary {
    font-size: 2em;
  }
}

@media (max-width: 768px) {
  #summary {
    font-size: 1.5em;
  }

  #description {
    font-size: 1em;
    text-align: center;
  }
}
</style>
