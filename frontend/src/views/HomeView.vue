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
        <strong>ExpressInHost_web</strong> tunes genes<br />for recombinant protein
        expression
      </p>
    </div>

    <div id="description">
      <p>
        In the rapidly evolving field of recombinant protein expression,
        <strong>ExpressInHost_web</strong> offers novel and innovative solutions for
        codon optimisation. <strong>ExpressInHost_web</strong> tunes mRNA sequences
        to enhance their expression while considering the challenge of proper
        protein folding in a host organism. Utilising tRNA gene copy number data
        from thousands of organisms, including customisable ones, it analyses
        and compares translational speed profiles between native and host
        systems. <strong>ExpressInHost_web</strong> also incorporates evolutionary
        perspectives and optimises mRNA structural features for efficient
        ribosome binding.
      </p>

      <p>
        <strong>ExpressInHost_web</strong> is a free, user-friendly, and fully
        <a href="https://github.com/tbonnabaud/expressinhost_web"
          >open source</a
        >
        web application that builds on our
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
    <p>Contact: <a href="mailto:Adelaide.Raguin@hhu.de">Adélaïde Raguin</a></p>

    <p>
      🅭 All results obtained from our softwares, e.g. by using this web
      interface, are licensed under the
      <a href="https://creativecommons.org/licenses/by/4.0/legalcode.txt"
        >CC BY 4.0 terms</a
      >, including that attribution must be given to the creator, the
      ExpressInHost_web team, by citing <a href="">our release article</a> and
      <a
        href="https://openresearchsoftware.metajnl.com/articles/10.5334/jors.385"
        >our previous work</a
      >.
    </p>
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
  font-size: clamp(1.25rem, 4vw, 2.5rem);
  backdrop-filter: blur(15px);
  padding: clamp(0.5rem, 2vw, 1rem);
  color: whitesmoke;
  text-shadow: 0 0 5px #303030;
  max-width: 90%;
}

#description {
  margin: clamp(1rem, 3vw, 2rem) clamp(1rem, 5vw, 5rem);
  font-size: clamp(1rem, 2vw, 1.5rem);
  text-align: justify;
}

#contact {
  text-align: center;
  font-style: italic;
}

@media (max-width: 768px) {
  #description {
    text-align: left;
  }
}

footer {
  margin: clamp(1rem, 3vw, 2rem) clamp(1rem, 5vw, 5rem);
}
</style>
