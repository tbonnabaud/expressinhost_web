<script setup lang="ts">
import RegisterForm from '@/components/RegisterForm.vue'
</script>

<template>
  <main class="container">
    <h1>Register</h1>

    <hr />

    <div class="centered">
      <RegisterForm id="register" />
    </div>
  </main>
</template>

<style scoped>
@media (min-width: 1024px) {
  .container {
    padding-left: 15em;
    padding-right: 15em;
  }
}

#register {
  width: 50%;
}

.centered {
  min-height: 50vh;
}
</style>
