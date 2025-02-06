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
.centered {
  min-height: 50vh;
}

@media (min-width: 1024px) {
  #register {
    width: 50%;
  }
}
</style>
