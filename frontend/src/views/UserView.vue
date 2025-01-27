<script setup lang="ts">
import UserProfileForm from '@/components/UserProfileForm.vue'
</script>

<template>
  <main class="container">
    <section>
      <h1>Profile</h1>
      <hr />

      <div class="centered">
        <UserProfileForm id="profileForm" />
      </div>
    </section>

    <!-- <section>
      <h1>Update password</h1>
      <hr />

      <div class="centered"></div>
    </section> -->
  </main>
</template>

<style scoped>
#profileForm {
  width: 50%;
}
</style>
