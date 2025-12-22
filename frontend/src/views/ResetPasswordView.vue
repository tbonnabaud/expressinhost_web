<script setup lang="ts">
import { useRouter } from 'vue-router'
import ResetPasswordForm from '@/components/ResetPasswordForm.vue'

defineProps<{ token: string }>()

const router = useRouter()

function redirectToHome() {
  router.push('/')
}
</script>

<template>
  <main class="container">
    <h1>Reset password</h1>

    <hr />

    <div class="centered">
      <ResetPasswordForm
        id="resetForm"
        :token="token"
        @update="redirectToHome"
      />
    </div>
  </main>
</template>

<style scoped>
.centered {
  min-height: 50vh;
}

@media (min-width: 1024px) {
  #resetForm {
    width: 50%;
  }
}
</style>
