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
#resetForm {
  width: 50%;
}

.centered {
  min-height: 50vh;
}
</style>
