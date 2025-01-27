<script setup lang="ts">
import RequiredAuth from '@/components/RequiredAuth.vue'
import UserProfileForm from '@/components/UserProfileForm.vue'
import ResetPasswordForm from '@/components/ResetPasswordForm.vue'
import { store } from '@/lib/store'

const user = store.currentUser

function passwordSuccessAlert() {
  alert('Password succesfully updated!')
}
</script>

<template>
  <main class="container">
    <template v-if="user">
      <section>
        <h1>Profile</h1>
        <hr />

        <div class="centered">
          <UserProfileForm id="profileForm" />
        </div>
      </section>

      <section>
        <h1>Update password</h1>
        <hr />

        <div class="centered">
          <ResetPasswordForm
            id="updatePasswordForm"
            @update="passwordSuccessAlert"
          />
        </div>
      </section>
    </template>

    <RequiredAuth v-else />
  </main>
</template>

<style scoped>
#profileForm,
#updatePasswordForm {
  width: 50%;
}
</style>
