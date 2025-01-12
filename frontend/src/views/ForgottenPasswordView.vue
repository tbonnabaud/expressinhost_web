<script setup lang="ts">
import { ref } from 'vue'
import { useRouter } from 'vue-router'
import { API } from '@/lib/api'

const router = useRouter()

const email = ref('')

async function handleSubmit() {
  const [data, error] = await API.auth.sendResetPasswordLink(email.value)

  if (!error) {
    console.log(data)
    // Redirect to home
    router.push('/')
  }
}
</script>

<template>
  <main class="container">
    <h1>Forgotten password</h1>

    <hr />

    <div class="centered">
      <form id="forgotForm" @submit.prevent="handleSubmit">
        <fieldset>
          <label>
            Please enter your e-mail address
            <input type="email" v-model="email" required />
          </label>
        </fieldset>

        <input type="submit" value="Submit" />
      </form>
    </div>
  </main>
</template>

<style scoped>
#forgotForm {
  width: 50%;
}

.centered {
  min-height: 50vh;
}
</style>
