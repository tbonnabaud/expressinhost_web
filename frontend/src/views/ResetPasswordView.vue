<script setup lang="ts">
import { ref } from 'vue'
import { useRouter } from 'vue-router'
import { API } from '@/lib/api'
import type { UserPasswordForm } from '@/lib/interfaces'

const props = defineProps<{ token: string }>()

const router = useRouter()

const password = ref('')

async function handleSubmit() {
  const [data, error] = await API.users.updatePassword({
    reset_token: props.token,
    password: password.value,
  } as UserPasswordForm)

  if (!error) {
    console.log(data)
    // Redirect to home
    router.push('/')
  }
}
</script>

<template>
  <main class="container">
    <h1>Reset password</h1>

    <hr />

    <div class="centered">
      <form id="resetForm" @submit.prevent="handleSubmit">
        <fieldset>
          <label>
            New password
            <input type="password" v-model="password" required />
          </label>
        </fieldset>

        <input type="submit" value="Submit" />
      </form>
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
