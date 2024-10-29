<script setup lang="ts">
import type { UserForm } from '@/lib/interfaces'
import { reactive } from 'vue'
import { useRouter } from 'vue-router'
import { API } from '@/lib/api'

const router = useRouter()

const form = reactive({
  full_name: '',
  email: '',
  password: '',
} as UserForm)

async function handleSubmit() {
  const [data, error] = await API.users.register(form)

  if (!error) {
    console.log(data)
    // Redirect to home
    router.push('/')
  }
}
</script>

<template>
  <form @submit.prevent="handleSubmit">
    <fieldset>
      <label
        >Full name
        <input type="text" v-model="form.full_name" required />
      </label>

      <label
        >Email
        <input type="email" v-model="form.email" required />
      </label>

      <label
        >Password
        <input type="password" v-model="form.password" required />
      </label>
    </fieldset>

    <input type="submit" value="Submit" />
  </form>
</template>
