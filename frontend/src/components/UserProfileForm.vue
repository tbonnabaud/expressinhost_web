<script setup lang="ts">
import type { UserProfileForm } from '@/lib/interfaces'
import { onMounted, reactive } from 'vue'
import { API, setCurrentUserInStore } from '@/lib/api'

const form = reactive({
  full_name: '',
  email: '',
  contact_consent: false,
} as UserProfileForm)

onMounted(async () => await getMe())

async function getMe() {
  const [data, error] = await API.users.me()

  if (!error) {
    form.full_name = data.full_name
    form.email = data.email
    form.contact_consent = data.contact_consent
  }
}

async function handleSubmit() {
  const [, error] = await API.users.updateProfile(form)

  if (!error) {
    alert('Profile edited succesfully!')
    await setCurrentUserInStore()
  } else if (error.code == 409) {
    alert('E-mail address already exists.')
  }
}
</script>

<template>
  <form @submit.prevent="handleSubmit">
    <fieldset>
      <label>
        Full name
        <input type="text" v-model="form.full_name" required />
      </label>

      <label>
        E-mail
        <input type="email" v-model="form.email" required />
      </label>

      <label>
        <input type="checkbox" v-model="form.contact_consent" role="switch" />
        I consent to receive e-mails regarding the application
      </label>
    </fieldset>

    <input type="submit" value="Update" />
  </form>
</template>
