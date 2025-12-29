<script setup lang="ts">
import { RouterLink } from 'vue-router'
import BaseModal from './BaseModal.vue'
import { reactive, useTemplateRef } from 'vue'
import { API } from '@/lib/api'

defineProps<{ open: boolean }>()
const emit = defineEmits(['close'])

const form = reactive({
  username: '',
  password: '',
})
const formRef = useTemplateRef('login-form')

async function handleSubmit() {
  if (formRef.value?.checkValidity()) {
    // Login submits a traditional form that causes page reload
    // The server will set cookie and redirect to home page
    await API.auth.login(form)
    // Note: If we reach here, login failed (otherwise page would have reloaded)
    // Error handling is done via server redirect or page reload
  } else {
    formRef.value?.reportValidity()
  }
}
</script>

<template>
  <BaseModal :open="open" title="User login" @close="$emit('close')">
    <form @keyup.enter="handleSubmit" ref="login-form">
      <fieldset>
        <label>
          Email
          <input
            name="email"
            type="email"
            placeholder="Email"
            autocomplete="email"
            v-model="form.username"
            required
          />
        </label>

        <label>
          Password
          <input
            name="password"
            type="password"
            placeholder="Password"
            v-model="form.password"
            required
          />
        </label>
      </fieldset>
    </form>

    <p>
      Not registered? Please register
      <RouterLink to="/register" @click="$emit('close')">here</RouterLink>.
    </p>

    <p>
      Forgot your password? Please click
      <RouterLink to="/forgotten-password" @click="$emit('close')">
        here</RouterLink
      >.
    </p>

    <footer>
      <button @click="handleSubmit">Submit</button>
    </footer>
  </BaseModal>
</template>
