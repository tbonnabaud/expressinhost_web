<script setup lang="ts">
import { RouterLink, useRouter } from 'vue-router'
import BaseModal from './BaseModal.vue'
import { reactive, useTemplateRef } from 'vue'
import { API, setCurrentUserInStore } from '@/lib/api'

defineProps<{ open: boolean }>()
const emit = defineEmits(['close'])

const router = useRouter()

const form = reactive({
  username: '',
  password: '',
})
const formRef = useTemplateRef('login-form')

async function handleSubmit() {
  if (formRef.value?.checkValidity()) {
    const error = await API.users.login(form)
    if (error === null) {
      await setCurrentUserInStore()
      emit('close')
      // Redirect to home
      router.push('/')
    }
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

    <footer>
      <button @click="handleSubmit">Submit</button>
    </footer>
  </BaseModal>
</template>
