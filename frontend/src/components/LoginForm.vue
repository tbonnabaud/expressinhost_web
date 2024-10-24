<script setup lang="ts">
import { RouterLink } from 'vue-router'
import BaseModal from './BaseModal.vue'
import { reactive, useTemplateRef } from 'vue'

defineProps<{ open: boolean }>()
const emit = defineEmits(['close', 'login'])

const form = reactive({
  email: '',
  password: '',
})
const formRef = useTemplateRef('login-form')

async function handleSubmit() {
  if (formRef.value?.checkValidity()) {
    console.log(JSON.stringify(form))
    emit('login')
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
            v-model="form.email"
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
