<script setup lang="ts">
import type { UserForm } from '@/lib/interfaces'
import { reactive, ref, watch } from 'vue'
import { useRouter } from 'vue-router'
import { API } from '@/lib/api'
import { checkPasswordConstraints } from '@/lib/checkers'
import WithAlertError from './WithAlertError.vue'

const router = useRouter()

const form = reactive({
  full_name: '',
  email: '',
  password: '',
  contact_consent: true,
} as UserForm)

const passwordError = ref('')

watch(
  () => form.password,
  password => {
    passwordError.value = checkPasswordConstraints(password)
  },
)

async function handleSubmit() {
  const [data, error] = await API.users.register(form)

  if (!error) {
    console.log(data)
    // Redirect to home
    router.push('/')
  } else if (error.code == 409) {
    alert('E-mail address already exists.')
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

      <WithAlertError :error="passwordError">
        <label
          >Password
          <i id="passwordIndications">
            (between 8 and 20 characters, with at least one letter and one
            number)
          </i>
          <input type="password" v-model="form.password" required />
        </label>
      </WithAlertError>

      <label>
        <input type="checkbox" role="switch" v-model="form.contact_consent" />
        I consent to receive e-mails regarding the application
      </label>
    </fieldset>

    <input type="submit" value="Submit" />
  </form>
</template>

<style scoped>
#passwordIndications {
  font-size: 80%;
}
</style>
