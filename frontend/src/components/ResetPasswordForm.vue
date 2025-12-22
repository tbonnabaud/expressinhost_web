<script setup lang="ts">
import { ref, watch } from 'vue'
import { API } from '@/lib/api'
import { store } from '@/lib/store'
import { checkPasswordConstraints } from '@/lib/checkers'
import type { UserPasswordForm } from '@/lib/interfaces'
import WithAlertError from '@/components/WithAlertError.vue'

const props = defineProps<{ token?: string }>()
const emit = defineEmits(['update'])

const password = ref('')
const passwordError = ref('')

watch(password, password => {
  passwordError.value = checkPasswordConstraints(password)
})

async function handleSubmit() {
  const isLoggedUser = Boolean(store.currentUser.value)
  const token = isLoggedUser ? localStorage.getItem('accessToken') : props.token

  const [, error] = await API.users.updatePassword(
    {
      reset_token: token,
      password: password.value,
    } as UserPasswordForm,
    !isLoggedUser, // If not logged user, use reset mode
  )

  if (!error) {
    emit('update')
  }
}
</script>

<template>
  <form @submit.prevent="handleSubmit">
    <fieldset>
      <WithAlertError :error="passwordError">
        <label>
          New password
          <i id="passwordIndications">
            (between 8 and 20 characters, with at least one letter and one
            number)
          </i>
          <input type="password" v-model="password" required />
        </label>
      </WithAlertError>
    </fieldset>

    <input type="submit" value="Submit" :disabled="Boolean(passwordError)" />
  </form>
</template>

<style scoped>
#passwordIndications {
  font-size: 80%;
}
</style>
