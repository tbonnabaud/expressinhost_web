<script setup lang="ts">
import { ref, watch } from 'vue'
import { useRouter } from 'vue-router'
import { API } from '@/lib/api'
import { checkPasswordConstraints } from '@/lib/checkers'
import type { UserPasswordForm } from '@/lib/interfaces'
import WithAlertError from '@/components/WithAlertError.vue'

const props = defineProps<{ token: string }>()

const router = useRouter()

const password = ref('')
const passwordError = ref('')

watch(password, password => {
  passwordError.value = checkPasswordConstraints(password)
})

async function handleSubmit() {
  const [, error] = await API.users.updatePassword({
    reset_token: props.token,
    password: password.value,
  } as UserPasswordForm)

  if (!error) {
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

#passwordIndications {
  font-size: 80%;
}
</style>
