<script setup lang="ts">
import { RouterLink } from 'vue-router'
import BaseModal from './BaseModal.vue'

defineProps<{ open: boolean }>()
defineEmits(['close'])
</script>

<template>
  <BaseModal :open="open" title="User login" @close="$emit('close')">
    <form>
      <fieldset>
        <label>
          Email
          <input
            name="email"
            type="email"
            placeholder="Email"
            autocomplete="email"
            required
          />
        </label>

        <label>
          Password
          <input
            name="password"
            type="password"
            placeholder="Password"
            required
          />
        </label>
      </fieldset>
    </form>

    <p>
      Not registered? Please register
      <RouterLink to="/sign-up" @click="$emit('close')">here</RouterLink>.
    </p>

    <footer>
      <button>Submit</button>
    </footer>
  </BaseModal>
</template>
