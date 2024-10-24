<script setup lang="ts">
import { RouterLink } from 'vue-router'
import BaseModal from './BaseModal.vue'

defineProps<{ open: boolean }>()
defineEmits(['close'])
</script>

<template>
  <BaseModal :open="open" title="Sign in" @close="$emit('close')">
    <form>
      <fieldset>
        <input
          name="email"
          type="email"
          placeholder="Email"
          autocomplete="email"
          required
        />
        <input
          name="password"
          type="password"
          placeholder="Password"
          required
        />
      </fieldset>
    </form>

    <p>
      Not registered? Please sign up
      <RouterLink to="/sign-up" @click="$emit('close')">here</RouterLink>.
    </p>

    <footer>
      <button>Submit</button>
    </footer>
  </BaseModal>
</template>
