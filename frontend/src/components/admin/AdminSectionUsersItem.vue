<script setup lang="ts">
import type { User } from '@/lib/interfaces'
import { API } from '@/lib/api'
import { ref } from 'vue'

const props = defineProps<{ user: User }>()
defineEmits<{
  askDelete: [user: User]
}>()

const userRole = ref(props.user.role)
const ariaInvalid = ref(undefined as boolean | undefined)

async function handleUpdate() {
  const [, error] = await API.users.updateUserRole(props.user.id, {
    role: userRole.value,
  })

  if (!error) {
    ariaInvalid.value = false
  } else {
    ariaInvalid.value = true
    alert('Failed to update the user role.')
  }

  await new Promise(resolve => setTimeout(resolve, 5000))
  ariaInvalid.value = undefined
}
</script>

<template>
  <article>
    <header>
      {{ user.full_name }}
      (<a :href="'mailto:' + user.email">{{ user.email }}</a
      >)
    </header>

    <div class="actions">
      <select
        class="role-select"
        v-model="userRole"
        @change="handleUpdate"
        :aria-invalid="ariaInvalid"
      >
        <option value="member">member</option>
        <option value="admin">admin</option>
      </select>

      <button class="danger" @click="$emit('askDelete', user)">
        Delete
      </button>
    </div>
  </article>
</template>

<style scoped>
@media (min-width: 768px) {
  .actions {
    padding-bottom: 0;
    display: grid;
    grid-template-columns: 2fr 1fr;
    grid-gap: 15px;
  }

  p {
    line-height: 3em;
  }
}

button,
select {
  height: fit-content;
  width: 100%;
}

header {
  font-size: large;
}
</style>
