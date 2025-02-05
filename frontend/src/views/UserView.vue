<script setup lang="ts">
import RequiredAuth from '@/components/RequiredAuth.vue'
import UserProfileForm from '@/components/UserProfileForm.vue'
import ResetPasswordForm from '@/components/ResetPasswordForm.vue'
import BaseModal from '@/components/BaseModal.vue'
import { store } from '@/lib/store'
import { API } from '@/lib/api'
import { ref } from 'vue'
import { useRouter } from 'vue-router'

const router = useRouter()
const user = store.currentUser
const openDeleteModal = ref(false)

function passwordSuccessAlert() {
  alert('Password succesfully updated!')
}

function askForDeletion() {
  openDeleteModal.value = true
}

async function deleteMe() {
  const [, error] = await API.users.deleteMe()

  // We close the modal
  openDeleteModal.value = false

  if (error) {
    alert('Fail to delete your account.')
  } else {
    API.auth.logout()
    store.emptyCurrentUser()
    router.push('/')
  }
}
</script>

<template>
  <BaseModal
    :open="openDeleteModal"
    title="Confirm the deletion"
    @close="openDeleteModal = false"
  >
    <p>Do you really want to delete your account?</p>

    <footer>
      <button class="secondary" @click="openDeleteModal = false">Cancel</button>
      <button class="danger" @click="deleteMe">Delete</button>
    </footer>
  </BaseModal>

  <main class="container">
    <template v-if="user">
      <section>
        <h2>Profile</h2>
        <hr />

        <div class="centered">
          <UserProfileForm id="profileForm" />
        </div>
      </section>

      <section>
        <h2>Update password</h2>
        <hr />

        <div class="centered">
          <ResetPasswordForm
            id="updatePasswordForm"
            @update="passwordSuccessAlert"
          />
        </div>
      </section>

      <section>
        <h2>Delete your account</h2>
        <hr />

        <div class="centered">
          <button id="askDelete" class="danger" @click="askForDeletion">
            Delete
          </button>
        </div>
      </section>
    </template>

    <RequiredAuth v-else />
  </main>
</template>

<style scoped>
#profileForm,
#updatePasswordForm,
#askDelete {
  width: 70%;
}

@media (max-width: 1024px) {
  #profileForm,
  #updatePasswordForm,
  #askDelete {
    width: 100%;
  }
}
</style>
