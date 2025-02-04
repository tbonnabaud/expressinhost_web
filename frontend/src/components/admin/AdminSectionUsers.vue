<script setup lang="ts">
import { API } from '@/lib/api'
import { store } from '@/lib/store'
import type { User } from '@/lib/interfaces'
import { computed, onMounted, ref } from 'vue'
import AdminSectionUsersItem from './AdminSectionUsersItem.vue'
import BaseModal from '../BaseModal.vue'
import PaginationWrapper from '../PaginationWrapper.vue'

const currentUser = store.currentUser

const userList = ref([] as Array<User>)
const openDeleteModal = ref(false)
const userToRemove = ref(null as User | null)

const userListWithoutMe = computed(() =>
  userList.value.filter(user => user.id !== currentUser.value?.id),
)

onMounted(async () => await fetchUserList())

function askForDeletion(user: User) {
  openDeleteModal.value = true
  userToRemove.value = user
}

async function deleteUser() {
  if (userToRemove.value) {
    const [, error] = await API.users.deleteUser(userToRemove.value.id)

    if (error) {
      alert(
        `Fail to delete user ${userToRemove.value.full_name} (${userToRemove.value.email}).`,
      )
    }
  }

  // We close the modal
  openDeleteModal.value = false
  await fetchUserList()
}

async function fetchUserList() {
  const [data, error] = await API.users.list()

  if (!error) {
    userList.value = data
  }
}
</script>

<template>
  <BaseModal
    :open="openDeleteModal"
    title="Confirm the deletion"
    @close="openDeleteModal = false"
  >
    <p>Do you really want to delete this result?</p>

    <footer>
      <button class="secondary" @click="openDeleteModal = false">Cancel</button>
      <button class="danger" @click="deleteUser">Delete</button>
    </footer>
  </BaseModal>

  <section>
    <h2>Management of users</h2>

    <div>
      <PaginationWrapper :items="userListWithoutMe" :per-page="8">
        <template #default="{ item }">
          <AdminSectionUsersItem :user="item" :ask-delete="askForDeletion" />
        </template>
      </PaginationWrapper>
    </div>
  </section>
</template>

<style scoped>
@media (min-width: 768px) {
  #userList {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 15px;
  }
}
</style>
