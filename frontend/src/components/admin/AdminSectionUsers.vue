<script setup lang="ts">
import { API } from '@/lib/api'
import { store } from '@/lib/store'
import type { User } from '@/lib/interfaces'
import { computed, onMounted, ref } from 'vue'
import AdminSectionUsersItem from './AdminSectionUsersItem.vue'

const currentUser = store.currentUser

const userList = ref([] as Array<User>)

onMounted(async () => await fetchUserList())

const userListWithoutMe = computed(() =>
  userList.value.filter(user => user.id !== currentUser.value?.id),
)

async function fetchUserList() {
  const [data, error] = await API.users.list()

  if (!error) {
    userList.value = data
  }
}
</script>

<template>
  <section>
    <h2>Management of users</h2>

    <div id="userList">
      <AdminSectionUsersItem
        v-for="user in userListWithoutMe"
        :user="user"
        :key="user.id"
      />
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
