import { ref } from 'vue'
import { API } from './api'
import type { User } from './interfaces'

export const store = {
  currentUser: ref(null as User | null),
  async setCurrentUser() {
    const [data, error] = await API.users.me()

    if (!error) {
      this.currentUser.value = data
    } else {
      this.currentUser.value = null
    }
  },
  emptyCurrentUser() {
    this.currentUser.value = null
  },
}
