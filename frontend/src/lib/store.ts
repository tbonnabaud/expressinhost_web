import { ref } from 'vue'
import type { User } from './interfaces'

export const store = {
  currentUser: ref(null as User | null),
  emptyCurrentUser() {
    this.currentUser.value = null
  },
  setCurrentUser(user: User) {
    this.currentUser.value = user
  },
}
