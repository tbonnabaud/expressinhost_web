<script setup lang="ts">
import { ref } from 'vue'

defineProps<{
  error?: string
  errors?: Array<string>
}>()

const show = ref(true)
</script>

<template>
  <div>
    <slot></slot>

    <Transition name="fade">
      <div
        v-if="error || errors"
        class="alert alert-danger"
        :class="{ 'alert-hidden': !show }"
      >
        <a class="close" aria-label="Close" @click="show = false">
          <span aria-hidden="true">&times;</span>
        </a>
        <span v-if="error">
          {{ error }}
        </span>

        <ul v-else-if="errors">
          <li v-for="(err, index) in errors" :key="index">{{ err }}</li>
        </ul>
      </div>
    </Transition>
  </div>
</template>

<style scoped>
.fade-enter-active,
.fade-leave-active {
  transition: opacity 0.5s ease;
}

.fade-enter-from,
.fade-leave-to {
  opacity: 0;
}

.alert {
  position: relative;
  padding: 0.75rem 1.25rem;
  margin-right: 15px;
  margin-bottom: 1rem;
  border: 1px solid transparent;
  border-radius: 0.25rem;
}

.alert-danger {
  color: #721c24;
  background-color: #f8d7da;
  border-color: #721c24;
}

.alert-danger li {
  color: #721c24;
}

.close {
  float: right;
  text-decoration: none;
  cursor: pointer;
  color: #721c24;
}

.alert-hidden {
  display: none;
}
</style>
