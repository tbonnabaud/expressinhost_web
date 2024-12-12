<script setup lang="ts">
import { ref, watch } from 'vue'

const props = defineProps<{
  error?: string | null
  errors?: Array<string> | null
}>()

const show = ref(true)

watch([() => props.error, () => props.errors], value => {
  show.value = value ? true : false
})
</script>

<template>
  <div>
    <slot></slot>

    <Transition name="fade">
      <div v-if="(error || errors?.length) && show" class="alert alert-danger">
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

.alert-danger ul {
  margin: 0;
  padding-left: 1em;
}

.alert-danger li {
  color: #721c24;
}

.close {
  float: right;
  text-decoration: none;
  cursor: pointer;
  color: #721c24;
  margin-left: 1em;
}
</style>
