<script setup lang="ts">
defineProps<{ show: boolean }>()

defineEmits(['close'])
</script>

<template>
  <div>
    <Transition name="fade">
      <div v-if="show" class="alert alert-danger">
        <a class="close" aria-label="Close" @click="$emit('close')">
          <span aria-hidden="true">&times;</span>
        </a>

        <span>
          <slot></slot>
        </span>
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
  margin-bottom: 0.75rem;
  border: 1px solid transparent;
  border-radius: 0.25rem;
}

.alert-danger {
  color: #721c24;
  background-color: #f8d7da;
  border-color: #721c24;
}

.close {
  float: right;
  text-decoration: none;
  cursor: pointer;
  color: #721c24;
  margin-left: 1em;
}
</style>
