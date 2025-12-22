<script setup lang="ts">
defineProps<{ open: boolean; title: string }>()
const emit = defineEmits(['close'])

function handleClose() {
  emit('close')
}
</script>

<template>
  <Transition>
    <dialog v-if="open" @click.self="handleClose" open>
      <article>
        <header>
          <button aria-label="Close" rel="prev" @click="handleClose"></button>
          <p>
            <strong class="modal-title">{{ title }}</strong>
          </p>
        </header>

        <slot></slot>
      </article>
    </dialog>
  </Transition>
</template>

<style scoped>
header .modal-title {
  font-size: 1.5em;
}

.v-enter-active,
.v-leave-active {
  transition: opacity 0.3s ease;
}

.v-enter-from,
.v-leave-to {
  opacity: 0;
}
</style>
