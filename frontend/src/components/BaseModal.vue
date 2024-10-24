<script setup lang="ts">
defineProps<{ open: boolean; title: string }>()
const emit = defineEmits(['close'])

function close() {
  emit('close')
}
</script>

<template>
  <dialog :open="open" @click.self="close">
    <article>
      <header>
        <button aria-label="Close" rel="prev" @click="close"></button>
        <p>
          <strong class="modal-title">{{ title }}</strong>
        </p>
      </header>

      <slot></slot>
    </article>
  </dialog>
</template>

<style scoped>
header .modal-title {
  font-size: 1.5em;
}
</style>
