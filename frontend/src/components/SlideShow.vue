<script setup lang="ts">
import { onMounted, onUnmounted, ref } from 'vue'

const images = ref(['/banners/banner0.png', '/banners/banner1.png'])
const currentIndex = ref(0)
const slideshowInterval = ref(0)

onMounted(() => {
  slideshowInterval.value = setInterval(nextImage, 7000)
})

onUnmounted(() => clearInterval(slideshowInterval.value))

function nextImage() {
  currentIndex.value = (currentIndex.value + 1) % images.value.length
}

function prevImage() {
  currentIndex.value =
    (currentIndex.value - 1 + images.value.length) % images.value.length
}
</script>

<template>
  <div class="slider">
    <img :src="images[currentIndex]" alt="Slideshow Image" />

    <a class="selector prev" @click="prevImage">&#10094;</a>
    <a class="selector next" @click="nextImage">&#10095;</a>
  </div>
</template>

<style scoped>
.slider {
  position: relative;
  overflow: hidden;
}

.slider img {
  width: 100%;
  height: auto;
}

.selector,
.selector:link,
.selector:visited,
.selector:hover,
.selector:active {
  text-decoration: none;
}

.selector {
  cursor: pointer;
  position: absolute;
  top: 50%;
  width: auto;
  padding: 16px;
  margin-top: -50px;
  color: white;
  font-weight: bold;
  font-size: 20px;
  border-radius: 0 3px 3px 0;
  user-select: none;
  -webkit-user-select: none;
}

.prev {
  left: 0;
  border-radius: 3px 0 0 3px;
}

.next {
  right: 0;
  border-radius: 3px 0 0 3px;
}

.prev:hover,
.next:hover {
  background-color: rgba(0, 0, 0, 0.8);
}
</style>
