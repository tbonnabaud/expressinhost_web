<script setup lang="ts">
import { onMounted, onUnmounted, ref } from 'vue'

const images = ref([...Array(6).keys()].map(i => `/banners/banner${i}.png`))
const currentIndex = ref(0)
const slideshowInterval = ref(0)

onMounted(() => {
  slideshowInterval.value = window.setInterval(nextImage, 7000)
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
    <Transition name="slide" mode="out-in">
      <img
        :key="currentIndex"
        :src="images[currentIndex]"
        :alt="'Slideshow Image ' + currentIndex"
      />
    </Transition>

    <a class="prev" @click="prevImage">&#10094;</a>
    <a class="next" @click="nextImage">&#10095;</a>
  </div>
</template>

<style scoped>
.slide-enter-active,
.slide-leave-active {
  transition: all 0.5s ease;
}

.slide-enter-from,
.slide-leave-to {
  opacity: 0;
}

.slider {
  position: relative;
  overflow: hidden;
  height: 320px;
  background-color: #000;
}

.slider img {
  width: 100%;
  height: 320px;
  object-fit: cover;
}

.prev,
.next {
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
  text-decoration: none;
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
