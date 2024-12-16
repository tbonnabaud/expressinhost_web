<template>
  <div class="tooltip">
    <slot></slot>
    <span class="tooltiptext">
      <slot name="tooltip">Lorem ipsum</slot>
    </span>
  </div>
</template>

<style scoped>
.tooltip {
  position: relative;
  display: inline-block;
  /* cursor: help; */
}

.tooltip .tooltiptext {
  visibility: hidden;
  width: max-content;
  max-width: 300px;
  background-color: #56637f;
  border: 2px solid #40454e;
  color: #fff;
  text-align: left;
  border-radius: 6px;
  padding: 10px;
  opacity: 0;
  transition: opacity 0.5s;

  position: absolute;
  z-index: 1;
  top: -150%;
  left: 105%;
}

.tooltip:hover .tooltiptext {
  visibility: visible;
  opacity: 0.9;
}
</style>
