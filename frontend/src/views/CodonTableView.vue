<script setup lang="ts">
import CodonTableForm from '@/components/codon-tables/CodonTableForm.vue'
import RequiredAuth from '@/components/RequiredAuth.vue'
import { store } from '@/lib/store'

const user = store.currentUser
</script>

<template>
  <main class="container-fluid">
    <template v-if="user">
      <h1>tRNA GCN tables</h1>
      <hr />

      <ul id="codonTableFormExplanations">
        <li>
          <strong>To create a table from scratch:</strong> select the table
          "(None)" &rarr; input organism name &rarr; input table name &rarr;
          modify the values in the amino acid boxes &rarr; click "Save as new".
        </li>
        <li>
          <strong>To create a table from an existing one:</strong> select the
          existing table &rarr; change the organism name and/or the table name
          &rarr; modify the values in the amino acid boxes &rarr; click "Save as
          new".
        </li>
        <li>
          <strong>To edit a custom table:</strong> select the custom table
          &rarr; modify the values in the amino acid boxes &rarr; click
          "Update". Each logged user can modify its own custom tables only.
          Default tables cannot be modified.
        </li>
        <li><strong>To reset all fields:</strong> click "Reset".</li>
        <li>
          <strong>To delete a custom table:</strong> select it and click
          "Delete".
        </li>
      </ul>

      <CodonTableForm />
    </template>

    <RequiredAuth v-else />
  </main>
</template>

<style scoped>
#codonTableFormExplanations li {
  font-style: italic;
  list-style-type: disc;
}
</style>
