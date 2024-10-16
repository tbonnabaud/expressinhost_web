<script setup lang="ts"></script>

<template>
  <form>
    <section>
      <h2>Data files</h2>

      <div class="flex-container">
        <div class="input-file">
          <label for="fasta">Sequence file (FASTA)</label>
          <input type="file" name="" id="fasta" />
        </div>

        <div class="input-file">
          <label for="clustal">Alignment file (CLUSTAL, optional)</label>
          <input type="file" name="" id="clustal" />
        </div>
      </div>
    </section>

    <section>
      <h2>Sequences</h2>

      <table>
        <tbody>
          <tr>
            <td>Sequence1</td>
            <td class="select-cell">
              <select name="" id="">
                <option value="">Table1</option>
                <option value="">Table2</option>
                <option value="">Table3</option>
              </select>
            </td>
          </tr>
          <tr>
            <td>Sequence2</td>
            <td class="select-cell">
              <select name="" id="">
                <option value="">Table1</option>
                <option value="">Table2</option>
                <option value="">Table3</option>
              </select>
            </td>
          </tr>
          <tr>
            <td>Sequence3</td>
            <td class="select-cell">
              <select name="" id="">
                <option value="">Table1</option>
                <option value="">Table2</option>
                <option value="">Table3</option>
              </select>
            </td>
          </tr>
        </tbody>
      </table>
    </section>

    <section>
      <h2>Host organism</h2>

      <div class="flex-container">
        <div id="organism">
          <label for="">Host organism</label>
          <select name="">
            <option value="">Escherichia_coli</option>
            <option value="">Danio_rerio</option>
          </select>
        </div>

        <div id="tableName">
          <label for="">Table name</label>
          <select name="">
            <option value="">Escherichia_coli_Default</option>
            <option value="">Escherichia_coli_Custom</option>
          </select>
        </div>
      </div>
    </section>

    <section>
      <h2>Mode</h2>

      <div id="mode-selector">
        <input type="radio" name="mode" value="direct_mapping" checked />
        <label>Direct mapping </label>
        <input
          type="radio"
          name="mode"
          value="optimisation_and_conservation_1"
        />
        <label>Optimisation and conservation 1 </label>
        <input
          type="radio"
          name="mode"
          value="optimisation_and_conservation_2"
        />
        <label>Optimisation and conservation 1 </label>
      </div>
    </section>

    <section>
      <h2>Thresholds</h2>

      <div class="flex-container">
        <div class="input-range">
          <label for="">Slow speed threshold</label>
          <input type="range" name="" id="" min="0" max="1" step="0.01" />
        </div>

        <div class="input-range">
          <label for="">Conservation threshold</label>
          <input type="range" name="" id="" min="0" max="1" step="0.01" />
        </div>
      </div>
    </section>

    <hr />

    <button id="run">Run tuning</button>
  </form>
</template>

<style scoped>
.input-file {
  width: 100%;
}

.input-range input {
  width: 90%;
  margin-bottom: 1em;
}

#organism,
#tableName {
  flex: auto;
  margin: 10px;
  width: 100%;
}

.input-range {
  flex: auto;
  text-align: center;
}

td select {
  width: 100%;
}

td {
  width: 50%;
}

#run {
  width: 100%;
  /* background-color: #294f29; */
}

#mode-selector {
  text-align: left;
}
</style>
