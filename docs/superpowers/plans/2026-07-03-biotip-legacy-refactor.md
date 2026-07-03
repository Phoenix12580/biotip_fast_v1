# BioTIP Legacy-Compatible Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor BioTIP internals for maintainability and speed while preserving legacy statistical outputs by default.

**Architecture:** Keep the package name and exported API compatible with BioTIP 1.26.0. Add characterization tests against the installed legacy package, split the monolithic implementation into focused R files, and add conservative speedups that preserve random draws and numeric results under the same seed.

**Tech Stack:** R 4.3.1, BioTIP 1.26.0 source, testthat 3, base R parallel, igraph, psych, Matrix.

---

### Task 1: Characterization Test Harness

**Files:**
- Create: `tests/testthat.R`
- Create: `tests/testthat/helper-legacy.R`
- Create: `tests/testthat/test-legacy-equivalence.R`
- Modify: `DESCRIPTION`

- [ ] **Step 1: Write failing test infrastructure**

Create `tests/testthat.R`:

```r
library(testthat)
library(BioTIP)

test_check("BioTIP")
```

Create `tests/testthat/helper-legacy.R`:

```r
legacy_ns <- new.env(parent = emptyenv())

load_legacy_functions <- function() {
  source(system.file("R", "BioTIP_update_06232025.R", package = "BioTIP.legacy"), local = legacy_ns)
}
```

Create `tests/testthat/test-legacy-equivalence.R` with a test that requires `simulationMCI()` to expose the new conservative acceleration arguments:

```r
test_that("simulationMCI exposes conservative acceleration controls", {
  expect_true("n_cores" %in% names(formals(simulationMCI)))
  expect_true("doParallel" %in% names(formals(simulationMCI)))
  expect_true("progress" %in% names(formals(simulationMCI)))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `D:/R/R-4.3.1/bin/Rscript.exe -e "devtools::test('E:/03_fig_copy/BioTIP_legacy_refactor')"`

Expected: FAIL because `simulationMCI()` does not have `n_cores`, `doParallel`, or `progress`.

- [ ] **Step 3: Add testthat to package metadata**

Modify `DESCRIPTION` so `Suggests` includes `testthat (>= 3.0.0)` and add `Config/testthat/edition: 3`.

- [ ] **Step 4: Commit**

Run:

```bash
git add DESCRIPTION tests
git commit -m "test: add legacy refactor characterization harness"
```

### Task 2: Conservative SimulationMCI Acceleration

**Files:**
- Modify: `R/BioTIP_update_06232025.R`
- Modify: `tests/testthat/test-legacy-equivalence.R`

- [ ] **Step 1: Expand the failing test**

Add this test:

```r
test_that("simulationMCI_fast path preserves seeded BioTIP output", {
  set.seed(12)
  df <- matrix(rnorm(24 * 15), nrow = 24)
  rownames(df) <- paste0("g", seq_len(nrow(df)))
  colnames(df) <- paste0("s", seq_len(ncol(df)))
  samplesL <- split(colnames(df), rep(paste0("state", 1:3), each = 5))
  M <- cor.shrink(df, Y = NULL, MARGIN = 1, shrink = TRUE, target = "zero")

  set.seed(102)
  expected <- simulationMCI(4, samplesL, df, B = 7, fun = "BioTIP", M = M,
                            doParallel = FALSE, progress = TRUE)
  set.seed(102)
  actual <- simulationMCI(4, samplesL, df, B = 7, fun = "BioTIP", M = M,
                          doParallel = TRUE, n_cores = 2, progress = FALSE)

  expect_equal(actual, expected)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `D:/R/R-4.3.1/bin/Rscript.exe -e "devtools::test('E:/03_fig_copy/BioTIP_legacy_refactor')"`

Expected: FAIL because `simulationMCI()` does not yet accept the new arguments.

- [ ] **Step 3: Implement minimal acceleration**

Change `simulationMCI()` to add `n_cores = 3, doParallel = FALSE, progress = TRUE`. Pre-generate `random_ids <- replicate(B, sample(...), simplify = FALSE)`, call a helper that accepts one random id, and use `parallel::parLapply()` only when `doParallel` is TRUE. Remove the hard `Sys.sleep(0.01)`.

- [ ] **Step 4: Run test to verify it passes**

Run: `D:/R/R-4.3.1/bin/Rscript.exe -e "devtools::test('E:/03_fig_copy/BioTIP_legacy_refactor')"`

Expected: PASS.

- [ ] **Step 5: Commit**

Run:

```bash
git add R/BioTIP_update_06232025.R tests/testthat/test-legacy-equivalence.R
git commit -m "refactor: add deterministic parallel simulationMCI path"
```

### Task 3: Network Calculation De-Duplication

**Files:**
- Modify: `R/BioTIP_update_06232025.R`
- Modify: `tests/testthat/test-legacy-equivalence.R`

- [ ] **Step 1: Add equivalence test**

Add a test that builds two toy state matrices, calls `getNetwork(fdr = 1)`, and verifies edge data frames are stable after the refactor.

- [ ] **Step 2: Run test**

Run: `D:/R/R-4.3.1/bin/Rscript.exe -e "devtools::test('E:/03_fig_copy/BioTIP_legacy_refactor')"`

Expected: PASS before implementation because it is a characterization test.

- [ ] **Step 3: Implement de-duplication**

Change `getNetwork()` to call `psych::corr.test(t(x), adjust = "fdr", ci = FALSE)` once per state, then extract `$r` and `$p`.

- [ ] **Step 4: Run test**

Run: `D:/R/R-4.3.1/bin/Rscript.exe -e "devtools::test('E:/03_fig_copy/BioTIP_legacy_refactor')"`

Expected: PASS.

- [ ] **Step 5: Commit**

Run:

```bash
git add R/BioTIP_update_06232025.R tests/testthat/test-legacy-equivalence.R
git commit -m "refactor: avoid duplicate network correlation tests"
```

### Task 4: Split Monolithic R File

**Files:**
- Create: `R/annotation.R`
- Create: `R/selection.R`
- Create: `R/network.R`
- Create: `R/mci.R`
- Create: `R/ic.R`
- Create: `R/correlation.R`
- Create: `R/plots.R`
- Create: `R/wrap.R`
- Modify: `R/BioTIP_update_06232025.R`

- [ ] **Step 1: Move functions by responsibility**

Move functions without changing bodies except already tested conservative speedups:

```text
annotation.R: getBiotypes, getReadthrough
selection.R: sd_selection, optimize.sd_selection
network.R: getNetwork, getCluster, getCluster_methods
mci.R: getMCI, getMCI_inner, simulationMCI, getMaxMCImember, getMaxStats, getTopMCI, getCTS, getNextMaxStats
ic.R: getIc, getIc.new, simulation_Ic, simulation_Ic_sample
correlation.R: avg.cor.shrink, cor.shrink
plots.R: plotBar_MCI, plotMaxMCI, plotIc, plot_Ic_Simulation, plot_MCI_Simulation, plot_SS_Simulation, plotIcSignificance
wrap.R: BioTIP.wrap
```

- [ ] **Step 2: Remove moved code from monolith**

Leave `R/BioTIP_update_06232025.R` as a short compatibility file with a comment explaining that the implementation has moved to focused modules.

- [ ] **Step 3: Run full tests**

Run: `D:/R/R-4.3.1/bin/Rscript.exe -e "devtools::test('E:/03_fig_copy/BioTIP_legacy_refactor')"`

Expected: PASS.

- [ ] **Step 4: Run package check**

Run: `D:/R/R-4.3.1/bin/Rscript.exe -e "devtools::check('E:/03_fig_copy/BioTIP_legacy_refactor', args='--no-manual')"`

Expected: No ERROR. Warnings/notes from old examples or docs are recorded but not mixed into behavior changes.

- [ ] **Step 5: Commit**

Run:

```bash
git add R tests DESCRIPTION
git commit -m "refactor: split BioTIP implementation into focused modules"
```
