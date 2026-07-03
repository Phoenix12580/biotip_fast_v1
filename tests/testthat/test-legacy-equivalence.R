test_that("simulationMCI exposes conservative acceleration controls", {
  expect_true("n_cores" %in% names(formals(simulationMCI)))
  expect_true("doParallel" %in% names(formals(simulationMCI)))
  expect_true("progress" %in% names(formals(simulationMCI)))
})

test_that("simulationMCI preserves seeded BioTIP output", {
  set.seed(12)
  df <- matrix(rnorm(24 * 15), nrow = 24)
  rownames(df) <- paste0("g", seq_len(nrow(df)))
  colnames(df) <- paste0("s", seq_len(ncol(df)))
  samplesL <- split(colnames(df), rep(paste0("state", 1:3), each = 5))
  M <- cor.shrink(df, Y = NULL, MARGIN = 1, shrink = TRUE, target = "zero")

  expected <- structure(c(
    3.73219917440204, 5.33879384891026, 3.95970636614022,
    3.6348302655973, 3.9712601793812, 4.61138064778495,
    4.44893995886252, 3.38712042929145, 3.52587108539964,
    2.83577209098599, 4.5361530713525, 3.89790472721299,
    3.60916537617761, 4.46293460787872, 6.17409571166059,
    3.88893441878673, 4.14506094711676, 3.57569278776616,
    3.80985585462572, 4.81240219440302, 3.58879607315715
  ), dim = c(3L, 7L), dimnames = list(c("state1", "state2", "state3"), NULL))

  set.seed(102)
  default <- simulationMCI(4, samplesL, df, B = 7, fun = "BioTIP", M = M, progress = FALSE)
  set.seed(102)
  parallel <- simulationMCI(4, samplesL, df, B = 7, fun = "BioTIP", M = M,
                            doParallel = TRUE, n_cores = 2, progress = FALSE)

  expect_equal(default, expected)
  expect_equal(parallel, expected)
})

test_that("getNetwork computes one correlation test per state", {
  set.seed(11)
  network_input <- list(
    state_a = matrix(rnorm(10 * 5), nrow = 10),
    state_b = matrix(rnorm(9 * 5), nrow = 9)
  )
  for (nm in names(network_input)) {
    rownames(network_input[[nm]]) <- paste0(nm, "_g", seq_len(nrow(network_input[[nm]])))
    colnames(network_input[[nm]]) <- paste0(nm, "_s", seq_len(ncol(network_input[[nm]])))
  }

  assign(".__biotip_corr_test_call_count__", 0L, envir = .GlobalEnv)
  trace("corr.test",
        where = asNamespace("BioTIP"),
        tracer = quote(assign(
          ".__biotip_corr_test_call_count__",
          get(".__biotip_corr_test_call_count__", envir = .GlobalEnv) + 1L,
          envir = .GlobalEnv
        )),
        print = FALSE)
  on.exit(untrace("corr.test", where = asNamespace("BioTIP")), add = TRUE)
  on.exit(rm(".__biotip_corr_test_call_count__", envir = .GlobalEnv), add = TRUE)

  suppressMessages(getNetwork(network_input, fdr = 1))

  call_count <- get(".__biotip_corr_test_call_count__", envir = .GlobalEnv)
  expect_equal(call_count, length(network_input))
})
