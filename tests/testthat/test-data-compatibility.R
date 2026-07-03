test_that("bundled GRanges datasets do not require the obsolete Seqinfo package", {
  for (dataset in c("gencode", "intron", "ILEF", "cod")) {
    data(list = dataset, package = "BioTIP")
    obj <- get(dataset)
    seqinfo_class <- attr(attributes(obj)$seqinfo, "class")

    expect_false(identical(attr(seqinfo_class, "package"), "Seqinfo"), info = dataset)
    expect_true(methods::is(GenomicRanges::GRanges(obj), "GRanges"), info = dataset)
  }
})
