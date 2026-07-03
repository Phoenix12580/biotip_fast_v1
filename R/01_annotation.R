getBiotypes <-
function (full_gr, gencode_gr, intron_gr = NULL, minoverlap = 1L) 
{
    if (all(is(full_gr) != "GRanges")) 
        stop("please give full_gr as a \"GRanges\" object")
    if (all(is(gencode_gr) != "GRanges")) 
        stop("pealse give gencode_gr as a \"GRanges\" object")
    if (all(is(intron_gr) != "GRanges" & !is.null(intron_gr))) 
        stop("please give intron_gr as a \"GRanges\" object")
    hits = findOverlaps(full_gr, gencode_gr, type = "within", 
        minoverlap = minoverlap)
    full = as.data.frame(full_gr)
    full$type.fullOverlap = "de novo"
    idx = as.data.frame(mcols(full_gr[queryHits(hits)]))
    if (nrow(idx) != 0) {
        idx$biotype = as.data.frame(mcols(gencode_gr[subjectHits(hits)]))[, 
            1]
        idx_collapse = aggregate(as.list(idx["biotype"]), idx["Row.names"], 
            FUN = function(X) paste(unique(X), collapse = ",  "))
        idx_full = match(idx_collapse$Row.names, full$Row.names)
        full[idx_full, ]$type.fullOverlap = idx_collapse$biotype
    }
    hits = findOverlaps(full_gr, gencode_gr, minoverlap = minoverlap)
    overlaps <- pintersect(full_gr[queryHits(hits)], gencode_gr[subjectHits(hits)])
    percentOverlap <- width(overlaps)/width(gencode_gr[subjectHits(hits)])
    idx = as.data.frame(mcols(full_gr[queryHits(hits)]))
    idx$biotype = as.data.frame(mcols(gencode_gr[subjectHits(hits)]))
    idx_collapse = aggregate(as.list(idx["biotype"]), idx["Row.names"], 
        FUN = function(X) paste(unique(X), collapse = ",  "))
    full$type.partialOverlap = "de novo"
    idx_partial = match(idx_collapse$Row.names, full$Row.names)
    full[idx_partial, ]$type.partialOverlap = idx_collapse$biotype
    idx$percentOverlap = percentOverlap
    idx_50 = subset(idx, percentOverlap >= 0.5)
    idx_50collapse = aggregate(as.list(idx_50["biotype"]), idx_50["Row.names"], 
        FUN = function(X) paste(unique(X), collapse = ",  "))
    full$type.50Overlap = "de novo"
    idx_50 = match(idx_50collapse$Row.names, full$Row.names)
    full[idx_50, ]$type.50Overlap = idx_50collapse$biotype
    if (!is.null(intron_gr)) {
        hits = findOverlaps(full_gr, intron_gr)
        idx = unique(as.data.frame(mcols(full_gr[queryHits(hits)])))
        full$hasIntron = "no"
        idx_intron = match(idx$Row.names, full$Row.names)
        if (length(idx_intron) != 0) 
            full[idx_intron, ]$hasIntron = "yes"
    }
    else (full$hasIntron = NA)
    full$type.toPlot = ifelse(full$hasIntron == "yes" & full$type.50Overlap == 
        "protein_coding", "protein_coding_intron", full$type.50Overlap)
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("protein_coding", 
        x) & grepl("antisense", x), "protein_coding_antisense", 
        x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("protein_coding, ", 
        x), "protein_coding_mixed", x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl(",  protein_coding", 
        x), "protein_coding_mixed", x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("lincRNA", 
        x), "lincRNA", x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("antisense, ", 
        x), "antisense", x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl(",  antisense", 
        x), "antisense", x))
    label = c("protein_coding", "protein_coding_mixed", "lincRNA", 
        "antisense", "pseudogene,  processed_pseudogene", "pseudogene,  unprocessed_pseudogene", 
        "de novo", "protein_coding_antisense", "protein_coding_intron", 
        "miRNA")
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(!x %in% 
        label, "other_noncoding", x))
    return(full)
}
getReadthrough <-
function (gr, cod_gr) 
{
    full_table = data.frame(gr)
    overlapcount = countOverlaps(gr, cod_gr)
    completeoverlap = unique(subjectHits(findOverlaps(cod_gr, 
        GRanges(full_table$ID), type = "within")))
    if (length(completeoverlap) == 0) {
        full_table$readthrough = ifelse(overlapcount > 2, 1, 
            0)
    }
    else {
        full_table$readthrough = ifelse(overlapcount > 2 & row.names(completeoverlap) %in% 
            completeoverlap, 1, 0)
    }
    gr = GRanges(subset(full_table, readthrough == 1))
    idx = subset(full_table, readthrough == 1)$ID
    overlaps = as.data.frame(findOverlaps(gr, cod_gr))
    splitoverlaps = split(overlaps, f = overlaps$queryHits)
    table(sapply(splitoverlaps, nrow) > 1)
    cod_grL = sapply(splitoverlaps, function(x) cod_gr[x$subjectHits])
    overlapL = sapply(cod_grL, function(x) findOverlaps(x))
    notoverlap = sapply(overlapL, function(x) identical(queryHits(x), 
        subjectHits(x)))
    tmp = rep(TRUE, nrow(full_table))
    tmp[full_table$readthrough == 1] = notoverlap
    full_table$readthrough = ifelse(full_table$readthrough == 
        1 & !tmp, 1, 0)
    return(full_table)
}
