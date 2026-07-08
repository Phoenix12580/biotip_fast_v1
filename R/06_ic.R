getIc <-
function (counts, sampleL, genes, output = c("Ic", "PCCg", "PCCs"), 
    fun = c("cor", "BioTIP"), shrink = TRUE, use = c("everything", 
        "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs"), 
    PCC_sample.target = 1) 
{
    if (class(genes) != "character") 
        stop("genes have to be a character of gene symbols, i.e. \n                                      genes have to be a subset of row.names(counts)")
    output <- match.arg(output)
    fun <- match.arg(fun)
    use <- match.arg(use)
    if (class(PCC_sample.target) == "numeric") 
        if ((PCC_sample.target < 0) | (PCC_sample.target > 1)) {
            stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
        }
        else if (class(PCC_sample.target) == "character") 
            if (!PCC_sample.target %in% c("none", "zero", "average", 
                "half")) {
                stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
            }
    PCC_gene.target = "zero"
    subsetC = subset(counts, row.names(counts) %in% genes)
    subsetC = lapply(sampleL, function(x) subsetC[, as.character(x)])
    if (fun == "BioTIP" & PCC_gene.target == "none") 
        warning("You are not really calling BioTIP function without a proper setting of PCC_gene.target !")
    get_PCCg <- function() {
        if (fun == "BioTIP") {
            res = lapply(subsetC, function(x) avg.cor.shrink(x, 
                MARGIN = 1, shrink = shrink, abs = TRUE, target = PCC_gene.target))
            res = unlist(res)
        }
        else {
            res = lapply(subsetC, function(x) abs(cor(t(x), use = use)))
            for (i in seq_along(res)) res[[i]][upper.tri(res[[i]], 
                diag = FALSE)]
            res = sapply(res, function(x) mean(x, na.rm = TRUE))
        }
        names(res) = names(sampleL)
        res
    }
    get_PCCs <- function() {
        if (fun == "BioTIP") {
            res = lapply(subsetC, function(x) avg.cor.shrink(x, 
                MARGIN = 2, shrink = shrink, abs = FALSE, target = PCC_sample.target))
            res = unlist(res)
        }
        else {
            res = lapply(subsetC, function(x) cor(x, use = use))
            for (i in seq_along(res)) res[[i]][upper.tri(res[[i]], 
                diag = FALSE)]
            res = sapply(res, function(x) mean(x, na.rm = TRUE))
        }
        names(res) = names(sampleL)
        res
    }
    if (output == "Ic") {
        PCCg = get_PCCg()
        PCCs = get_PCCs()
        toplot = PCCg/PCCs
        names(toplot) = names(sampleL)
        return(toplot)
    }
    else if (output == "PCCg") {
        return(get_PCCg())
    }
    else if (output == "PCCs") {
        return(get_PCCs())
    }
}
simulation_Ic_sample <-
function (counts, sampleNo, Ic = NULL, genes, B = 1000, ylim = NULL, 
    main = "simulation of samples", fun = c("cor", "BioTIP"), 
    shrink = TRUE, use = c("everything", "all.obs", "complete.obs", 
        "na.or.complete", "pairwise.complete.obs"), output = c("Ic", 
        "PCCg", "PCCs"), plot = FALSE, PCC_sample.target = 1) 
{
    output <- match.arg(output)
    fun <- match.arg(fun)
    use <- match.arg(use)
    if (class(PCC_sample.target) == "numeric") 
        if ((PCC_sample.target < 0) | (PCC_sample.target > 1)) {
            stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero', 'average', 'half'")
        }
        else if (class(PCC_sample.target) == "character") 
            if (!PCC_sample.target %in% c("none", "zero", "average", 
                "half")) {
                stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero', 'average', 'half'")
            }
    PCC_gene.target = "zero"
    sampleL = lapply(1:B, function(x) sample(colnames(counts), 
        sampleNo))
    tmp = sapply(1:B, function(x) getIc(counts, sampleL = sampleL[x], 
        genes = genes, output = output, fun = fun, shrink = shrink, 
        use = use, PCC_sample.target = PCC_sample.target))
    if (plot) {
        p_v = length(tmp[tmp > Ic])/B
        den = density(tmp)
        xmin = min(Ic, den$x)
        xmax = max(Ic, den$x)
        plot(den, main = main, xlim = c(xmin, xmax), ylim = ylim)
        abline(v = Ic, col = "red", lty = 2)
        x = max(den$x) - 0.20000000000000001 * diff(range(den$x))
        if (p_v == 0) 
            p_v = paste("<", 1/B)
        text(x, max(den$y) - 0.050000000000000003, paste("P = ", 
            p_v))
    }
    return(tmp)
}
simulation_Ic <-
function (obs.x, sampleL, counts, B = 1000, fun = c("cor", "BioTIP"), 
    shrink = TRUE, use = c("everything", "all.obs", "complete.obs", 
        "na.or.complete", "pairwise.complete.obs"), output = c("Ic", 
        "PCCg", "PCCs"), PCC_sample.target = 1, n_cores = 3, 
    doParallel = TRUE) 
{
    fun <- match.arg(fun)
    use <- match.arg(use)
    if (n_cores > 0 & doParallel == TRUE) {
        n_cores_avail <- parallel::detectCores()
        if (n_cores > n_cores_avail - 1) 
            n_cores = n_cores_avail - 1
        cat("nCore: ", n_cores, "\n")
    }
    if (class(PCC_sample.target) == "numeric") 
        if ((PCC_sample.target < 0) | (PCC_sample.target > 1)) {
            stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
        }
        else if (class(PCC_sample.target) == "character") 
            if (!PCC_sample.target %in% c("none", "zero", "average", 
                "half")) {
                stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
            }
    PCC_gene.target = "zero"
    output <- match.arg(output)
    random = sapply(1:B, function(x) sample(row.names(counts), 
        obs.x))
    m <- matrix(nrow = length(sampleL), ncol = B)
    if (!doParallel) {
        pb <- txtProgressBar(min = 0, max = B, style = 3)
        for (i in 1:B) {
            setTxtProgressBar(pb, i)
            m[, i] <- BioTIP::getIc(counts, sampleL = sampleL, 
                genes = random[, i], output = output, fun = fun, 
                shrink = shrink, use = use, PCC_sample.target = PCC_sample.target)
            if (i == B) 
                cat("Done!\n")
        }
        close(pb)
    }
    if (doParallel) {
        cluster <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cluster)
        m <- foreach(i = 1:B, .combine = cbind) %dopar% {
            BioTIP::getIc(counts, sampleL = sampleL, genes = random[, 
                i], output = output, fun = fun, shrink = shrink, 
                use = use, PCC_sample.target = PCC_sample.target)
        }
        parallel::stopCluster(cl = cluster)
    }
    row.names(m) = names(sampleL)
    return(m)
}
getIc.new <-
function (X, method = c("BioTIP", "Ic"), PCC_sample.target = 1, 
    output = c("Ic", "PCCg", "PCCs")) 
{
    PCC_gene.target = "zero"
    method = match.arg(method)
    if (class(PCC_sample.target) == "numeric") 
        if ((PCC_sample.target < 0) | (PCC_sample.target > 1)) {
            stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
        }
        else if (class(PCC_sample.target) == "character") 
            if (!PCC_sample.target %in% c("none", "zero", "average", 
                "half")) {
                stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
            }
    output <- match.arg(output)
    shrink = (method == "BioTIP")
    if (output == "Ic") {
        numerator = avg.cor.shrink(X, MARGIN = 1, shrink = shrink, 
            abs = TRUE, target = PCC_gene.target)
        denominator = avg.cor.shrink(X, MARGIN = 2, shrink = shrink, 
            abs = FALSE, target = PCC_sample.target)
        return(numerator/denominator)
    }
    if (output == "PCCg") 
        return(avg.cor.shrink(X, MARGIN = 1, shrink = shrink, 
            abs = TRUE, target = PCC_gene.target))
    if (output == "PCCs") 
        return(avg.cor.shrink(X, MARGIN = 2, shrink = shrink, 
            abs = FALSE, target = PCC_sample.target))
}
