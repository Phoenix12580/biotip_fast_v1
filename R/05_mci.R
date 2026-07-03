getMaxMCImember <-
function (membersL, MCIl, minsize = 1, n = 1) 
{
    if (n < 1 | class(n) != "numeric") 
        stop("please provide a >= 1 numeric for n")
    if (is.null(names(membersL))) 
        names(membersL) <- 1:length(membersL)
    n <- round(n)
    listn = names(membersL)
    if (minsize >= 1) {
        minsize = minsize - 1
        CIl = lapply(seq_along(membersL), function(x) ifelse(table(membersL[[x]]) > 
            minsize, MCIl[[x]], NA))
        module_keep = lapply(seq_along(membersL), function(x) names(table(membersL[[x]])[table(membersL[[x]]) > 
            (minsize)]))
        membersL = lapply(seq_along(membersL), function(x) membersL[[x]][membersL[[x]] %in% 
            module_keep[[x]]])
    }
    else {
        stop("please provide a minimum number of the cluster of interest,  \n         which should be an integer that is larger than 0")
    }
    if (n >= 1) {
        idx = lapply(CIl, which.max)
        maxCI = lapply(seq_along(idx), function(x) names(membersL[[x]][membersL[[x]] == 
            idx[x]]))
        names(maxCI) = listn
        names(idx) = listn
        results <- list(idx = idx, members = maxCI)
    }
    names(CIl) = listn
    if (n > 1) {
        for (j in 2:n) {
            x <- unlist(lapply(idx, length))
            x <- names(x)[x > 0]
            if (length(x) > 0) {
                for (i in x) {
                  CIl[[i]][idx[[i]]] <- NA
                  idx[[i]] = c(idx[[i]], which.max(CIl[[i]]))
                }
            }
            results[[j + 1]] <- lapply(seq_along(idx), function(x) unlist(names(membersL[[x]][membersL[[x]] == 
                idx[[x]][j]])))
            names(results[[j + 1]]) = listn
            names(results)[j + 1] <- paste0(j, "topest.members")
        }
    }
    results[["idx"]] <- idx
    return(results)
}
getMaxStats <-
function (membersL, idx) 
{
    if (any(is.null(names(idx))) | any(!names(idx) %in% names(membersL))) 
        stop("please make sure \"idx\" has names and all of its names is included in names of \"membersL\"")
    member_max = lapply(names(idx), function(x) membersL[[x]][idx[[x]]])
    names(member_max) = names(idx)
    member_max = member_max[lengths(member_max) > 0]
    member_max = sapply(member_max, function(x) mean(x[[1]]))
    return(member_max)
}
getTopMCI <-
function (modulesL, MCI1, membersL, min, n = 1) 
{
    maxMCIms <- getMaxMCImember(modulesL, MCI1, min)
    topMCI = getMaxStats(membersL, maxMCIms[[1]])
    topMCI = topMCI[order(topMCI, decreasing = TRUE)]
    topMCI = topMCI[1:n]
    return(topMCI)
}
getCTS <-
function (maxMCI, maxMCIms) 
{
    if (is.null(names(maxMCI))) {
        stop("No names for maxMCI. Please provide names.")
    }
    if (is.null(names(maxMCIms))) {
        stop("No names for maxMCIms. Please provide names.")
    }
    if (!all(names(maxMCI) %in% names(maxMCIms))) {
        stop("Names of maxMCI has to be in maxMCIms.")
    }
    CTS_list = vector(mode = "list", length = length(maxMCI))
    for (i in 1:length(maxMCI)) {
        CTS_list[[i]] <- maxMCIms[[names(maxMCI)[i]]]
        message(paste0("Length: ", length(CTS_list[[i]])))
    }
    names(CTS_list) <- names(maxMCI)
    return(CTS_list)
}
getMCI <-
function (groups, countsL, adjust.size = FALSE, fun = c("cor", 
    "BioTIP"), df = NULL) 
{
    fun <- match.arg(fun)
    PCC_gene.target = "zero"
    if (all(is.na(groups))) {
        warning("no loci in any of the state in the list given,  \n            please rerun getCluster_methods with a larger cutoff \n            or provide a list of loci")
    }
    else {
        if (all(sapply(groups, class) == "communities")) {
            membersL = lapply(groups, membership)
        }
        else if (any(is.na(groups)) & any(sapply(groups, class) == 
            "communities")) {
            removed = groups[is.na(groups)]
            groups = groups[!is.na(groups)]
            membersL = lapply(groups, membership)
        }
        else {
            membersL = groups
        }
        CIl = PCCol = PCCl = sdl = list()
        names(membersL) = names(groups)
        if (is.null(names(groups))) 
            warning("No names provided for \"groups\"")
        if (is.null(names(countsL))) 
            warning("No names provided for \"countsL\"")
        loop = names(membersL)
        for (i in loop) {
            test = membersL[[i]]
            if (all(is.na(test))) {
                CI = sdL = PCC = PCCo = NA
            }
            else {
                test.counts = countsL[[i]]
                m = lapply(1:max(test), function(x) subset(test.counts, 
                  row.names(test.counts) %in% names(test[test == 
                    x])))
                comple = lapply(1:max(test), function(x) subset(test.counts, 
                  !row.names(test.counts) %in% names(test[test == 
                    x])))
                names(m) = names(comple) = 1:max(test)
                if (fun == "cor") {
                  PCCo = lapply(names(comple), function(x) abs(cor(t(comple[[x]]), 
                    t(m[[x]]))))
                  PCCo_avg = sapply(PCCo, function(x) mean(x, 
                    na.rm = TRUE))
                  PCC = lapply(m, function(x) abs(cor(t(x))))
                  PCC_avg = sapply(PCC, function(x) (sum(x, na.rm = TRUE) - 
                    nrow(x))/(nrow(x)^2 - nrow(x)))
                }
                if (fun == "BioTIP") {
                  if (is.null(df)) {
                    PCCo_avg = lapply(names(comple), function(x) avg.cor.shrink(comple[[x]], 
                      Y = m[[x]], abs = TRUE, MARGIN = 1, target = PCC_gene.target))
                    PCCo_avg = unlist(PCCo_avg)
                    PCC_avg = lapply(m, function(x) avg.cor.shrink(x, 
                      Y = NULL, abs = TRUE, MARGIN = 1, target = PCC_gene.target))
                    PCC_avg = unlist(PCC_avg)
                  }
                  else {
                    M <- cor.shrink(df, Y = NULL, MARGIN = 1, 
                      target = PCC_gene.target)
                    PCCo_avg <- array(dim = length(m))
                    names(PCCo_avg) <- names(m)
                    for (j in 1:length(m)) {
                      PCCo_avg[j] <- mean(abs(M[rownames(comple[[j]]), 
                        rownames(m[[j]])]))
                    }
                    PCC_avg <- array(dim = length(m))
                    names(PCC_avg) <- names(m)
                    for (j in 1:length(m)) {
                      tmp <- M[rownames(m[[j]]), rownames(m[[j]])]
                      U <- upper.tri(tmp, diag = FALSE)
                      PCC_avg[j] <- mean(abs(U))
                    }
                  }
                }
                sdL = lapply(m, function(x) apply(x, 1, sd))
                if (adjust.size) {
                  MCI = mapply(function(x, y, z, w) mean(x) * 
                    (y/z) * sqrt(nrow(w)), sdL, PCC_avg, PCCo_avg, 
                    m)
                }
                else {
                  MCI = mapply(function(x, y, z) mean(x) * (y/z), 
                    sdL, PCC_avg, PCCo_avg)
                }
            }
            CIl[[i]] = MCI
            sdl[[i]] = sdL
            PCCl[[i]] = PCC_avg
            PCCol[[i]] = PCCo_avg
        }
        names(CIl) = names(sdl) = names(PCCl) = names(PCCol) = names(membersL)
        return(list(members = membersL, MCI = CIl, sd = sdl, 
            PCC = PCCl, PCCo = PCCol))
    }
}
simulationMCI <-
function (len, samplesL, df, adjust.size = FALSE, B = 1000, fun = c("cor", 
    "BioTIP"), M = NULL, n_cores = 3, doParallel = FALSE, progress = TRUE) 
{
    fun <- match.arg(fun)
    PCC_gene.target = "zero"
    if (is.null(names(samplesL))) 
        stop("please provide names for list countsL")
    countsL = lapply(samplesL, function(x) df[, as.character(x)])
    if (is.null(names(countsL))) 
        names(countsL) = names(samplesL)
    if (fun == "BioTIP") {
        if (is.null(M)) 
            M <- cor.shrink(df, Y = NULL, MARGIN = 1, shrink = TRUE, 
                target = PCC_gene.target)
    }
    else M = NULL
    random_ids <- replicate(B, sample(1:nrow(countsL[[1]]), len), 
        simplify = FALSE)
    mci_inner <- getMCI_inner
    worker <- function(i) {
        mci_inner(len, countsL, adjust.size, fun = fun, PCC_gene.target = PCC_gene.target, 
            M = M, random_id = random_ids[[i]])
    }
    run_serial <- function() {
        if (progress) 
            pb <- txtProgressBar(min = 0, max = B, style = 3)
        m <- matrix(nrow = length(samplesL), ncol = B)
        for (i in 1:B) {
            if (progress) 
                setTxtProgressBar(pb, i)
            m[, i] <- worker(i)
            if (i == B && progress) 
                cat("Done!\n")
        }
        if (progress) 
            close(pb)
        m
    }
    if (doParallel && B > 1 && n_cores > 1) {
        n_cores <- min(n_cores, max(1, parallel::detectCores() - 
            1), B)
        cluster <- tryCatch(parallel::makeCluster(n_cores), error = function(e) e)
        if (inherits(cluster, "error")) {
            warning("parallel cluster could not be created; falling back to serial: ", 
                conditionMessage(cluster), call. = FALSE)
            m <- run_serial()
        }
        else {
            on.exit(parallel::stopCluster(cluster), add = TRUE)
            cols <- parallel::parLapply(cluster, seq_len(B), worker)
            m <- do.call(cbind, cols)
        }
    }
    else {
        m <- run_serial()
    }
    row.names(m) = names(countsL)
    return(m)
}
getMCI_inner <-
function (members, countsL, adjust.size, fun = c("cor", "BioTIP"), 
    PCC_gene.target = "zero", M = NULL, random_id = NULL) 
{
    fun <- match.arg(fun)
    if (is.null(random_id)) 
        random_id = sample(1:nrow(countsL[[1]]), members)
    randomL = lapply(names(countsL), function(x) countsL[[x]][random_id, 
        ])
    comple = lapply(names(countsL), function(x) subset(countsL[[x]], 
        !row.names(countsL[[x]]) %in% row.names(randomL[[x]])))
    names(randomL) = names(comple) = names(countsL)
    if (fun == "BioTIP") {
        PCCo_avg = array(dim = length(countsL))
        names(PCCo_avg) = names(countsL)
        for (i in 1:length(PCCo_avg)) {
            PCCo_avg[i] <- mean(abs(M[rownames(comple[[i]]), 
                rownames(randomL[[i]])]))
        }
        PCC_avg = array(dim = length(countsL))
        names(PCC_avg) = names(countsL)
        for (i in 1:length(PCC_avg)) {
            PCC_avg[i] = mean(abs(M[rownames(randomL[[i]]), rownames(randomL[[i]])]))
        }
    }
    else if (fun == "cor") {
        PCCo = lapply(names(comple), function(x) abs(cor(t(comple[[x]]), 
            t(randomL[[x]]))))
        PCCo_avg = sapply(PCCo, function(x) mean(x, na.rm = TRUE))
        PCC = lapply(randomL, function(x) abs(cor(t(x))))
        PCC_avg = sapply(PCC, function(x) (sum(x, na.rm = TRUE) - 
            nrow(x))/(nrow(x)^2 - nrow(x)))
    }
    sdL = lapply(randomL, function(x) apply(x, 1, sd))
    if (adjust.size) {
        MCI = mapply(function(x, y, z, w) mean(x) * (y/z) * sqrt(members), 
            sdL, PCC_avg, PCCo_avg, members)
    }
    else {
        MCI = mapply(function(x, y, z) mean(x) * (y/z), sdL, 
            PCC_avg, PCCo_avg)
    }
    return(MCI)
}
getNextMaxStats <-
function (membersL, idL = NULL, whoisnext, which.next = 2) 
{
    if (is.null(idL)) 
        stop("please provide idL, usually the 'idx' element returned by getMaxMCImember")
    score <- array(dim = length(whoisnext))
    names(score) <- whoisnext
    idx <- lapply(whoisnext, function(x) idL[[x]][which.next])
    for (i in 1:length(whoisnext)) {
        score[[i]] <- membersL[[whoisnext[i]]][idx[[i]]]
    }
    return(score)
}
