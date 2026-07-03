sd_selection <-
function (df, samplesL, cutoff = 0.01, method = c("other", "reference", 
    "previous", "itself", "longitudinal reference"), control_df = NULL, 
    control_samplesL = NULL) 
{
    method = match.arg(method)
    if (is.null(names(samplesL))) 
        stop("please provide name to samplesL")
    if (any(!do.call(c, lapply(samplesL, as.character)) %in% 
        colnames(df))) 
        stop("please check if all sample names provided in \"samplesL\" are in colnames of \"df\"")
    if (any(lengths(samplesL) < 2)) 
        stop("please make sure there are at least one sample in every state")
    tmp = names(samplesL)
    samplesL = lapply(samplesL, as.character)
    test2 = sapply(tmp, function(x) apply(df[, as.character(samplesL[[x]])], 
        1, sd, na.rm = TRUE))
    if (method == "reference") {
        ref = as.character(samplesL[[1]])
        sdref = apply(df[, ref], 1, sd, na.rm = TRUE)
        sds = sapply(tmp, function(x) test2[, x]/sdref)
        names(sds) = tmp
    }
    else if (method == "other") {
        othersample = lapply(1:length(samplesL), function(x) do.call(c, 
            samplesL[-x]))
        names(othersample) = tmp
        sdother = sapply(tmp, function(x) apply(df[, as.character(othersample[[x]])], 
            1, sd, na.rm = TRUE))
        sds = lapply(tmp, function(x) test2[, x]/sdother[, x])
        names(sds) = tmp
    }
    else if (method == "previous") {
        warning("Using method \"previous\",  make sure sampleL is in the right order")
        sds = lapply(2:ncol(test2), function(x) test2[, x]/test2[, 
            x - 1])
        tmp = tmp[-1]
        names(sds) = tmp
    }
    else if (method == "itself") {
        if (cutoff > 1) 
            stop("Using method \"itself\",  cutoff must be smaller or equal to 1")
        sds = lapply(tmp, function(x) test2[, x])
        names(sds) = tmp
    }
    else if (method == "longitudinal reference") {
        if (is.null(control_df) | is.null(control_samplesL)) 
            stop("Using method \"longitudinal reference\",  \n           make sure \"control_df\" and \"sampleL\" are assigned")
        if (nrow(df) != nrow(control_df) | !all(row.names(df) %in% 
            row.names(control_df))) 
            stop("please make sure the row numbers of \"control_df\" \n           is the same as \"df\" and all transcripts in \"df\" are also in \"control_df\".")
        control = sapply(tmp, function(x) apply(control_df[, 
            as.character(control_samplesL[[x]])], 1, sd, na.rm = TRUE))
        sds = lapply(tmp, function(x) test2[, x]/control[, x])
        names(sds) = tmp
    }
    else {
        stop("method need to be selected from 'reference', 'other', 'previous',  'itself',  \n         or 'longitudinal reference' ")
    }
    if (cutoff <= 1) {
        topdf = nrow(df) * cutoff
        sdtop = lapply(tmp, function(x) names(sds[[x]][order(sds[[x]], 
            decreasing = TRUE)[1:topdf]]))
    }
    else {
        sdtop = lapply(tmp, function(x) names(sds[[x]][sds[[x]] > 
            cutoff]))
    }
    names(sdtop) = tmp
    subdf = lapply(tmp, function(x) df[, as.character(samplesL[[x]])])
    names(subdf) = tmp
    subm = lapply(names(subdf), function(x) subset(subdf[[x]], 
        row.names(subdf[[x]]) %in% sdtop[[x]]))
    names(subm) = tmp
    for (i in seq_along(subm)) {
        if (any(is.na(subm[[i]]))) {
            a <- apply(subm[[i]], 1, function(x) sum(x, na.rm = TRUE))
            tmp <- which(is.na(a))
            if (length(tmp) > 0) 
                subm[[i]] <- subm[[i]][-tmp, ]
            b <- apply(subm[[i]], 2, function(x) sum(x, na.rm = TRUE))
            tmp <- which(is.na(b))
            if (length(tmp) > 0) 
                subm[[i]] <- subm[[i]][, -tmp]
        }
    }
    return(subm)
}
optimize.sd_selection <-
function (df, samplesL, B = 100, percent = 0.80000000000000004, 
    times = 0.80000000000000004, cutoff = 0.01, method = c("other", 
        "reference", "previous", "itself", "longitudinal reference"), 
    control_df = NULL, control_samplesL = NULL, n_cores = 3, 
    doParallel = TRUE) 
{
    if (n_cores > 0 & doParallel == TRUE) {
        n_cores_avail <- parallel::detectCores()
        if (n_cores > n_cores_avail - 1) 
            n_cores = n_cores_avail - 1
        cat("nCore: ", n_cores, "\n")
    }
    method = match.arg(method)
    if (is.null(names(samplesL))) 
        stop("please provide name to samplesL")
    if (any(!do.call(c, lapply(samplesL, as.character)) %in% 
        colnames(df))) 
        stop("please check if all sample names provided in \"samplesL\" are in colnames of \"df\"")
    if (any(lengths(samplesL) < 2)) 
        stop("please make sure there are at least one sample in every state")
    N.random = lapply(seq_along(samplesL), function(x) matrix(0, 
        nrow = nrow(df), ncol = B))
    for (i in seq_along(N.random)) {
        row.names(N.random[[i]]) = row.names(df)
    }
    N = lengths(samplesL)
    k = round(N * percent)
    clusterID = names(samplesL)
    names(N.random) = clusterID
    if (method == "other") {
        samplesL = lapply(samplesL, as.character)
        othersample = lapply(seq_along(clusterID), function(x) do.call(c, 
            samplesL[-x]))
        names(othersample) = clusterID
        sdother = sapply(clusterID, function(x) apply(df[, othersample[[x]]], 
            1, function(y) sd(y, na.rm = TRUE)))
    }
    if (!doParallel) {
        pb <- txtProgressBar(min = 0, max = 100, style = 3)
        for (i in c(1:B)) {
            setTxtProgressBar(pb, i)
            random_sample_id = lapply(seq_along(k), function(x) sample(1:N[[x]], 
                k[[x]]))
            names(random_sample_id) = names(samplesL)
            selected_matrix = lapply(names(samplesL), function(x) df[, 
                samplesL[[x]][random_sample_id[[x]]]])
            test2 = sapply(selected_matrix, function(x) apply(x, 
                1, sd, na.rm = TRUE))
            clusterID = names(samplesL)
            colnames(test2) = clusterID
            if (method == "reference") {
                sdref = test2[, 1]
                sds = sapply(clusterID, function(x) test2[, x]/sdref)
                names(sds) = clusterID
            }
            else if (method == "other") {
                sds = lapply(clusterID, function(x) test2[, x]/sdother[, 
                  x])
                names(sds) = clusterID
            }
            else if (method == "previous") {
                cat("Using method \"previous\",  make sure sampleL is in the right order")
                sds = lapply(2:ncol(test2), function(x) test2[, 
                  x]/test2[, x - 1])
                clusterID <- clusterID[-1]
                names(sds) = clusterID
            }
            else if (method == "itself") {
                if (cutoff > 1) 
                  stop("Using method \"itself\",  cutoff must be smaller or equal to 1")
                sds = lapply(clusterID, function(x) test2[, x])
                names(sds) = clusterID
            }
            else if (method == "longitudinal reference") {
                if (is.null(control_df) | is.null(control_samplesL)) 
                  stop("Using method \"longitudinal reference\",  \n\t\t\t\t\t make sure \"control_df\" and \"sampleL\" are assigned")
                if (nrow(df) != nrow(control_df) | !all(row.names(df) %in% 
                  row.names(control_df))) 
                  stop("please make sure the row numbers of \"control_df\" \n\t\t\t\t\t is the same as \"df\" and all transcripts in \"df\" are also in \"control_df\".")
                control = sapply(clusterID, function(x) apply(control_df[, 
                  as.character(control_samplesL[[x]])], 1, sd, 
                  na.rm = TRUE))
                sds = lapply(clusterID, function(x) test2[, x]/control[, 
                  x])
                names(sds) = clusterID
            }
            else {
                stop("method need to be selected from 'reference', \n\t\t\t\t   'other', 'previous', 'itself', 'longitudinal reference'")
            }
            if (cutoff <= 1) {
                topdf = nrow(selected_matrix[[1]]) * cutoff
                sdtop = lapply(clusterID, function(x) names(sds[[x]][order(sds[[x]], 
                  decreasing = TRUE)[1:topdf]]))
            }
            else {
                sdtop = lapply(clusterID, function(x) names(sds[[x]][sds[[x]] > 
                  cutoff]))
            }
            names(sdtop) = clusterID
            names(N.random) = clusterID
            for (j in clusterID) {
                N.random[[j]][sdtop[[j]], i] = 1
            }
        }
        close(pb)
    }
    if (doParallel) {
        cluster <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cluster)
        n_iterations <- B
        sdtop_results <- list()
        sdtop_results <- foreach(i = 1:n_iterations) %dopar% 
            {
                random_sample_id = lapply(seq_along(k), function(x) sample(1:N[[x]], 
                  k[[x]]))
                names(random_sample_id) = names(samplesL)
                selected_matrix = lapply(names(samplesL), function(x) df[, 
                  samplesL[[x]][random_sample_id[[x]]]])
                test2 = sapply(selected_matrix, function(x) apply(x, 
                  1, sd, na.rm = TRUE))
                colnames(test2) = clusterID
                if (method == "reference") {
                  sdref = test2[, 1]
                  sds = sapply(clusterID, function(x) test2[, 
                    x]/sdref)
                  names(sds) = clusterID
                }
                else if (method == "other") {
                  sds = lapply(clusterID, function(x) test2[, 
                    x]/sdother[, x])
                  names(sds) = clusterID
                }
                else if (method == "previous") {
                  cat("Using method \"previous\",  make sure sampleL is in the right order")
                  sds = lapply(2:ncol(test2), function(x) test2[, 
                    x]/test2[, x - 1])
                  clusterID <- clusterID[-1]
                  names(sds) = clusterID
                }
                else if (method == "itself") {
                  if (cutoff > 1) 
                    stop("Using method \"itself\",  cutoff must be smaller or equal to 1")
                  sds = lapply(clusterID, function(x) test2[, 
                    x])
                  names(sds) = clusterID
                }
                else if (method == "longitudinal reference") {
                  if (is.null(control_df) | is.null(control_samplesL)) 
                    stop("Using method \"longitudinal reference\",  \n             make sure \"control_df\" and \"sampleL\" are assigned")
                  if (nrow(df) != nrow(control_df) | !all(row.names(df) %in% 
                    row.names(control_df))) 
                    stop("please make sure the row numbers of \"control_df\" \n             is the same as \"df\" and all transcripts in \"df\" are also in \"control_df\".")
                  control = sapply(clusterID, function(x) apply(control_df[, 
                    as.character(control_samplesL[[x]])], 1, 
                    sd, na.rm = TRUE))
                  sds = lapply(clusterID, function(x) test2[, 
                    x]/control[, x])
                  names(sds) = clusterID
                }
                else {
                  stop("method need to be selected from 'reference', \n           'other', 'previous', 'itself', 'longitudinal reference'")
                }
                if (cutoff <= 1) {
                  topdf = nrow(df) * cutoff
                  sdtop = lapply(clusterID, function(x) names(sds[[x]][order(sds[[x]], 
                    decreasing = TRUE)[1:topdf]]))
                }
                else {
                  sdtop = lapply(clusterID, function(x) names(sds[[x]][sds[[x]] > 
                    cutoff]))
                }
                names(sdtop) = clusterID
                sdtop_results[[i]] <- sdtop
            }
        parallel::stopCluster(cl = cluster)
        for (i in c(1:B)) {
            sdtop = sdtop_results[[i]]
            for (j in clusterID) {
                N.random[[j]][sdtop[[j]], i] = 1
            }
        }
    }
    times = times * B
    stable = lapply(N.random, function(x) row.names(x[rowSums(x) > 
        times, ]))
    names(stable) = clusterID
    subdf = lapply(clusterID, function(x) df[, as.character(samplesL[[x]])])
    names(subdf) = clusterID
    subm = lapply(names(subdf), function(x) subset(subdf[[x]], 
        row.names(subdf[[x]]) %in% stable[[x]]))
    names(subm) = clusterID
    return(subm)
}
