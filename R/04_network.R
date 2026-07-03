getCluster <-
function (igraphL, steps = 4) 
{
    if (length(steps) == 1 & steps%%1 == 0) {
        steps = rep(steps, length(igraphL))
    }
    else if (length(steps) != 1 | length(steps) != length(igraphL)) {
        stop("check step: must be postive integer(s) of length 1 or length of igraphL")
    }
    groups = list()
    for (i in seq_along(igraphL)) {
        if (nrow(as_data_frame(igraphL[[i]])) != 0) {
            groups[[i]] = cluster_walktrap(igraphL[[i]], weight = abs(E(igraphL[[i]])$weight), 
                steps = steps[i])
        }
        else {
            groups[[i]] = NA
        }
    }
    names(groups) = names(igraphL)
    return(groups)
}
getCluster_methods <-
function (igraphL, method = c("rw", "hcm", "km", "pam", "natural"), 
    cutoff = NULL) 
{
    method <- match.arg(method)
    if (method == "rw") {
        if (all(sapply(igraphL, class) != "igraph")) 
            stop("random walk clustering needs a list of igraph object \n           which can be obtained using getNetwork")
        if (!is.null(cutoff)) 
            if (cutoff%%1 != 0) 
                warning("Please provide a integer as \"cutoff\" for the cluster method random walk")
        if (is.null(cutoff)) 
            cutoff = 4
        groups = getCluster(igraphL, cutoff)
    }
    else if (method == "hcm") {
        if (all(!sapply(igraphL, class) %in% c("matrix", "data.frame"))) 
            stop("hierarchical clustering needs a list of matrix or data.frame as the 1st argument")
        if (is.null(cutoff)) 
            stop("hierarchical clustering needs \"cutoff\" \n                             to be assigned as the number of clusters wanted")
        testL = lapply(igraphL, function(x) corr.test(t(x), adjust = "fdr", 
            ci = FALSE)$r)
        groupsL = lapply(seq_along(testL), function(x) hclust(dist(testL[[x]]), 
            method = "complete"))
        par(mfrow = c(1, length(groupsL)))
        sapply(groupsL, function(x) plot(x))
        groups = lapply(groupsL, function(x) cutree(x, cutoff))
    }
    else if (method %in% c("km", "pam")) {
        if (all(!sapply(igraphL, class) %in% c("matrix", "data.frame"))) 
            stop("k-mediods or PAM clustering needs a list of matrix or data.frame as the 1st argument")
        if (is.null(cutoff)) 
            stop("hierarchical clustering needs \"cutoff\" \n                             to be assigned as the number of clusters wanted")
        testL = lapply(igraphL, function(x) corr.test(t(x), adjust = "fdr", 
            ci = FALSE)$r)
        groups = lapply(seq_along(testL), function(x) pam(testL[[x]], 
            cutoff, metric = "euclidean")$clustering)
    }
    else if (method == "natrual") {
        warning("selecting \"natural\" which will not use \"cutoff\" parameter")
        if (all(sapply(igraphL, class) != "igraph")) 
            stop("selecting \"natural\" which needs a list of igraph object \n           as the 1st argument which can be obtained using getNetwork")
        groups = lapply(seq_along(igraphL), function(x) components(igraphL[[x]])$membership)
    }
    else (stop("please select from \"rw\",  \"hcm\", \"km\",  \"pam\",  \"natrual\" as method"))
    return(groups)
}
getNetwork <-
function (optimal, fdr = 0.050000000000000003) 
{
    corrL = lapply(optimal, function(x) corr.test(t(x), adjust = "fdr", 
        ci = FALSE))
    rL = lapply(corrL, function(x) x$r)
    names(rL) = names(optimal)
    pL = lapply(corrL, function(x) x$p)
    if (is.null(names(rL))) 
        stop("give names to the input list")
    igraphL = list()
    for (i in names(rL)) {
        test = rL[[i]]
        test.p = pL[[i]]
        row.names(test) = gsub(".", "+", row.names(test), perl = FALSE, 
            fixed = TRUE)
        row.names(test.p) = gsub(".", "+", row.names(test.p), 
            perl = FALSE, fixed = TRUE)
        test[lower.tri(test, diag = TRUE)] = NA
        tmp = lapply(1:nrow(test), function(x) test[x, test.p[x, 
            ] < fdr])
        tmp_name = lapply(1:nrow(test), function(x) which(test.p[x, 
            ] < fdr))
        idx = which(lengths(tmp_name) == 1)
        for (j in idx) {
            names(tmp[[j]]) = names(tmp_name[[j]])
        }
        names(tmp) = row.names(test)
        edges = stack(do.call(c, tmp))
        edges = subset(edges, !is.na(edges$values))
        edges$node1 = str_split_fixed(edges$ind, "\\.", 2)[, 
            1]
        edges$node2 = str_split_fixed(edges$ind, "\\.", 2)[, 
            2]
        edges$node1 = gsub("+", ".", edges$node1, fixed = TRUE)
        edges$node2 = gsub("+", ".", edges$node2, fixed = TRUE)
        edges = edges[, c("node1", "node2", "values")]
        edges$weight = abs(edges$values)
        nodes = data.frame(unique(c(edges$node1, edges$node2)))
        message(paste0(i, ":", nrow(nodes), " nodes"))
        routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, 
            directed = FALSE)
        igraphL[[i]] = routes_igraph
    }
    return(igraphL)
}
