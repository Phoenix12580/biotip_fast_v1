plotBar_MCI <-
function (MCIl, ylim = NULL, nr = 1, nc = NULL, order = NULL, 
    minsize = 3, states = NULL, title.size = 30) 
{
    membersL = MCIl[[1]]
    MCI = MCIl[[2]]
    if (is.null(order)) {
        loop = names(membersL)
    }
    else {
        if (any(!order %in% names(membersL))) 
            stop("make sure all names provided in \"order\" is included in names of \"countsL\"")
        loop = order
    }
    if (!is.null(states)) 
        loop = states
    if (is.null(nc)) 
        nc = length(loop)
    par(mfrow = c(nr, nc))
    for (i in loop) {
        if (!i %in% names(MCI)) {
            mci = m = 0
        }
        else {
            mci = MCI[[i]]
            m = membersL[[i]]
            tmp = names(mci[is.na(mci)])
            if (length(tmp) != 0) 
                m = m[!m %in% tmp]
            nmembers = sapply(names(table(m)), function(x) length(m[m == 
                x]))
            mci = mci[!is.na(mci)]
            if (!minsize < 0 & minsize != 1) {
                mci = mci[!nmembers < minsize]
                nmembers = nmembers[!nmembers < minsize]
            }
            else {
                warning("\"minisize\" need to be a non")
            }
            if (length(mci) == 0) 
                mci = 0
        }
        mci[is.na(mci)] = 0
        bar = barplot(mci, col = rainbow(length(mci), alpha = 0.29999999999999999), 
            main = "", ylab = "DNB score", xlab = "modules", 
            ylim = ylim, cex.axis = 1.5, cex.names = 1.5, cex.main = 1.5, 
            cex.lab = 1.5)
        title(main = paste0(i, "\n", max(m), " modules"), ps = title.size)
        if (all(mci != 0)) 
            text(bar, mci, nmembers, cex = 1.5)
    }
}
plotMaxMCI <-
function (maxMCIms, MCIl, las = 0, order = NULL, states = NULL) 
{
    if (any(is.null(names(maxMCIms[[1]]))) | any(is.null(names(maxMCIms[[2]])))) 
        stop("Please give names for the 1st and 2nd element of the \"maxMCIms\" as well as \"MCIl\"")
    if (is.null(order)) {
        CI = sapply(names(maxMCIms[[1]]), function(x) MCIl[[x]][maxMCIms[[1]][[x]]])
        ln = names(maxMCIms[[1]])
        names(CI) = ln
    }
    else {
        if (any(!order %in% names(maxMCIms[[2]]))) 
            stop("make sure all names in \"order\" are in names of the 2nd element of \"maxMCIms\"")
        if (any(!names(maxMCIms[[2]]) %in% order)) 
            warning("not every state in \"simulation\" is plotted,  make sure \"order\" is complete")
        CI = sapply(order, function(x) MCIl[[x]][maxMCIms[[1]][[x]]])
        ln = order
    }
    if (any(is(CI) == "list")) {
        warning("changing NA CI score(s) to 0")
        idx = sapply(CI, function(x) length(x) == 0)
        CI[idx] = 0
        CI = do.call(c, CI)
        names(CI) = ln
    }
    if (!is.null(states)) {
        CI = CI[states]
        CI[is.na(CI)] = 0
        ln = names(CI) = states
    }
    matplot(CI, type = "l", ylab = "DNB score", axes = FALSE)
    len = sapply(ln, function(x) length(maxMCIms[[2]][[x]]))
    len[is.na(len)] = 0
    names(len) = ln
    text(seq_along(CI), CI + 0.01, paste0("(", len, ")"))
    axis(2)
    axis(side = 1, at = seq_along(CI), labels = ln, las = las)
}
plotIc <-
function (Ic, las = 0, order = NULL, ylab = "Ic.shrink", col = "black", 
    main = NULL, add = FALSE, ylim = NULL, lty = 1:5, lwd = 1) 
{
    if (!is.null(order)) {
        if (any(!order %in% names(Ic))) 
            stop("make sure \"Ic\" is named using names in \"order\"")
        if (any(!names(Ic) %in% order)) 
            warning("not every state in \"Ic\" is plotted,  make sure \"order\" is complete")
        Ic = Ic[order]
    }
    matplot(Ic, type = "l", ylab = ylab, axes = FALSE, col = col, 
        main = main, add = add, ylim = ylim, lty = lty, lwd = lwd)
    axis(2)
    stages = names(Ic)
    axis(side = 1, at = seq_along(Ic), labels = stages, las = las)
}
plot_Ic_Simulation <-
function (Ic, simulation, las = 0, ylim = NULL, order = NULL, 
    main = NULL, ylab = "Ic.shrink", fun = c("matplot", "boxplot"), 
    which2point = NULL) 
{
    fun <- match.arg(fun)
    if (any(is.null(names(Ic)))) 
        stop("Please provide name for vector \"Ic\" ")
    if (any(is.null(rownames(simulation)))) 
        stop("Please provide rowname for vectors of \"simulation\" ")
    if (length(Ic) != nrow(simulation)) 
        stop("Please provide the same length of \"Ic\" and vectors of \"simulation\" ")
    if (!identical(names(Ic), row.names(simulation))) 
        Ic = Ic[match(row.names(simulation), names(Ic))]
    if (!is.null(which2point)) {
        if (!(which2point %in% rownames(simulation) | which2point %in% 
            1:nrow(simulation))) 
            stop("which2point must be a state name of integer indicating the state of interested.")
    }
    if (fun == "matplot") {
        toplot = cbind(simulation, Ic)
    }
    else {
        toplot = simulation
    }
    if (!is.null(order)) {
        if (any(!names(Ic) %in% order)) 
            warning("not all states in Ic is plotted")
        if (any(!order %in% names(Ic))) 
            stop("make sure \"Ic\" is named using names in \"order\"")
        toplot = toplot[order, ]
    }
    if (is.null(ylim)) {
        if (is.null(which2point)) {
            ylim = c(min(c(Ic, simulation)), max(Ic, 2 * (max(simulation) - 
                min(simulation))))
        }
        else {
            ylim = c(min(c(Ic, simulation)), max(max(Ic[which2point], 
                2 * (max(simulation[which2point, ]) - min(simulation[which2point, 
                  ])), simulation)))
        }
    }
    if (fun == "matplot") {
        matplot(toplot, type = "l", col = c(rep("grey", ncol(toplot) - 
            1), "red"), lty = 1, ylab = ylab, axes = FALSE, ylim = ylim, 
            main = main)
    }
    else {
        boxplot(t(toplot), col = c(rep("grey", ncol(toplot) - 
            1), "red"), ylab = ylab, axes = FALSE, ylim = ylim, 
            main = main)
        points(Ic, col = "red", type = "b")
        x <- lapply(Ic, function(x) table(toplot > x))
        y <- unlist(lapply(x, function(X) X[2]/sum(X)))
        if (any(is.na(y))) {
            y[which(is.na(y))] = 0
        }
        sig <- which(y < 0.050000000000000003)
        if (length(sig) > 0) 
            mtext(round(y[sig], 3), line = -5, at = (1:length(Ic))[sig])
    }
    axis(2)
    stages = row.names(toplot)
    axis(side = 1, at = seq_along(stages), labels = stages, las = las)
    if (is.null(which2point)) {
        abline(h = min(simulation), col = "grey", lty = 3)
        abline(h = max(simulation), col = "grey", lty = 3)
        abline(h = min(simulation) + 2 * (max(simulation) - min(simulation)), 
            col = "grey", lty = 2)
    }
    else {
        abline(h = min(simulation[which2point, ]), col = "grey", 
            lty = 3)
        abline(h = max(simulation[which2point, ]), col = "grey", 
            lty = 3)
        abline(h = min(simulation[which2point, ]) + 2 * (max(simulation[which2point, 
            ]) - min(simulation[which2point, ])), col = "grey", 
            lty = 2)
    }
}
plot_MCI_Simulation <-
function (MCI, simulation, las = 0, order = NULL, ylim = NULL, 
    main = NULL, which2point = NULL, ...) 
{
    if (is.null(names(MCI))) 
        stop("make sure elements in \"MCI\" have names")
    if (!is.null(order)) {
        if (any(!order %in% row.names(simulation))) 
            stop("make sure \"simulation\" has row.names which are in \"order\"")
        if (any(!row.names(simulation) %in% order)) 
            warning("not every state in \"simulation\" is plotted,  make sure \"order\" is complete")
        simulation = simulation[order, ]
    }
    maxpt = max(simulation, MCI, na.rm = TRUE)
    tmp = c(min(simulation, MCI, na.rm = TRUE), maxpt)
    if (is.null(ylim)) {
        if (min(simulation, na.rm = TRUE) < maxpt) {
            ylim = tmp
        }
        else {
            ylim = rev(tmp)
        }
    }
    boxplot(t(simulation), col = "grey", ylab = "DNB score", 
        axes = FALSE, ylim = ylim, main = main, pch = 20, ...)
    x = which.max(MCI)
    maxCI = MCI[x]
    if (!is.null(order)) {
        if (is.null(names(MCI))) 
            stop("make sure \"MCI\" is named using names in \"order\"")
    }
    axis(2)
    if (is.null(order)) {
        stages = row.names(simulation)
    }
    else {
        stages = order
    }
    x = which(stages == names(x))
    axis(side = 1, at = 1:nrow(simulation), labels = stages, 
        las = las)
    points(x, maxCI, col = "red", pch = 16)
    if (is.null(which2point)) {
        abline(h = min(simulation), col = "grey", lty = 3)
        abline(h = max(simulation), col = "grey", lty = 3)
        abline(h = min(simulation) + 2 * (max(simulation) - min(simulation)), 
            col = "grey", lty = 2)
    }
    else {
        abline(h = min(simulation[which2point, ]), col = "grey", 
            lty = 3)
        abline(h = max(simulation[which2point, ]), col = "grey", 
            lty = 3)
        abline(h = min(simulation[which2point, ]) + 2 * (max(simulation[which2point, 
            ]) - min(simulation[which2point, ])), col = "grey", 
            lty = 2)
    }
}
plot_SS_Simulation <-
function (Ic, simulation, las = 0, xlim = NULL, ylim = NULL, 
    order = NULL, main = "1st max - 2nd max", ylab = "Density", 
    na.rm = TRUE) 
{
    if (any(is.null(names(Ic)))) 
        stop("Please provide name for vector \"Ic\" ")
    if (!identical(names(Ic), row.names(simulation))) 
        Ic = Ic[match(row.names(simulation), names(Ic))]
    if (!is.null(order)) {
        if (any(!names(Ic) %in% order)) 
            warning("not all states in Ic is plotted")
        if (any(!order %in% names(Ic))) 
            stop("make sure \"Ic\" is named using names in \"order\"")
        toplot = toplot[order, ]
    }
    diff_Ic <- apply(simulation, MARGIN = 2, function(x) sort(x, 
        decreasing = TRUE)[2:1])
    diff_Ic <- apply(diff_Ic, MARGIN = 2, diff)
    density_diff = density(diff_Ic, na.rm = na.rm)
    if (is.null(ylim)) 
        ylim = c(0, 0.10000000000000001 + max(density_diff$y))
    if (is.null(xlim)) {
        xlim = c(-0.050000000000000003, 0.050000000000000003) + 
            c(min(density_diff$x), max(density_diff$x))
    }
    plot(density_diff, type = "l", lwd = 2, col = "black", main = main, 
        ylab = paste(ylab, "Density"), xlim = xlim, ylim = ylim, 
        cex.main = 1.2, cex.lab = 1.2)
    v <- diff(sort(Ic, decreasing = TRUE)[2:1])
    abline(v = v, col = "red", lwd = 2, lty = 2)
    P.value <- round(mean(diff_Ic >= v), 3)
    legend("topright", legend = paste0("p = ", P.value), text.col = "red", 
        bty = "n", cex = 1.5)
    return(P.value)
}
plotIcSignificance <-
function (filename, BioTIP_scores, CTS.candidate, SimResults_g, 
    width = 10, height = 6, fixedylim = FALSE, nc = NULL) 
{
    require(BioTIP)
    n = length(BioTIP_scores)
    P_delta_Ic = array(dim = n)
    names(P_delta_Ic) = names(BioTIP_scores)
    if (is.null(names(CTS.candidate))) 
        stop("pls give a named list of \"CTS.candidate\" ")
    if (all(names(BioTIP_scores) %in% names(CTS.candidate))) 
        names(BioTIP_scores) = names(CTS.candidate)
    else print("names(BioTIP_scores) != names(CTS.candidate)")
    if (all(names(SimResults_g) %in% names(CTS.candidate))) 
        names(SimResults_g) = names(CTS.candidate)
    else print("names(SimResults_g) != names(CTS.candidate)")
    if (any(grepl(".", names(BioTIP_scores), fixed = TRUE))) {
        names(BioTIP_scores) = lapply(names(BioTIP_scores), function(x) unlist(strsplit(x, 
            split = ".", fixed = T))[1]) %>% unlist
        names(CTS.candidate) = lapply(names(CTS.candidate), function(x) unlist(strsplit(x, 
            split = ".", fixed = T))[1]) %>% unlist
        names(SimResults_g) = lapply(names(SimResults_g), function(x) unlist(strsplit(x, 
            split = ".", fixed = T))[1]) %>% unlist
    }
    plot_simulation <- function(i) {
        x = length(CTS.candidate[[i]])
        if (fixedylim) 
            ylim = c(0, max(unlist(BioTIP_scores)))
        else ylim = c(0, max(BioTIP_scores[[i]]))
        plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], 
            ylim = ylim, las = 2, ylab = "Ic.shrink", main = paste("Cluster", 
                names(BioTIP_scores)[i], "\n", x, "genes"), fun = "boxplot", 
            which2point = names(BioTIP_scores)[i])
        P = plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], 
            main = paste("Delta Ic.shrink", x, "genes"), ylab = NULL, 
            xlim = range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]], 
                SimResults_g[[i]])))
        return(P)
    }
    pdf(file = filename, width = width, height = height)
    if (is.null(nc)) 
        nc <- ifelse(n > 8, 8, n)
    newPageFLAG <- ifelse(n > 8, TRUE, FALSE)
    par(mfrow = c(2, nc))
    for (i in 1:nc) P_delta_Ic[1:nc] = plot_simulation(i)
    if (newPageFLAG) {
        if (n > 8 & n <= 16) {
            for (i in (nc + 1):n) P_delta_Ic[i] = plot_simulation(i)
        }
        if (n > 16 & n <= 24) {
            for (i in (nc + 1):(2 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (2 * nc + 1):n) P_delta_Ic[i] = plot_simulation(i)
        }
        if (n > 24 & n <= 32) {
            for (i in (nc + 1):(2 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (2 * nc + 1):(3 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (3 * nc + 1):n) P_delta_Ic[i] = plot_simulation(i)
        }
        if (n > 32 & n <= 40) {
            for (i in (nc + 1):(2 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (2 * nc + 1):(3 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (3 * nc + 1):(4 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (4 * nc + 1):n) P_delta_Ic[i] = plot_simulation(i)
        }
        if (n > 40 & n <= 48) {
            for (i in (nc + 1):(2 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (2 * nc + 1):(3 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (3 * nc + 1):(4 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (4 * nc + 1):(5 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (5 * nc + 1):n) P_delta_Ic[i] = plot_simulation(i)
        }
        if (n > 48 & n <= 56) {
            for (i in (nc + 1):(2 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (2 * nc + 1):(3 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (3 * nc + 1):(4 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (4 * nc + 1):(5 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (5 * nc + 1):(6 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (6 * nc + 1):n) P_delta_Ic[i] = plot_simulation(i)
        }
        if (n > 56 & n <= 64) {
            for (i in (nc + 1):(2 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (2 * nc + 1):(3 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (3 * nc + 1):(4 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (4 * nc + 1):(5 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (5 * nc + 1):(6 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (6 * nc + 1):(7 * nc)) P_delta_Ic[i] = plot_simulation(i)
            for (i in (7 * nc + 1):n) P_delta_Ic[i] = plot_simulation(i)
        }
        if (n > 64) 
            warning("only the first 65 cases were plotted")
    }
    dev.off()
    return(P_delta_Ic)
}
