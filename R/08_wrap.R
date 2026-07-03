BioTIP.wrap <-
function (sce, samplesL, subDir = "newrun", smallest.population.size = 20, 
    globa.HVG.select = FALSE, dec.pois = NULL, n.getTopHVGs = round(nrow(sce) * 
        0.40000000000000002), getTopMCI.n.states = ncol(sce), 
    getTopMCI.gene.minsize = 30, getTopMCI.gene.maxsize = NULL, 
    MCIbottom = 2, local.HVG.optimize = TRUE, localHVG.preselect.cut = 0.10000000000000001, 
    localHVG.runs = 100, logmat.local.HVG.testres = NULL, getNetwork.cut.fdr = 0.050000000000000003, 
    n.getMaxMCImember = 2, getMCI.adjust.size = FALSE, n.permutation = 100, 
    empirical.MCI.p.cut = 0.01, empirical.IC.p.cut = 0.050000000000000003, 
    local.IC.p = TRUE, IC.rank = 1, M = NULL, permutation.method = c("gene", 
        "both", "sample"), verbose = FALSE, plot = TRUE) 
{
    assay_names <- assayNames(sce)
    if (!"logcounts" %in% assay_names & !"counts" %in% assay_names) 
        stop("no 'counts' nor logcounts' in names(assays(sce))")
    if (!permutation.method %in% c("gene", "both", "sample")) {
        permutation.method = "gene"
        warning("No permutation.method found for Ic signicance, gene permutation is performed")
    }
    mainDir = getwd()
    if (!file.exists(subDir)) {
        dir.create(file.path(mainDir, subDir))
    }
    outputpath = paste0(subDir, "/BioTIP_top", localHVG.preselect.cut, 
        "FDR", getNetwork.cut.fdr, "_")
    (tmp = lengths(samplesL))
    if (any(tmp < smallest.population.size)) 
        samplesL <- samplesL[-which(tmp < smallest.population.size)]
    sce = sce[, unlist(samplesL)]
    myplotIc <- function(filename, BioTIP_scores, CTS.candidate, 
        SimResults_g) {
        pdf(file = filename, width = 10, height = 6)
        n = length(BioTIP_scores)
        nn <- ifelse(n > 8, 8, n)
        par(mfrow = c(2, nn))
        for (i in 1:nn) {
            x = length(CTS.candidate[[i]])
            plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], 
                ylim = c(0, max(unlist(BioTIP_scores))), las = 2, 
                ylab = "Ic.shrink", main = paste("Cluster", names(BioTIP_scores)[i], 
                  "\n", x, " genes"), fun = "boxplot", which2point = names(BioTIP_scores)[i])
            TEXT = ifelse(local.IC.p, "p.Local=", "p.Global=")
            text(1, 0.089999999999999997, paste(TEXT, p.IC[i]), 
                cex = 1.5)
        }
        for (i in 1:nn) {
            x = length(CTS.candidate[[i]])
            plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], 
                main = paste("Delta Ic.shrink", x, "genes"), 
                ylab = NULL, xlim = range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]], 
                  SimResults_g[[i]])))
        }
        if (n > 8) {
            for (i in (nn + 1):n) {
                x = length(CTS.candidate[[i]])
                plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], 
                  ylim = c(0, max(unlist(BioTIP_scores))), las = 2, 
                  ylab = "Ic.shrink", main = paste("Cluster", 
                    names(BioTIP_scores)[i], "\n", x, " genes"), 
                  fun = "boxplot", which2point = names(BioTIP_scores)[i])
                TEXT = ifelse(local.IC.p, "p.Local=", "p.Global=")
                text(1, 0.089999999999999997, paste(TEXT, p.IC[i]), 
                  cex = 1.5)
            }
            for (i in (nn + 1):n) {
                x = length(CTS.candidate[[i]])
                plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], 
                  main = paste("Delta Ic.shrink", x, "genes"), 
                  ylab = NULL, xlim = range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]], 
                    SimResults_g[[i]])))
            }
        }
        dev.off()
    }
    if (globa.HVG.select) {
        if (is.null(dec.pois)) {
            dec.pois <- scran::modelGeneVarByPoisson(sce, block = sce$batch)
        }
        else dec.pois = dec.pois[rownames(sce), ]
        libsf.var <- getTopHVGs(dec.pois, n = n.getTopHVGs)
        length(libsf.var)
        table(libsf.var %in% rownames(sce))
        libsf.var <- intersect(libsf.var, rownames(sce))
        length(libsf.var)
        dat <- sce[libsf.var, ]
    }
    else dat <- sce
    if (grepl("cell_data_set", class(dat))) {
        if (!requireNamespace("monocle3", quietly = TRUE)) 
            stop("BioTIP.wrap requires the monocle3 package for cell_data_set input")
        logmat <- as.matrix(getExportedValue("monocle3", "normalized_counts")(dat))
    }
    else if ((grepl("SingleCellExperiment", class(dat)))) 
        logmat <- as.matrix(logcounts(dat))
    dim(logmat)
    rm(dat)
    if (local.HVG.optimize) {
        if (is.null(logmat.local.HVG.testres)) {
            cat("running local HVG optimalization ...")
            testres <- optimize.sd_selection(logmat, samplesL, 
                B = localHVG.runs, cutoff = localHVG.preselect.cut, 
                times = 0.75, percent = 0.80000000000000004)
            if (verbose) 
                save(testres, file = paste0(outputpath, "optimized_local.HVG_selection.RData"), 
                  compress = TRUE)
        }
        else {
            testres <- logmat.local.HVG.testres
            if (!all(rownames(testres) %in% rownames(logmat))) {
                testres <- testres[intersect(rownames(testres), 
                  rownames(logmat)), ]
                warning("Some of genes in the given logmat.local.HVG.testres are outside the golbal HVG!")
            }
        }
    }
    else {
        testres <- sd_selection(logmat, samplesL, cutoff = localHVG.preselect.cut)
    }
    names(testres)
    cat("HVG selection is done", "\n")
    igraphL <- getNetwork(testres, fdr = getNetwork.cut.fdr)
    cluster <- getCluster_methods(igraphL)
    names(cluster)
    cat("Network partition is done", "\n")
    membersL <- getMCI(cluster, testres, adjust.size = getMCI.adjust.size, 
        fun = "BioTIP")
    names(membersL)
    cat("CTS score calculation is done", "\n")
    topMCI = getTopMCI(membersL[["members"]], membersL[["MCI"]], 
        membersL[["MCI"]], min = getTopMCI.gene.minsize, n = getTopMCI.n.states)
    topMCI
    if (class(topMCI) == "numeric" & length(topMCI) > 0) {
        if (plot) {
            w <- ifelse(length(membersL$MCI) > 10, 80, 40)
            h <- ifelse(length(membersL$MCI) > 10, 10, 5)
            pdf(file = paste0(outputpath, "MCIBar_", localHVG.preselect.cut, 
                "_fdr", getNetwork.cut.fdr, "_minsize", getTopMCI.gene.minsize, 
                ".pdf"), width = w, height = h)
            plotBar_MCI(membersL, ylim = c(0, ceiling(max(topMCI, 
                na.rm = TRUE))), minsize = getTopMCI.gene.minsize)
            if (!is.null(MCIbottom)) 
                abline(h = MCIbottom, lty = 2, col = "red")
            plotBar_MCI(membersL, ylim = c(0, ceiling(max(topMCI, 
                na.rm = TRUE)) * 2))
            if (!is.null(MCIbottom)) 
                abline(h = MCIbottom, lty = 2, col = "red")
            dev.off()
        }
        x <- which(topMCI >= MCIbottom)
        if (length(x) > 0) {
            topMCI = topMCI[x]
            if (getTopMCI.n.states > length(x)) 
                warning(paste("less number of states have the highest CTS score larger than", 
                  MCIbottom))
            getTopMCI.n.states = length(x)
            CTS.candidate.ms <- getMaxMCImember(membersL[["members"]], 
                membersL[["MCI"]], minsize = getTopMCI.gene.minsize, 
                n = n.getMaxMCImember)
            names(CTS.candidate.ms)
            CTS.candidate = getCTS(topMCI, CTS.candidate.ms[["members"]][names(topMCI)])
            if (n.getMaxMCImember > 1) {
                tmp <- unlist(lapply(CTS.candidate.ms[["idx"]][names(topMCI)], 
                  length))
                (whoistop2nd <- names(tmp[tmp >= 2]))
                if (length(whoistop2nd) > 0) {
                  nextMCI = getNextMaxStats(membersL[["MCI"]], 
                    idL = CTS.candidate.ms[["idx"]], whoistop2nd, 
                    which.next = 2)
                  nextMCI = nextMCI[order(nextMCI, decreasing = TRUE)]
                  x <- which(nextMCI >= MCIbottom)
                  if (length(x) > 0) {
                    nextMCI = nextMCI[x]
                    whoistop2nd = names(nextMCI)
                    CTS.candidate = append(CTS.candidate, CTS.candidate.ms[["2topest.members"]][whoistop2nd])
                    topMCI = append(topMCI, nextMCI)
                    if (n.getMaxMCImember > 2) {
                      whoistop3rd = names(tmp[tmp >= 3])
                      if (length(whoistop3rd) > 0) {
                        nextMCI = getNextMaxStats(membersL[["MCI"]], 
                          idL = CTS.candidate.ms[["idx"]], whoistop3rd, 
                          which.next = 3)
                        nextMCI = nextMCI[order(nextMCI, decreasing = TRUE)]
                        x <- which(nextMCI >= MCIbottom)
                        if (length(x) > 0) {
                          CTS.candidate = append(CTS.candidate, 
                            CTS.candidate.ms[["3topest.members"]][whoistop3rd])
                          topMCI = append(topMCI, nextMCI)
                          if (n.getMaxMCImember > 3) 
                            warning("Please manually modify the BioTIP() accordingly to extract 4th and more CTS candididate per cell cluster")
                        }
                      }
                    }
                  }
                }
            }
            else warning(paste("All 2nd highest CTS scores are lower than", 
                MCIbottom))
            if (verbose) 
                save(CTS.candidate, topMCI, file = paste0(outputpath, 
                  "CTS.candidate.RData"))
            rm(sce)
            dim(logmat)
            if (is.null(M)) {
                cat("calcualting M .... ")
                M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, 
                  shrink = TRUE)
                save(M, file = paste0(outputpath, "CTS_ShrinkM.RData"), 
                  compress = TRUE)
            }
            dim(M)
            M = M[rownames(logmat), rownames(logmat)]
            simuMCI = list()
            for (i in 1:length(CTS.candidate)) {
                n <- length(CTS.candidate[[i]])
                simuMCI[[i]] <- simulationMCI(n, samplesL, logmat, 
                  adjust.size = getMCI.adjust.size, B = n.permutation, 
                  fun = "BioTIP", M = M)
            }
            names(simuMCI) = names(CTS.candidate)
            if (verbose) 
                save(simuMCI, file = paste0(outputpath, "SimuMCI_", 
                  n.permutation, "_", localHVG.preselect.cut, 
                  "_fdr", getNetwork.cut.fdr, "_minsize", getTopMCI.gene.minsize, 
                  ".RData"))
            n = length(CTS.candidate)
            if (plot & n > 0 & class(topMCI) == "numeric" & length(topMCI) > 
                0) {
                w <- ifelse(n > 3, 4 * n, 7)
                h <- ifelse(n > 3, n, 7)
                pdf(file = paste0(outputpath, "barplot_MCI_Sim_RandomGene.pdf"), 
                  width = w, height = h)
                par(mfrow = c(1, n))
                for (i in 1:n) {
                  plot_MCI_Simulation(topMCI[i], simuMCI[[i]], 
                    las = 2, ylim = c(0, max(c(topMCI, simuMCI[[i]]), 
                      na.rm = TRUE)), main = paste("Cluster", 
                      names(CTS.candidate)[i], ";", length(CTS.candidate[[i]]), 
                      "genes", "\n", "vs. ", n.permutation, "times of gene-permutation"), 
                    which2point = names(CTS.candidate)[i])
                }
                dev.off()
            }
            n <- length(CTS.candidate)
            p <- array(dim = n)
            for (i in 1:n) {
                p[i] = length(which(simuMCI[[i]] >= topMCI[i]))/n.permutation/length(samplesL)
            }
            dropoff <- which(p >= empirical.MCI.p.cut)
            if (any(dropoff)) {
                CTS.candidate <- CTS.candidate[-dropoff]
            }
            cat("CTS.dandidate is done", "\n")
            if (length(CTS.candidate) > 0) {
                BioTIP_scores <- SimResults_g <- SimResults_s <- SimResults_b <- list()
                for (i in 1:length(CTS.candidate)) {
                  CTS <- CTS.candidate[[i]]
                  n <- length(CTS)
                  BioTIP_scores[[i]] <- getIc(logmat, samplesL, 
                    CTS, fun = "BioTIP", shrink = TRUE, PCC_sample.target = "none")
                  if (permutation.method == "both") {
                    SimResults_b[[i]] <- matrix(nrow = length(samplesL), 
                      ncol = n.permutation)
                    rownames(SimResults_b[[i]]) = names(samplesL)
                    for (j in 1:length(samplesL)) {
                      ns <- length(samplesL[[j]])
                      CTS.sim <- sample(rownames(logmat), n)
                      SimResults_b[[i]][j, ] <- simulation_Ic_sample(logmat, 
                        ns, genes = CTS.sim, B = n.permutation, 
                        fun = "BioTIP", shrink = TRUE, PCC_sample.target = "none")
                    }
                  }
                  else {
                    if (permutation.method == "sample") {
                      SimResults_s[[i]] <- matrix(nrow = length(samplesL), 
                        ncol = n.permutation)
                      rownames(SimResults_s[[i]]) = names(samplesL)
                      for (j in 1:length(samplesL)) {
                        ns <- length(samplesL[[j]])
                        SimResults_s[[i]][j, ] <- simulation_Ic_sample(logmat, 
                          ns, genes = CTS, B = n.permutation, 
                          fun = "BioTIP", shrink = TRUE, PCC_sample.target = "none")
                      }
                    }
                    else {
                      SimResults_g[[i]] <- simulation_Ic(n, samplesL, 
                        logmat, B = n.permutation, fun = "BioTIP", 
                        shrink = TRUE, PCC_sample.target = "none")
                    }
                  }
                }
                names(BioTIP_scores) <- names(CTS.candidate)
                if (length(SimResults_g) > 0) 
                  names(SimResults_g) <- names(BioTIP_scores)
                if (length(SimResults_s) > 0) 
                  names(SimResults_s) <- names(BioTIP_scores)
                if (length(SimResults_b) > 0) 
                  names(SimResults_b) <- names(BioTIP_scores)
                if (verbose) 
                  if (permutation.method == "both") {
                    save(SimResults_b, BioTIP_scores, file = paste0(outputpath, 
                      "IC_sim.PermutateBoth.RData"), compress = TRUE)
                  }
                  else {
                    save(SimResults_g, SimResults_s, BioTIP_scores, 
                      file = paste0(outputpath, "IC_sim.Permutation.RData"), 
                      compress = TRUE)
                  }
                n <- length(CTS.candidate)
                p.IC <- rep(1, n)
                if (length(SimResults_g) > 0) {
                  SimResults_4_p.IC = SimResults_g
                }
                else if (length(SimResults_s) > 0) {
                  SimResults_4_p.IC = SimResults_s
                }
                else if (length(SimResults_b) > 0) {
                  SimResults_4_p.IC = SimResults_b
                }
                if (length(SimResults_4_p.IC) > 0) {
                  for (i in 1:n) {
                    interesting = names(BioTIP_scores[i])
                    p = length(which(SimResults_4_p.IC[[i]][interesting, 
                      ] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
                    p = p/n.permutation
                    p2 = length(which(SimResults_4_p.IC[[i]] >= 
                      BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
                    p2 = p2/n.permutation
                    p.IC[i] = p2/nrow(SimResults_4_p.IC[[i]])
                  }
                }
                cat("BioTIP is done", "\n")
                if (plot) {
                  if (length(SimResults_g) > 0 & length(CTS.candidate) > 
                    0) {
                    pdf.file.name = paste0(outputpath, "IC_Delta_SimresultGene.pdf")
                    myplotIc(pdf.file.name, BioTIP_scores, CTS.candidate, 
                      SimResults_g)
                  }
                  if (length(SimResults_s) > 0) {
                    pdf.file.name = paste0(outputpath, "IC_Delta_SimresultSample.pdf")
                    myplotIc(pdf.file.name, BioTIP_scores, CTS.candidate, 
                      SimResults_s)
                  }
                  if (length(SimResults_b) > 0) {
                    pdf.file.name = paste0(outputpath, "IC_Delta_SimresultBoth.pdf")
                    myplotIc(pdf.file.name, BioTIP_scores, CTS.candidate, 
                      SimResults_b)
                  }
                }
                n = length(BioTIP_scores)
                x <- p.IC < empirical.IC.p.cut
                x2 <- rep(FALSE, n)
                for (i in 1:n) {
                  if (IC.rank == 1) 
                    x2[i] = ifelse(names(which.max(BioTIP_scores[[i]])) == 
                      names(BioTIP_scores)[i], TRUE, FALSE)
                  else {
                    n2 <- length(BioTIP_scores[[i]])
                    which.high <- rank(BioTIP_scores[[i]])[n2:(n2 - 
                      IC.rank + 1)]
                    x2[i][which.high] <- TRUE
                  }
                }
                significant = x & x2
                names(significant) <- names(CTS.candidate)
            }
            else {
                BioTIP_scores = NULL
                significant = rep(FALSE, length(CTS.candidate))
                names(significant) <- names(CTS.candidate)
                cat("All CTS.candidate are insignificant in DMB model")
            }
            if (!is.null(getTopMCI.gene.maxsize)) {
                x <- lengths(CTS.candidate)
                dropoff <- which(x > getTopMCI.gene.maxsize)
                if (length(dropoff) > 0) {
                  CTS.candidate = CTS.candidate[-dropoff]
                  topMCI = topMCI[-dropoff]
                  BioTIP_scores = BioTIP_scores[-dropoff]
                  significant = significant[-dropoff]
                }
            }
            return(list(CTS.candidate = CTS.candidate, CTS.score = topMCI, 
                Ic.shrink = BioTIP_scores, significant = significant))
        }
        else {
            warning(paste("No CTS has a score higher than", MCIbottom, 
                "!"))
            return(CTS.candidate = NULL)
        }
    }
    else {
        warning(paste("No module has", getTopMCI.gene.minsize, 
            " gene members!"))
        return(CTS.candidate = NULL)
    }
}
