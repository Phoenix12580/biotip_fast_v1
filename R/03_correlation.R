avg.cor.shrink <-
function (X, Y = NULL, MARGIN = c(1, 2), shrink = TRUE, abs = FALSE, 
    target = 0) 
{
    if (target == "none") 
        shrink = FALSE
    if (class(target) == "numeric") 
        if ((target < 0) | (target > 1)) {
            stop("Argument `target` must be a value between 0 and 1, or a choice of 'zero',  'average', 'half', 'none'")
        }
        else if (class(target) == "character") 
            if (!target %in% c("none", "zero", "average", "half")) {
                stop("Argument `target` must be a value between 0 and 1, or a choice of 'zero',  'average', 'half','none'")
            }
    if (MARGIN != 1 & MARGIN != 2) 
        stop("MARGIN must be a choice of 1 or 2.")
    X_cor_shrink = cor.shrink(X = X, Y = Y, MARGIN = MARGIN, 
        shrink = shrink, target = target)
    if (is.null(Y)) {
        U <- upper.tri(X_cor_shrink, diag = FALSE)
        X_cor_shrink <- X_cor_shrink[U]
    }
    if (abs == TRUE) {
        res = mean(abs(X_cor_shrink), na.rm = TRUE)
    }
    else {
        res = mean(X_cor_shrink, na.rm = TRUE)
    }
    return(res)
}
cor.shrink <-
function (X, Y = NULL, MARGIN = c(1, 2), shrink = TRUE, target = 0) 
{
    if (target == "average") 
        shrink = FALSE
    if (target == "none") 
        shrink = FALSE
    if (class(target) == "numeric") 
        if ((target < 0) | (target > 1)) {
            stop("Argument `target` must be a value between 0 and 1, or a choice of 'zero',  'average', 'half', 'none'")
        }
        else if (class(target) == "character") 
            if (!target %in% c("zero", "average", "half", "none")) {
                stop("Argument `target` must be a value between 0 and 1, or a choice of 'zero',  'average', 'half', 'none")
            }
    dim_X = dim(X)
    dim_Y = dim(Y)
    Y.exist = FALSE
    if (!is.null(Y)) {
        if (MARGIN == 1) {
            X <- rbind(X, Y)
            Y.exist = TRUE
        }
        else {
            X <- cbind(X, Y)
        }
    }
    X_means = apply(X, MARGIN = MARGIN, mean, na.rm = TRUE)
    X_sds = apply(X, MARGIN = MARGIN, sd, na.rm = TRUE)
    X_sds[X_sds == 0] = 1
    X_std = sweep(sweep(X, MARGIN = MARGIN, STATS = X_means, 
        FUN = "-"), MARGIN = MARGIN, STATS = X_sds, FUN = "/")
    rm(Y)
    Y = !is.na(X_std)
    X_std[!Y] = 0
    if (MARGIN == 1) {
        XtX = Matrix::tcrossprod(X_std)
        X2tX2 = Matrix::tcrossprod(X_std^2)
        YtY = Matrix::tcrossprod(Y)
    }
    else {
        XtX = Matrix::crossprod(X_std)
        X2tX2 = Matrix::crossprod(X_std^2)
        YtY = Matrix::crossprod(Y)
    }
    X_cor = XtX/(pmax(1, YtY - 1, na.rm = TRUE))
    X_cor_shrink = X_cor
    U = upper.tri(XtX, diag = FALSE)
    if (shrink) {
        numerator = sum(((YtY[U] * X2tX2[U]) - (XtX[U])^2)/((pmax(1, 
            YtY[U] - 1))^3), na.rm = TRUE)
        if (class(target) == "character") {
            if (target == "zero") {
                target = 0
            }
            else {
                if (target == "half") {
                  target = 0.5
                }
                else {
                  target = mean(X_cor[U])
                }
            }
        }
        denominator = sum((X_cor[U] - target)^2)
        lambda = ifelse(shrink, max(0, min(1, numerator/denominator), 
            na.rm = TRUE), 0)
        target_cor = matrix(target, nrow = nrow(X_cor), ncol = ncol(X_cor))
        diag(target_cor) = 1
        X_cor_shrink = (lambda * target_cor) + ((1 - lambda) * 
            X_cor)
    }
    return(X_cor_shrink)
}
