

perturb <- function(a, b) {
    p <- length(a)
    out <- rep(NA, p)
    for (i in 1:p) {
        out[i] <- a[i] * b[i] / t(a) %*% b
    }
    return(out)
}


normalize <- function(X,cen=TRUE,scale=TRUE) {
    # centering
    n <- nrow(X)
    if(cen)
        X <- apply(X, 2, function(x) x - mean(x))
    if(scale)
        X <- apply(X, 2, function(x) sqrt(n) * x / sd(x))
    return(X)
}

normalize_X <- function(X){
    X <- log(X)
    X <- normalize(X,scale=F)
    X <- exp(X)
    return(X)
}


get_data <- function(X,y,y1, y_mode,min_perturb,max_perturb, K=5){
    Z  <- X
    for(i in 1:nrow(X)){
        U <- runif(p, min = min_perturb, max = max_perturb) ## measurement error factors
        Z[i, ] <- perturb(X[i, ], U)
    }
    if(y_mode==1)
        y <- as.matrix(y) ## BMI
    else if(y_mode==2)
        y <- as.matrix(y1) ## cov-adjusted BMI
    X <- log(as.matrix(X))
    Z <- log(as.matrix(Z))
    B <- Z - X
    Sig_B <- t(B) %*% B / nrow(B)
    return(list(X=X,Z=Z,y=y,Sig_B=Sig_B))
}

get_norm_data <- function(X,y,y1, y_mode,min_perturb,max_perturb, K=5){
    Z  <- X
    X <- normalize_X(X)
    for(i in 1:nrow(X)){
        U <- runif(p, min = min_perturb, max = max_perturb) ## measurement error factors
        # U_c <- U * runif(p, min = 0.5, max = 2)
        Z[i, ] <- perturb(X[i, ], U)
    }
    if(y_mode==1)
        y <- as.matrix(y) ## BMI
    else if(y_mode==2)
        y <- as.matrix(y1) ## cov-adjusted BMI
    X <- log(as.matrix(X))
    Z <- log(as.matrix(Z))
    B <- Z - X
    Sig_B <- t(B) %*% B / nrow(B)
    return(list(X=X,Z=Z,y=y,Sig_B=Sig_B))
}