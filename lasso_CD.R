soft_threshold <- function(beta, lambda) {
    sign(beta) * max(abs(beta) - lambda, 0)
}

lasso_cd <- function(x, y, lambda, beta = NULL, tol = 1e-10, max_iter = 500,
                     patience = 10) {

    if(is.null(beta)) {
        beta = rep(0, ncol(x))
    }
    
    loss <- rep(0, max_iter)
    
    for (k in 1:max_iter) {
        residual <- y - x %*% beta
        loss[k] <- mean(residual * residual)
        
        for (j in 1:ncol(x)) {
            residual <- residual + x[, j] * beta[j]
            
            beta_ols_j <- mean(residual * x[, j])
            
            beta[j] <- soft_threshold(beta_ols_j, lambda)
            
            residual <- residual - x[, j] * beta[j]
        }
        
        if (k > patience) {
            if (loss[k] <= tol) {
                break
            }
        }
    }
    beta
}

######## logistic version ########

# calculate z, px, w each turn, then update betaj

getP <- function(X, betavec){
    Px <- 1/(1 + exp(-(X %*% betavec)))
    return(Px)
}

getW <- function(Px){
    W <- Px*(1-Px)
    return(Px)
}

getZ <- function(X, y, betavec, Px, W){
    Z <- X %*% betavec + (y - Px)/W
    return(Z)
}

lassoCD <- function(
        X, y, lambda = 0.5, init_beta = NULL, max_iter = 1e3, tol = 1e-5
){
    betavec <- init_beta
    betares <- init_beta
    N <- length(y)
    i <- 0
    cont <- TRUE
    while(i <= max_iter && cont){
        i <- i+1
        betaold <- betavec
        for(j in 1:length(betavec)){
            Px <- getP(X, betavec)
            W <- getW(Px)
            Z <- getZ(X, y, betavec, Px, W)
            betaresj <- betavec
            betaresj[j] <- 0
            Pres <- getP(X, betaresj)
            Wres <- getW(Pres)
            Zres <- getZ(X, y, betaresj, Pres, Wres)
            Wdiag <- diag(c(W), N, N)
            betaj <- (t(X[,j]) %*% Wdiag %*% (Z-Pres))
            betaj <- soft_threshold(betaj/N, lambda)
            betaj <- betaj/(t(X[,j]) %*% Wdiag %*% X[,j])
            betavec[j] <- betaj
        }
        if(sum(abs(betaold - betavec)) <= tol){
            cont <- FALSE
        }
        betares <- rbind(betares, betavec)
    }
    return(betares)
}

sim_logit_dat <- function(n, p, k, beta){
    x <- matrix(rnorm(n*p, 0, 1), n, p)
    betavec <- rep(0, p)
    betavec[1:k] <- beta
    p <- 1/(1+exp(-(x %*% betavec)))
    y <- rbinom(length(p), size = 1, prob = p)
    ret <- list(x, y)
    names(ret) <- c("X", "y")
    return(ret)
}

sim.dat <- sim_logit_dat(500, 10, 7, c(5,6,7,8,9,10,1))
sim.beta.res <- lassoCD(sim.dat$X,sim.dat$y, lambda = 0.02, init_beta = rep(3,10), max_iter = 1e3, tol = 1e-5)

#lassoCD(x1, y, lambda = 0.02, init_beta = rep(3,10), max_iter = 1e3, tol = 1e-5)

path <- function(x, y, grid){
  betas <- NULL
  for (i in 1:length(grid)){
    cor.result <- lassoCD(x, y, lambda = grid[i], init_beta = rep(3, dim(x)[2]), max_iter = 1e3, tol = 1e-5)
    lasbeta <- cor.result[nrow(cor.result),3:dim(cor.result)[2]]
    start <- lasbeta
    betas <- rbind(betas,c(lasbeta))
  }
  return(data.frame(cbind(grid,betas)))
}
path1 = path(x1, y, grid=exp(seq(-1,-10, length=100)))

path.plot <- path1 %>%
  gather(key = par, value = estimate, c(2:dim(path1)[2])) %>% 
  ggplot(aes(x = log(grid), y = estimate, group = par, col = par)) +
  geom_line()+
  ggtitle("path of solutions with a sequence of descending lambda's") +
  xlab("log(Lambda)") + 
  ylab("Estimate") +
  theme(legend.position = "right", legend.text = element_text(size = 6))

path.plot

