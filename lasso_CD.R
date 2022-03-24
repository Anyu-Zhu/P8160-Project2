soft_threshold <- function(beta, lambda) {
    sign(beta) * max(abs(beta) - lambda, 0)
}

lasso_cd <- function(x, y, lambda, beta = NULL, tol = 1e-10, max_iter = 5000,
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
        X, y, lambda = 0.5, init_beta = NULL, max_iter = 1e4, tol = 1e-8
){
    betavec <- init_beta
    N <- length(y)
    i <- 0
    loss <- 1e5
    prevloss <- Inf
    res <- c(0, loss, betavec)
    cont <- TRUE
    while(i <= max_iter && cont){
        i <- i + 1
        prevloss <- loss
        for(j in 1:length(betavec)){
            Px <- getP(X, betavec)
            W <- getW(Px)
            W <- ifelse(abs(W-0) < 1e-5, 1e-5, W)
            Z <- getZ(X, y, betavec, Px, W)
            betaresj <- betavec
            betaresj[j] <- 0
            Zresj <- X %*% betaresj
            betaj <- 
                soft_threshold(mean(W * X[,j] * (Z - Zresj)), lambda)/mean(W * X[,j] * X[,j])
            betavec[j] <- betaj
            loss <- (1/(2*N))*sum(W * (Z - X %*% betavec)^2) + lambda * sum(abs(betavec))
        }
        print(loss)
        if(abs(prevloss - loss) < tol || loss > Inf){
            cont <- FALSE
        }
        res <- rbind(res, c(i, loss, betavec))
    }
    return(res)
}

## matrix version ##
lassoCD.matrix <- function(
        X, y, lambda = 0.5, init_beta = NULL, max_iter = 1e4, tol = 1e-8
){
    betavec <- init_beta
    N <- length(y)
    i <- 0
    loss <- 1e5
    prevloss <- Inf
    res <- c(0, loss, betavec)
    cont <- TRUE
    while(i <= max_iter && cont){
        i <- i + 1
        prevloss <- loss
        for(j in 1:length(betavec)){
            Px <- getP(X, betavec)
            W <- getW(Px)
            W <- ifelse(abs(W-0) < 1e-5, 1e-5, W)
            Wdiag <- diag(c(W), N, N)
            Z <- getZ(X, y, betavec, Px, W)
            betaresj <- betavec
            betaresj[j] <- 0
            Zresj <- X %*% betaresj
            betaj <- 
                soft_threshold((1/N) * (t(X[,j]) %*% Wdiag %*% (Z-(X %*% betaresj))),
                               lambda)/(1/N)*(t(X[,j]) %*% Wdiag %*% X[,j])
            betavec[j] <- betaj
            loss <- (1/(2*N))*sum(Wdiag %*% (Z - X %*% betavec)^2) + lambda * sum(abs(betavec))
        }
        print(loss)
        if(abs(prevloss - loss) < tol || loss > Inf){
            cont <- FALSE
        }
        res <- rbind(res, c(i, loss, betavec))
    }
    return(res)
}

##########
coordinatelasso <- function(lambda, dat, s, tol=1e-10, maxiter = 200){
    i <- 0 
    pp <- length(s)
    n <- length(dat$y)
    betavec <- s
    loglik <- 1e6
    res <- c(0, loglik, betavec)
    prevloglik <- Inf # To make sure it iterates 
    while (i < maxiter && abs(loglik - prevloglik) > tol && loglik < Inf) {
        i <- i + 1 
        prevloglik <- loglik
        for (j in 1:pp) {
            u <- dat$X %*% betavec
            expu <- exp(u) 
            prob <- expu/(expu+1)
            w <- prob*(1-prob) # weighted
            # avoid coeffcients diverging in order to achieve fitted  probabilities of 0 or 1.
            w <- ifelse(abs(w-0) < 1e-5, 1e-5, w)
            z <- u + (dat$y-prob)/w
            # calculate noj
            znoj <- dat$X[,-j] %*% betavec[-j]
            # revise the formula to be z
            betavec[j] <- soft_threshold(mean(w*(dat$X[,j])*(z - znoj)), lambda)/(mean(w*dat$X[,j]*dat$X[,j]))
        }
        loglik <- sum(w*(z-dat$X %*% betavec)^2)/(2*n) + lambda * sum(abs(betavec))
        res <- rbind(res, c(i, loglik, betavec))}  
    return(res)
}

dat = list(X = data.cancer, y = cancer$diagnosis)
coordlasso.res <- coordinatelasso(0.5, dat, rep(0,31))
##########

library(tidyverse)
cancer = read.csv("breast-cancer.csv") %>% 
    mutate(diagnosis = factor(diagnosis))

mean(cancer$radius_mean)
sqrt(var(cancer$radius_mean))

cancer["intercept"] = 1
cancer$diagnosis = case_when(cancer$diagnosis == "M" ~ 1,
                             TRUE ~ 0)
data.cancer = data.matrix(cancer[3:33])

cancer.lassoCD <- lassoCD(data.cancer, cancer$diagnosis, init_beta = c(rep(0,31)), lambda = 0.5)
cancer.lassoCD.mtx <- lassoCD.matrix(data.cancer, cancer$diagnosis, init_beta = c(rep(20,31)), lambda = 0.5)
lassoCD.predict <- function(betavec, X_new, y){
    Py <- 1/(1 + exp(-(X_new %*% betavec)))
    Pn <- 1-Py
    res <- list(P_pos = Py, P_neg = Pn, res = Py > Pn)
    res$res = as.numeric(res$res)
    return(res)
}

cancer.lassoCD.predict <- lassoCD.predict(cancer.lassoCD[6,3:33],data.cancer, cancer$diagnosis)

pred.obs <- data.frame(pred = cancer.lassoCD.predict$res, obs = cancer$diagnosis)
pred.obs$mis = pred.obs$pred == pred.obs$obs
sum(pred.obs$mis)

library(glmnet)
glmnet.cancer = glmnet(data.cancer, cancer$diagnosis, family = "binomial", lambda = 0, intercept = FALSE)
predict.glmnet(glmnet.cancer, s = 0.5, type = "coefficient")


library(pROC)
roc(pred.obs$pred, pred.obs$obs)



###################################################
#                 for path                        #
###################################################
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





