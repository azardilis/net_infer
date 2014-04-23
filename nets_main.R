source("nets_functions.R")
dag0 <- NULL
dag0$parents <- list(numeric(0),numeric(0),c(1,2),c(2,3))
dag0$n <- length(dag0$parents)
dag0$levels <- c(2,2,2,2)
dag0$changed <- 1:4
dag0$loglik.families <- NULL

dag0$sorted <- topsort.parents(dag0$parents)

CPD1 <- array(c(0.5),dim=c(1))  # prob that node 1 is 1
CPD2 <- array(c(0.5),dim=c(1))  # prob that node 2 is 1
CPD3 <- array(c(0.1,0.1,
                0.2,0.9),dim=c(2,2,1))
CPD4 <- array(c(0.1,0.2,
               0.1,0.9),dim=c(2,2,1))

dag0$CPDs <- list(CPD1,CPD2,CPD3,CPD4)

arc.priors <- make.arc.priors(dag0$n,
                              list(c(3,4)),
                              list(c(0.99, 0.005, 0.005))
                              )

data <- generate.data.dag(dag0, 300)
CPD.prior.value <- 0.03
dag0$changed <- 1:(dag0$n)
dag0 <- loglik.dag(dag0, CPD.prior.value, data, arc.priors)
dag0$loglik


dag <- search.dag(data, 30, 30, CPD.prior.value, arc.priors)
dag$loglik; dag0$loglik






dat <- generate.data.dag(dag0, 200)
CPD.prior.value <- 0.03
dag0$changed <- 1:(dag0$n)
dag0 <- loglik.dag(dag0, CPD.prior.value, dat, arc.priors)

n.iter <- 2000
mcmc <- mcmc.dag(dat, n.iter, CPD.prior.value, arc.priors)

mean(mcmc$accepts)
logliks <- sapply(mcmc$dags,function(dag) dag$loglik)
plot(logliks,ty="l")
abline(h=dag0$loglik,col="red")


