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


# question 2
test.search <- function(dag0, n.datapoints) {
    data <- generate.data.dag(dag0, n.datapoints)
    CPD.prior.value <- 0.03
    dag0$changed <- 1:(dag0$n)
    dag0 <- loglik.dag(dag0, CPD.prior.value, data, arc.priors)
    dag <- search.dag(data, 30, 30, CPD.prior.value, arc.priors)

    return(dag$loglik == dag0$loglik)
}

arc.priors <- make.arc.priors(dag0$n,
                              list(c(3,4)),
                              list(c(0.99,0.005,0.005))
                              )


n.datapoints <- seq(10, 500, 50)
res1 <- lapply(n.datapoints, function(n) { replicate(100, test.search(dag0, n)) })
succ <- sapply(res, function(r) { mean(r*1) })


## question 3
dag2 <- NULL
dag2$parents <- list(numeric(0), numeric(0))
dag2$n <- length(dag1$parents)
dag2$levels <- c(2,2)
dag2$changed <- 1:2
dag2$loglik.families <- NULL

dag2$sorted <- topsort.parents(dag1$parents)

CPD1 <- array(c(0.5),dim=c(1))
CPD2 <- array(c(0.5), dim=c(1))

dag2$CPDs <- list(CPD1,CPD2)

#question 4
dag0 <- NULL
dag0$parents <- list(numeric(0), numeric(0), c(1, 2), c(2, 3), c(3, 4), c(5))
dag0$n <- length(dag0$parents)
dag0$levels <- rep(2, 6)
dag0$changed <- 1:6
dag0$sorted <- topsort.parents(dag0$parents)

CPD1 <- array(c(0.5),dim=c(1))  # prob that node 1 is 1
CPD2 <- array(c(0.5),dim=c(1))  # prob that node 2 is 1
CPD3 <- array(c(0.1,0.1,
                0.2,0.9),dim=c(2,2,1))
CPD4 <- array(c(0.1,0.2,
               0.1,0.9),dim=c(2,2,1))
CPD5 <- array(c(0.1,0.2,
               0.1,0.9),dim=c(2,2,1))
CPD6 <- matrix(c(0.4, 0.5), nrow=2, ncol=1)

dag0$CPDs <- list(CPD1,CPD2,CPD3,CPD4, CPD5, CPD6)

test.priors <- function(data, CPD.prior.value) {
    dag <- search.dag(data, 30, 30, CPD.prior.value)

    return(mean(sapply(dag$parents, function(par) { length(par) })))
}

data <- generate.data.dag(dag0, 300)
CPD.prior.values <- seq(0.0001, 10, 0.2)
res <- lapply(CPD.prior.values, function(pr) { replicate(50, test.priors(data, pr)) })
sapply(res, function(r) { mean(r) })















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






dag1 <- NULL
dag1$parents <- list(numeric(0), c(1))
dag1$n <- length(dag1$parents)
dag1$levels <- c(2,2)
dag1$changed <- 1:2
dag1$loglik.families <- NULL

dag1$sorted <- topsort.parents(dag1$parents)

CPD2 <- array(c(0.8),dim=c(1))
CPD1 <- matrix(c(0.4, 0.5), nrow=2, ncol=1)

dag1$CPDs <- list(CPD1,CPD2)

data <- generate.data.dag(dag1, 100)
dag1$changed <- 1:(dag1$n)
dag1 <- loglik.dag(dag1, 0.03, data, NULL)






