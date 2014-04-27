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
    dag0 <- loglik.dag(dag0, CPD.prior.value, data)
    dag <- search.dag(data, 30, 30, CPD.prior.value, NULL)

    return(dag$loglik == dag0$loglik)
}

arc.priors <- make.arc.priors(dag0$n,
                              list(c(3,4)),
                              list(c(0.99,0.005,0.005))
                              )


n.datapoints <- c(10, seq(50, 500, 50))
res <- lapply(n.datapoints, function(n) { replicate(100, test.search(dag0, n)) })
succ1<- sapply(res, function(r) { mean(r*1) })


## question 3
dag2 <- NULL
dag2$parents <- list(numeric(0), c(1))
dag2$n <- length(dag2$parents)
dag2$levels <- c(2,2)
dag2$changed <- 1:2
dag2$loglik.families <- NULL

dag2$sorted <- topsort.parents(dag2$parents)

CPD1 <- array(c(0.3),dim=c(1))
CPD2 <- matrix(c(0.5, 0.65), nrow=2, ncol=1)
#CPD2 <- array(c(0.5), dim=c(1))

dag2$CPDs <- list(CPD1,CPD2)
data <- generate.data.dag(dag2, 300)
dag2$changed <- 1:2
dag2 <- loglik.dag(dag2, 0.03, data)

dag3 <- NULL
dag3$parents <- list(numeric(0), numeric(0))
dag3$n <- length(dag3$parents)
dag3$levels <- c(2,2)
dag3$changed <- 1:2
dag3$loglik.families <- NULL

dag3 <- loglik.dag(dag3, 0.03, data)

test.CPD <- function(dag1, dag2) {
    data <- generate.data.dag(dag1, 300)
    dag1$changed <- 1:2
    dag1 <- loglik.dag(dag1, 0.03, data)
    dag2$changed <- 1:2
    dag2 <- loglik.dag(dag2, 0.03, data)

    if (dag1$loglik == dag2$loglik) { return(0) }
    else { return(which.max(c(dag1$loglik, dag2$loglik))) }

}

res <- replicate(100, test.CPD(dag2, dag3))

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
    dag <- search.dag(data, 30, 30, CPD.prior.value, NULL)

    return(mean(sapply(dag$parents, function(par) { length(par) })))
}

data <- generate.data.dag(dag0, 300)
CPD.prior.values <- c(0.0001, seq(0.2, 10, 0.2))
res <- lapply(CPD.prior.values, function(pr) { replicate(20, test.priors(data, pr)) })
sapply(res, function(r) { mean(r) })


#question 5
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


#question 7
dag <- NULL
dag$parents <- list(c(2), numeric(0), c(1, 2))
dag$n <- length(dag$parents)
dag$levels <- c(2,2,2)
dag$changed <- 1:3
dag$loglik.families <- NULL

dag$sorted <- topsort.parents(dag$parents)

#CPD1 <- array(c(0.5),dim=c(1))  # prob that node 1 is 1
CPD1 <- matrix(c(0.2, 0.7), nrow=2, ncol=1)
CPD2 <- array(c(0.5),dim=c(1))  # prob that node 2 is 1
CPD3 <- array(c(0.1,0.1,
                0.2,0.9),dim=c(2,2,1))

dag$CPDs <- list(CPD1,CPD2, CPD3)

data <- generate.data.dag(dag, 300)

test.search <- function(dag, data) {
    inf.dag <- search.dag(data, 30, 30, 0.03, NULL)
    return(inf.dag$loglik == dag$loglik)
}

dag <- loglik.dag(dag, 0.03, data)
res <- replicate(100, test.search(dag, data))












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












