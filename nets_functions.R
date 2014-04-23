topsort.parents <- function(parents) {
  ## parents: list of parent vectors
  ## out: vector of topsorted nodes, NULL if cyclic
  k <- length(parents)
  children <- vector("list",k)
  for (i in 1:k) {
    par <- parents[[i]]
    children[par] <- lapply(children[par],"c",i)
  }
  indeg <- sapply(parents,length)
  indeg0 <- which(indeg == 0)
  sorted <- NULL
  while (length(indeg0) > 0) {
    sorted <- c(sorted,indeg0)
    indeg[indeg0] <- -1
    for (n in indeg0) {
      ch <- children[[n]]
      indeg[ch] <- indeg[ch] - 1
    }
    indeg0 <- which(indeg == 0)
  }
  if (any(indeg > 0)) NULL # cyclic
else sorted
}


generate.data.dag <- function(dag,n) {
  ## dag, n: sample size
  ## out: data frame n x dag$n of samples
  m <- matrix(0,ncol=dag$n,nrow=n)
  for (i in 1:n) {
    m[i,] <- generate.sample.dag(dag)
  }
  ## convert columns to factors:
  df <- data.frame(m)
  for (i in 1:ncol(m)) {
    df[,i] <- factor(m[,i],levels=1:dag$levels[i])
  }
  df
}

generate.sample.dag <- function(dag) {
    x <- rep(0, dag$n)
    for (i in 1:dag$n) {
        CPD <- dag$CPDs[[dag$sorted[i]]]
        par <- dag$parents[[dag$sorted[i]]]
        probs <- get.CPD.probs(CPD, x[par])
        x[i] <- sample(1:length(probs), size=1, prob=probs)
    }
    return(x)
}

get.CPD.probs <- function(CPD,index.vector) {
  ## CPD: CDP array with child node at last dimension,
  ## index.vector: vector of values of parent nodes
  ## out: probability vector for child node values
  d <- dim(CPD)
  vl <- c(index.vector,list(1:d[length(d)]))
  p <- do.call("[",c(list(CPD),vl))
  c(p,1-sum(p))
}


loglik.dag <- function(dag,CPD.prior.value,data, arc.priors) {
  ## dag: a dag structure
  ## CPD.prior.value: a constant pseudocount eta for all CPDs
  ## data: rows are samples, columns are variables
  if (is.null(dag$changed))
      return(dag)
  for (i in dag$changed) {
    ## retrieve the parents
    parents <- dag$parents[[i]]
    ## create the DPD
    DPD <- table(data[, c(parents, i)])
    dag$loglik.families[i] <- loglik.family(DPD, CPD.prior.value)
  }
  prior <- prior.prob(dag, arc.priors)
  dag$loglik <- sum(dag$loglik.families) + prior
  dag$changed <- NULL
  return(dag)
}

loglik.family <- function(DPD,CPD.prior.value) {
  d <- dim(DPD)
  ## DPD.mat: rows - parent configs, cols - child levels
  DPD.mat <- matrix(DPD,ncol=d[length(d)])
  n.levels <- ncol(DPD.mat)
  n.pconfigs <- nrow(DPD.mat)
  CPD.prior.value <- CPD.prior.value / (n.levels*n.pconfigs)

  eta <- lgamma(CPD.prior.value * n.levels)
  npeta <- lgamma(apply(DPD.mat+CPD.prior.value, 1, sum))
  fp <- eta - npeta

  sp <- apply(lgamma(DPD.mat+CPD.prior.value) - lgamma(CPD.prior.value * n.levels), 1, sum)

  loglik.mat <- sum(fp + sp)

     ## use DPD.mat and CDP.prior.value and the apply function
     ## scale CDP.prior.value down by n.levels and n.pconfigs for each cell
     ## to implement the family part of equation 2
     ## and store it in loglik.mat
  return(loglik.mat)
}

max.arc.round <- 10
resample <- function(x, ...) x[sample.int(length(x), ...)]
add.arc <- function(dag) {
  ## adds a random arc to the dag keeping it acyclic
  sorted <- NULL
  round <- 0
  ## try until we get an acyclic graph or too many trials
  while (is.null(sorted) && round < max.arc.round) {
    round <- round + 1
    sorted <- NULL
    parents <- dag$parents
    ## sample a random node:
    node <- sample(dag$n,1)
    opar <- parents[[node]]
    if (length(opar) < dag$n - 1) {
      ## there is space for one more parent
      ## sample from the nonparents excluding the node itself:
      new.par <- resample((1:dag$n)[-c(opar,node)],1)
      parents[[node]] <- c(opar,new.par)
      sorted <- topsort.parents(parents)
    }
  }
  if (round == max.arc.round) {
    ## we give up, return old DAG
      dag$changed <- NULL
  } else {
    dag$parents <- parents
    dag$sorted <- sorted
    dag$changed <- node
  }
  dag
}


remove.arc <- function(dag) {
    round <- 0

    while(round < max.arc.round) {
        round <- round + 1
        parents <- dag$parents
        ##sample a random node
        node <- sample(dag$n, 1)
        opar <- parents[[node]]
        if (length(opar) > 0) {
            rpar <- sample(opar, size=1)
            opar <- opar[-which(opar == rpar)]
            parents[[node]] <- opar
            sorted <- topsort.parents(parents)
            break
        }
    }

    if (round == max.arc.round) {
        dag$changed <- NULL
    } else {
        dag$parents <- parents
        dag$sorted <- sorted
        dag$changed <- node
    }

    dag
}

revert.arc <- function(dag) {
    round <- 0
    sorted <- NULL

    while(is.null(sorted) &&
          round < max.arc.round) {
        round <- round + 1
        parents <- dag$parents

        node <- sample(dag$n, 1)
        opar <- parents[[node]]
        if (length(opar) > 0) {
            rpar <- sample(opar, size=1)
            opar <- opar[-which(opar==rpar)]
            apar <- parents[[rpar]]
            apar <- c(apar, node)
            parents[[node]] <- opar
            parents[[rpar]] <- apar
            sorted <- topsort.parents(parents)
        }

    }

    if (round == max.arc.round) {
        dag$changed <- NULL
    } else {
        dag$parents <- parents
        dag$sorted <- sorted
        dag$changed <- c(node, rpar)
    }

    dag
}

all.operations <- c("add.arc","remove.arc","revert.arc")

random.dag <- function(n,operations=all.operations,density=2) {
  ## n: number of arc moves from an empty network
  ## operations: allowed operations, default: add, revert, remove
  ## return a random DAG
  dag <- NULL
  ## define an empty DAG to start with
  dag$n <- n
  dag$parents <- rep(list(numeric(0)),n)
  dag$sorted <- 1:n
  for (i in 1:(density*n*(n-1)/2)) {
    ## choose a random operation
    op <- sample(length(operations),1)
    dag <- do.call(operations[op],list(dag))
  }
  dag$changed <- 1:n
  dag$loglik.families <- NULL
  dag$sorted <- topsort.parents(dag$parents)
  dag
}


search.dag <- function(data,n.starts=50,n.iterations=10,
                       prior.value = 0.1, arc.priors,
                       operations=all.operations) {
  ## data: matrix for variables in columns and samples in rows
  ## n.starts: number of random restarts
  ## n.iterations: number of local improvement steps allowed
  ## prior.value: the prior eta or also CPD.prior.value
  ## operations: allowed operations, default: add, revert, remove
  ## returns a DAG with the largest likelihood among those visited
  overall.dag <- NULL
  overall.dag$loglik <- -Inf
  logliks <- NULL
  for (j in 1:n.starts) { # global search
    start.dag <- random.dag(ncol(data),all.operations)
    start.dag <- loglik.dag(start.dag,prior.value,data, arc.priors)
    best.dag <- start.dag
    for (rounds in 1:n.iterations) {  # local search
      logliks <- c(logliks,best.dag$loglik)
      ## sample a random arc operation
      op <- sample(length(operations),1)
      ## apply it to the best.dag so far
      new.dag <- do.call(operations[op],list(best.dag))
      new.dag <- loglik.dag(new.dag,prior.value,data, arc.priors)
      ## if the new loglikelihood is better than best.dag
      ## keep the new dag as best.dag
      if (new.dag$loglik > best.dag$loglik) {
          best.dag <- new.dag
      }

   }
   if (best.dag$loglik > overall.dag$loglik) {
      overall.dag <- best.dag
   }
  }
  overall.dag$logliks <- logliks
  overall.dag
}


add.arc.mcmc <- function(dag) {
  n <- length(dag$parents)
  sorted <- NULL
  round <- 0
  node <- sample(n,1)
  parents <- dag$parents
  opar <- parents[[node]]
  if (length(opar) < n - 1) {
    new.par <- resample((1:n)[-c(opar,node)],1)
    parents[[node]] <- c(opar,new.par)
    sorted <- topsort.parents(parents)
  }
  if (is.null(sorted)) {
    dag$changed <- NULL
    dag$ratio <- 1
  } else {
    dag$parents <- parents
    dag$sorted <- sorted
    dag$changed <- node
    dag$ratio <- (n-length(opar)-1)/(length(opar) + 1)
  }
  dag
}

remove.arc.mcmc <- function(dag) {
    sorted <- NULL
    parents <- dag$parents
    n <- length(dag$parents)
    node <- sample(dag$n, 1)
    opar <- parents[[node]]
    k <- length(opar)
    if (length(opar) > 0) {
        rpar <- sample(opar, size=1)
        opar <- opar[-which(opar==rpar)]
        parents[[node]] <- opar
        sorted <- topsort.parents(parents)
    }

    if (is.null(sorted)) {
        dag$changed <- NULL
        dag$ratio <- 1
    } else {
        dag$parents <- parents
        dag$sorted <- sorted
        dag$changed <- node
        dag$ratio <- k / (n-k)
    }
    dag
}

revert.arc.mcmc <- function(dag) {
    sorted <- NULL
    n <- length(dag$parents)
    parents <- dag$parents

    node <- sample(dag$n, 1)
    opar <- parents[[node]]
    k <- length(opar)
    if (length(opar) > 0 ) {
        rpar <- sample(opar, size=1)
        opar <- opar[-which(opar==rpar)]
        apar <- parents[[rpar]]
        apar <- c(apar, node)
        x <- length(apar)
        parents[[node]] <- opar
        parents[[rpar]] <- apar
        sorted <- topsort.parents(parents)
    }

    if (is.null(sorted)) {
        dag$changed <- NULL
        dag$ratio <- 1
    } else {
        dag$parents <- parents
        dag$sorted <- sorted
        dag$changed <- c(node, rpar)
        dag$ratio <- k / x

    }
    dag
}

all.mcmc.operations <- c("add.arc.mcmc", "revert.arc.mcmc", "remove.arc.mcmc")

mcmc.dag <- function(data,n.iterations=2000,
                       prior.value = 0.01,arc.priors=NULL,
                       operations=all.mcmc.operations) {
  ## n.iterations: number of MCMC steps
  ## returns a list of sample dags from the MH simulation, as well
  ## as a vector of acceptance indicators
  dags <- list()
  accepts <- NULL
  old.dag <- random.dag(ncol(data),all.operations)
  old.dag <- loglik.dag(old.dag,prior.value,data)
  for (rounds in 1:n.iterations) {
      accept <- 0
      dags <- c(dags,list(old.dag))
      op <- sample(length(operations),1)
      dag <- do.call(operations[op],list(old.dag))
      if (!is.null(dag$changed)) {
          ## get the loglikelihood of the dag
          dag <- loglik.dag(dag,prior.value, data)
          ## calculate the acceptanced ratio p
          if (dag$loglik > old.dag$loglik) {
              old.dag <- dag
          } else {
              u <- runif(1)
              if (u < exp(dag$loglik-old.dag$loglik)) {
                  old.dag <- dag
              }
          }
      ## and set old.dag to the dag if it is smaller than p
      }
      accepts <- c(accepts,accept)
  }
  list(dags=dags,accepts=accepts)
}

mcmc.dag1 <- function(data,n.iterations=2000,
                       prior.value = 0.01,arc.priors=NULL,
                       operations=all.mcmc.operations) {
  ## n.iterations: number of MCMC steps
  ## returns a list of sample dags from the MH simulation, as well
  ## as a vector of acceptance indicators
  dags <- list()
  accepts <- NULL
  old.dag <- random.dag(ncol(data),all.operations)
  old.dag <- loglik.dag(old.dag,prior.value,data, arc.priors)
  for (rounds in 1:n.iterations) {
      accept <- 0
      dags <- c(dags,list(old.dag))
      op <- sample(length(operations),1)
      dag <- do.call(operations[op],list(old.dag))
      if (!is.null(dag$changed)) {
          ## get the loglikelihood of the dag
          dag <- loglik.dag(dag,prior.value, data, arc.priors)
          accept <- min(1, exp(dag$loglik-old.dag$loglik+log(dag$ratio)))
          u <- runif(1)
          if (u < accept) {
              old.dag <- dag
          }
      }
      accepts <- c(accepts,accept)
  }
  list(dags=dags,accepts=accepts)
}

make.arc.priors <- function(n,arc.lst,prob.lst) {
  arcs <-  no.arcs <- matrix(1/3,nrow=n,ncol=n)
  diag(arcs) <- diag(no.arcs) <- 1
  for (i in 1:length(arc.lst)) {
    arc <- arc.lst[[i]]
    probs <- prob.lst[[i]]
    arcs[arc[1],arc[2]] <- probs[1]
    arcs[arc[2],arc[1]] <- probs[2]
    no.arcs[arc[1],arc[2]] <- probs[3]
    no.arcs[arc[2],arc[1]] <- probs[3]
  }
  list(arcs=log(arcs),no.arcs=log(no.arcs),no.arcs.total=sum(log(no.arcs)))
}

prior.prob <- function(dag, arc.priors) {
    p <- arc.priors$no.arcs.total
    for (node in 1:dag$n) {
        parents <- dag$parents[[node]]
        p <- p + sum(vapply(parents, function(p) { arc.priors$arcs[p, node] },
                            FUN.VALUE = c(1)))
        p <- p - sum(vapply(parents, function(p) { arc.priors$no.arcs[p, node] },
                            FUN.VALUE = c(1)))
    }

    return(p)
}




