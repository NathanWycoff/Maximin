library(nloptr)
library(Rcpp)

myminphi <- function(n, m, p=2, X=NULL, T=10000, save.prog=FALSE) 
{  
  ## for timing purposes
  tic <- proc.time()[3]
  
  ## initial design and sanity check
  if(is.null(X)) X <- matrix(runif(n*m), ncol=m)
  else if(nrow(X) != n || ncol(X) != m) stop("n and m don't match X")
  
  ## progress meter
  if(save.prog) {
    prog <- data.frame(phi=rep(NA, T), md=rep(NA, T))
  } else prog <- NULL  
  
  ## random swapping
  for(t in 1:T) {
    
    ## propose swapped row
    row.out <- sample(1:n, 1)  ## select any row
    
    ## calculate distances to that location
    d <- drop(distance(X, X[row.out,,drop=FALSE]))
    phi <- sum(d[-row.out]^(-p))^(1/p)
    
    ## perform the swap with a random new row
    xold <- X[row.out,] 
    X[row.out,] <- runif(m)
    
    ## updated distance and phi for new row
    dnew <- drop(distance(X, X[row.out,,drop=FALSE]))
    phiprime <- sum(dnew[-row.out]^(-p))^(1/p)
    
    ## accept or reject; nothing to do if accept
    if(phiprime >= phi) { ## corresponds to reject
      X[row.out,] <- xold 
    }   
    
    ## saving progress requires full distance matrix
    if(save.prog) {
      if(phiprime < phi || t == 1) { ## accept above
        d <- distance(X)
        dvec <- d[upper.tri(d)]
        phi <- sum(dvec^(-p))^(1/p)
        prog$phi[t] <- phi
        prog$md[t] <- min(dvec)
      } else {
        prog$phi[t] <- prog$phi[t-1]
        prog$md[t] <- prog$md[t-1]
      }
    }
  }
  
  ## end timing
  toc <- proc.time()[3]
  
  return(list(X=X, prog=prog, time=as.numeric(toc-tic)))
}

mymaximin.rf <- function(n, m, X=NULL, T=10000) 
{  
  ## for timing purposes
  tic <- proc.time()[3]
  
  ## initial design and sanity check
  if(is.null(X)) X <- matrix(runif(n*m), ncol=m)
  else if(nrow(X) != n || ncol(X) != m) stop("n and m don't match X")
  
  ## initial distance calculations
  d <- distance(X)
  dvec <- d[upper.tri(d)]
  wmd <- which.min(dvec)
  md <- dvec[wmd]
  
  ## calculate pair (row, col) involved in minimum dist
  rs <- 0:(n-1)
  crs <- cumsum(rs)
  col <- which.min(wmd/crs > 1)
  row <- wmd - crs[col-1]
  
  ## progress meter
  prog <- rep(NA, T)
  
  ## random swapping
  for(t in 1:T) {
    
    ## propose swapped row
    row.out <- sample(c(row,col), 1)  ## select only from pair
    xold <- X[row.out,]               ## save old rows
    X[row.out,] <- runif(m)           ## random new row
    
    ## @@ modifications to compute only new row/col of d
    dnew <- drop(distance(X, X[row.out,,drop=FALSE]))
    mdprime <- min(dnew[-row.out])
    
    ## accept or reject
    if(mdprime > md) {                ## accept
      
      ## @@ update distance matrix (only new row) and row/col
      d[row.out,] <- d[,row.out] <- dnew
      dvec <- d[upper.tri(d)]
      wmd <- which.min(dvec)
      md <- dvec[wmd]                 
      col <- which.min(wmd/crs > 1)
      row <- wmd - crs[col-1]
      
    } else {                          ## reject
      X[row.out,] <- xold 
    }   
    
    ## update progress
    prog[t] <- md
  }
  
  ## end timing
  toc <- proc.time()[3]
  
  return(list(X=X, prog=prog, time=as.numeric(toc-tic)))
}

#' Math Programming based Maximin
#'
#' This script formulates the Maximin problem as the following math program:
#'
#' max t
#'
#' S.T.
#'
#' t <= dist(xi,xj) for all 1 <= i <= j <= N
#'
#' With xi the i'th design point.
#'
#' @param N The desired number of design points.
#' @param P The dimension of the design space.
#' @param design_start The initialization for optimization routine; random if not provided.
opt_maximin <- function(N, P, design_start) {
    if (missing(N)) {
        N <- nrow(design_start)
    }
    if (missing(P)) {
        P <- ncol(design_start)
    }

    # Define distance matrices
    D <- list()
    for (n1 in 1:(N-1)) {
        D[[n1]] <- list()
        for (n2 in (n1+1):N) {
            Dij <- matrix(0, nrow = N*P, ncol = N*P)

            start1 <- (n1-1)*P + 1
            start2 <- (n2-1)*P + 1
            Dij[(start1):(start1+P-1),(start1):(start1+P-1)] = diag(P);
            Dij[(start2):(start2+P-1),(start2):(start2+P-1)] = diag(P);
            Dij[(start1):(start1+P-1),(start2):(start2+P-1)] = -diag(P);
            Dij[(start2):(start2+P-1),(start1):(start1+P-1)] = -diag(P);

            D[[n1]][[length(D[[n1]])+1]] <- Dij
        }
    }

    # Define objective and constraints
    n_constr <- choose(N,2)
    objective <- function(par) -par[N*P+1]
    gradient <- function(par) c(rep(0, N*P), -1)

    hin <- function(par) {
        x <- par[1:(N*P)]
        t <- par[(N*P)+1]

        ret <- rep(NA, n_constr)
        i <- 1
        for (n1 in 1:(N-1)) {
            for (n2 in (n1+1):N) {
                ret[i] <- as.numeric(t(x) %*% D[[n1]][[n2-n1]] %*% x) - t
                i <- i + 1
            }
        }

        return(ret)
    }

    hinjac <- function(par) {
        x <- par[1:(N*P)]
        t <- par[(N*P)+1]

        RET <- matrix(NA, nrow = n_constr, ncol = N*P+1)
        i <- 1
        for (n1 in 1:(N-1)) {
            for (n2 in (n1+1):N) {
                RET[i,] <- c(2 * D[[n1]][[n2-n1]] %*% x, -1)
                i <- i + 1
            }
        }

        return(RET)
    }

    if (missing(design_start)) {
        ides <- c(runif(N*P), 0)
    } else {
        if (nrow(design_start) != N || ncol(design_start != P)) {
            stop("Starting design does not agree with provided N, P")
        }
        ides <- as.numeric(design_start)
    }
    my_time <- system.time(my_obj <- slsqp(x0 = ides, fn = objective, gr = gradient, 
                                           hin = hin, hinjac = hinjac, lower = rep(0, N*P+1), upper = rep(1, N*P+1)))[3]

    return(list(opt = my_obj, time=my_time))
}

cppFunction('NumericMatrix c_hinjac(NumericVector par, int N, int P, int n_constr) {

            NumericMatrix RET(n_constr, N*P+1);

            NumericVector dist2(P);

            // Initialize t derivatives:
            for (int i = 0; i < n_constr; i++) {
                RET(i,N*P) = -1;
            }

            int i = 0;
            for (int n1 = 0; n1 < N-1; n1++) {
                for (int n2 = n1+1; n2 < N; n2++) {
                    for (int p = 0; p < P; p++) {
                        dist2(p) = 2*(par(n1*P+p) - par(n2*P+p));
                        RET(i,n1*P+p) =  dist2(p);
                        RET(i,n2*P+p) =  -dist2(p);
                    }
                    i++;
                }
            }

            return RET;
}')

#' Math Programming based Maximin
#'
#' This script formulates the Maximin problem as the following math program:
#'
#' max t
#'
#' S.T.
#'
#' t <= dist(xi,xj) for all 1 <= i <= j <= N
#'
#' With xi the i'th design point.
#'
#' @param N The desired number of design points.
#' @param P The dimension of the design space.
#' @param design_start The initialization for optimization routine; random if not provided.
#' @return A list giving the optimal design and the time elapsed in creating the design.
opt_maximin_fast <- function(N, P, design_start) {
    if (missing(N)) {
        N <- nrow(design_start)
    }
    if (missing(P)) {
        P <- ncol(design_start)
    }

    # Define objective and constraints
    n_constr <- choose(N,2)
    objective <- function(par) -par[N*P+1]
    gradient <- function(par) c(rep(0, N*P), -1)

    hin <- function(par) {
        X <- matrix(par[1:(N*P)], ncol = P, byrow=T)
        t <- par[(N*P)+1]
        return(as.numeric(dist(X))^2 - t)
    }

    #hinjac <- function(par) {
    #    X <- matrix(par[1:(N*P)], ncol = P, byrow=T)
    #    t <- par[(N*P)+1]

    #    RET <- matrix(0, nrow = n_constr, ncol = N*P+1)
    #    RET[,N*P+1] <- -1
    #    i <- 1
    #    for (n1 in 1:(N-1)) {
    #        for (n2 in (n1+1):N) {
    #            diff <- 2*(X[n1,] - X[n2,])
    #            #RET[i,] <- c(rep(0, P*(n1-1)), diff, rep(0, P*(n2-n1-1)), -diff, rep(0, P*(N-n2)), -1)
    #            RET[i,(P*(n1-1)+1):(P*(n1))] <- diff
    #            RET[i,(P*(n2-1)+1):(P*(n2))] <- -diff
    #            i <- i + 1
    #        }
    #    }

    #    return(RET)
    #}

    if (missing(design_start)) {
        ides <- c(runif(N*P), 0)
    } else {
        if (nrow(design_start) != N || ncol(design_start != P)) {
            stop("Starting design does not agree with provided N, P")
        }
        ides <- as.numeric(design_start)
    }

    hinc_wrap <- function(par) c_hinjac(par, N, P, choose(50,2))

    my_time <- system.time(my_obj <- slsqp(x0 = ides, fn = objective, gr = gradient, 
                                           hin = hin, hinjac = hinc_wrap, lower = rep(0, N*P+1), upper = rep(1, N*P+1)))[3]

    return(list(opt = my_obj, time=my_time))
}


