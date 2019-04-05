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