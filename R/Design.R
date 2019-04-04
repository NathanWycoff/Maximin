library(nloptr)

N = 300
P = 2
initialize = function(N, P = 2){
  #' we only need to do squares that could be smaller than N
  sN = N^(1/P)
  #' the two designs I analytically know the solution for are square numbers
  #' side note: not necessarily squares, but I wrote this before I remembered designs
  #' existed in more than 2 dimensions
  squares = (1:floor(sN))^P
  #' and sums of consecutive square numbers
  sums = squares[-length(squares)] + squares[-1]
  designs = sort(c(squares, sums), decreasing = T)
  #' find the biggest number that is less than the desired size
  d = designs[which(designs<N)[1]]
  
  if(d %in% squares[-1]){
    #' if it is a square number, make the grid
    s = d^(1/P)
    # init = (expand.grid(1:s,1:s)-1)/(s-1)
    init = matrix(1:s, ncol = 1)
    if(P>1){
      for(i in 2:P){
        t1 = cbind(init, 1)
        for(j in 2:s){
          t2 = cbind(init, j)
          t1 = rbind(t1, t2)
        }
        init = t1
      }
    }
    init = (init-1)/(s-1)
  }else{
    #'the extra check makes sure N isn't < 4 in which case the grids will break because we
    #'divide by 0
    if(d %in% sums){
      #' if it isn't a square number, make the interlocking grids
      s = which(sums == d)
      s2 = s+1
      i1 = matrix(1:s, ncol = 1)
      if(P>1){
        for(i in 2:P){
          t1 = cbind(i1, 1)
          for(j in 2:s){
            t2 = cbind(i1, j)
            t1 = rbind(t1, t2)
          }
          i1 = t1
        }
      }
      i1 = (i1-0.5)/(s)
      i2 = matrix(1:s2, ncol = 1)
      if(P>1){
        for(i in 2:P){
          t1 = cbind(i2, 1)
          for(j in 2:s2){
            t2 = cbind(i2, j)
            t1 = rbind(t1, t2)
          }
          i2 = t1
        }
      }
      i2 = (i2-1)/(s2-1)
      init = rbind(i1, i2)
    }
  }
  #' generate extra points
  init = rbind(init, matrix(runif(P*(N-d)), ncol = P))
  return(init %>% t %>% c)
}


maximin_opt = function(N, P, iters){
  time = md = 0
  for (iter in 1:iters) {
    cat("iter:", iter, "\n")
    
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
    
    ##Optimization Problem:
    ##
    ## max t
    ##
    ## S.T.
    ##
    ## t <= t(x) %*% A[[i]] %*% x
    ##
    ## A is a list of constraints.
    ##
    ## x = vec(X), and t is put at the end for optim purposes.
    
    # Define objective and constraints
    n_constr <- sum(sapply(D, length))
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
    
    #fb_time <- system.time(fb_obj <- max(maximin(N, P, 10*N)$mi))[3]
    lower = c(rep(0, N*P+1))
    upper = c(rep(1, N*P+1))
    start = c(initialize(N, P), 0)
    
    my_time <- system.time(my_obj <- slsqp(x0 = start, fn = objective, gr = gradient,
                                           hin = hin, hinjac = hinjac,
                                           lower = lower, upper = upper)$par)[3]
    time = time + my_time
    mdprime = my_obj[length(my_obj)]
    if(mdprime>md){
      md = mdprime
      d = my_obj[-length(my_obj)] %>% matrix(ncol = P, byrow = T)
    }
  }
  return(list(design = d, time = time, min = md))
}

check = maximin_opt(N, P, iters = 1)
