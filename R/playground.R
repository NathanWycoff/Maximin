#!/usr/bin/Rscript
#  R/playground.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.21.2019
library(nloptr)
library(maximin)

set.seed(123)
Ns <- seq(5,15,by=10)
P <- 2
iters <- 30

time_res <- as.data.frame(matrix(NA, nrow = length(Ns)*P*2, ncol = 3))
obj_res <- as.data.frame(matrix(NA, nrow = length(Ns)*P*2, ncol = 3))
colnames(time_res) <- c("time", "method", "N")
colnames(obj_res) <- c("obj", "method", "N")

it <- 0
for (N_i in 1:length(Ns)) {
    cat("N:", Ns[N_i], "\n")
    for (iter in 1:iters) {
        cat("iter:", iter, "\n")
        it <- it + 1
        N <- Ns[N_i]

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

        fb_time <- system.time(fb_obj <- max(maximin(N, P, 10*N)$mi))[3]

        my_time <- system.time(my_obj <- -slsqp(x0 = c(runif(N*P), 0), fn = objective, gr = gradient, 
              hin = hin, hinjac = hinjac,
              lower = rep(0, N*P+1), upper = rep(1, N*P+1))$value)[3]

        #cat("Furong/Bobby Time:", fb_time, "\n")
        #cat("My time:", my_time, "\n")
        #cat("Furong/Bobby Time:", fb_obj, "\n")
        #cat("Furong/Bobby Time:", my_obj, "\n")

        time_res[2*(it-1)+1,] <- c(my_time, 'Me', N)
        time_res[2*(it-1)+2,] <- c(fb_time, 'FB', N)
        obj_res[2*(it-1)+1,] <- c(my_obj, 'Me', N)
        obj_res[2*(it-1)+2,] <- c(fb_obj, 'FB', N)
    }
}

# Fix type issues
time_res$time <- as.numeric(time_res$time)
time_res$N <- as.numeric(time_res$N)
obj_res$obj <- as.numeric(obj_res$obj)
obj_res$N <- as.numeric(obj_res$N)

#save(time_res, obj_res, file = 'data/cool_sim_results.RData')
