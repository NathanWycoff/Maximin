library(microbenchmark)

N <- 10
P <- 2
par <- runif(N*P+1)

microbenchmark(hin(par))
microbenchmark(hin_fast(par))

A <- matrix(0,2,2)

microbenchmark(hinjac(par))
microbenchmark(hinjac_fast(par))

hinjac <- function(par) {
    X <- matrix(par[1:(N*P)], ncol = P, byrow=T)
    t <- par[(N*P)+1]

    RET <- matrix(0, nrow = n_constr, ncol = N*P+1)
    RET[,N*P+1] <- -1
    i <- 1
    for (n1 in 1:(N-1)) {
        for (n2 in (n1+1):N) {
            diff <- 2*(X[n1,] - X[n2,])
            #RET[i,] <- c(rep(0, P*(n1-1)), diff, rep(0, P*(n2-n1-1)), -diff, rep(0, P*(N-n2)), -1)
            RET[i,(P*(n1-1)+1):(P*(n1))] <- diff
            RET[i,(P*(n2-1)+1):(P*(n2))] <- -diff
            i <- i + 1
        }
    }

    return(RET)
}

library(nloptr)
library(Rcpp)


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
N <- 10
P <- 2
par <- runif(N*P+1)
n_constr <- choose(N,2)

microbenchmark(c_hinjac(par, N, P, n_constr))
microbenchmark(hinjac(par))
