#!/usr/bin/Rscript
#  R/speed_lab.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 05.28.2019
source('R/Benchmark.R')

N <- 50
P <- 2

set.seed(123)
opt_maximin(50,2)
set.seed(123)
opt_maximin_fast(50,2)

# The inequality evals and gradients are 100 and 20 times faster respectively, why don't we see a commensurate increase in the overall optimization?
