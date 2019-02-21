#!/usr/bin/Rscript
#  R/analyze_playground.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.21.2019

load("data/cool_sim_results.RData")

pdf('images/playground.pdf', width = 12)
par(mfrow=c(1,2))
myplot=boxplot(obj ~ method*N , data=obj_res, ylab = "Minimum Distance", 
               main = "Effectiveness")
myplot=boxplot(time ~ method*N , data=time_res, ylab = "Execution Time (Seconds)", 
               main = "Efficiency")
dev.off()
