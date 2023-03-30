### Ananlysis bifactor model ###
library(data.table)
library(lavaan)

setwd("C:/Users/wqemi/OneDrive - Florida State University/MyFSU_OneDrive/Documents/FSU/moderated mediation/sensitivity analysis/apology/")
ap <- fread("ap_analysis_average.dat")
names(ap) <- c("X", paste0("M", 1:6), paste0("YG", 1:7), paste0("YS1", 1:4), paste0("YS2", 1:6))

model <- "
M =~ 1 * M1 + M2 + M3 + M4 + M5 + M6
Y =~ 1 * YG1 + YG2 + YG3 + YG4 + YG5 + YG6 + YG7 + YS11 + YS12 + YS13 + YS14 + YS21 + YS22 + YS23 + YS24 + YS25 + YS26
YS1 =~ 1 * YS11 + YS12 + YS13 + YS14
YS2 =~ 1 * YS21 + YS22 + YS23 + YS24 + YS25 + YS26

M ~ a * X
Y ~ b * M

Y ~~ 0 * YS1
Y ~~ 0 * YS2
YS1 ~~ 0 * YS2

ind:= a * b
"

fit <- sem(model, data = ap, se = "bootstrap", bootstrap = 1000)
summary(fit)
