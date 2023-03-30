### Ananlysis bifactor model ###
library(data.table)
library(lavaan)

setwd("C:/Users/wqemi/OneDrive - Florida State University/MyFSU_OneDrive/Documents/FSU/moderated mediation/sensitivity analysis/JOBS II data/")
job <- fread("job_analysis_average.dat")
names(job) <- c("X", paste0("CMY", 1:18), paste0("M", 1:6), paste0("Y", 1:6))

model <- "
M =~ NA * M1 + M2 + M3 + M4 + M5 + M6
Y =~ NA * Y1 + Y2 + Y3 + Y4 + Y5 + Y6

M ~ a * X + CMY1 + CMY2 + CMY3 + CMY4 + CMY5 + CMY6 + CMY7 + CMY8 + CMY9 + CMY10 + CMY11 + CMY12
    + CMY13 + CMY14 + CMY15 + CMY16 + CMY17 + CMY18
Y ~ b * M + c * X + CMY1 + CMY2 + CMY3 + CMY4 + CMY5 + CMY6 + CMY7 + CMY8 + CMY9 + CMY10 + CMY11 + CMY12
    + CMY13 + CMY14 + CMY15 + CMY16 + CMY17 + CMY18

M ~~ 1 * M
Y ~~ 1 * Y

ind:= a * b
"

fit <- sem(model, data = job, se = "bootstrap", bootstrap = 1000)
tail(parameterEstimates(fit, boot.ci.type = "perc", level = 0.95, ci = TRUE))