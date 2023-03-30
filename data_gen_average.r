library(MASS)
library(mediation)
library(lavaan)

### Extract variables from JOBS II data ###
### Read data from mediation package and export as a dat file, recode some variables ###
data("jobs", package = "mediation")
path <- "C:/Users/wqemi/OneDrive - Florida State University/MyFSU_OneDrive/Documents/FSU/moderated mediation/sensitivity analysis/JOBS II data/jobsII.dat"
write.table(jobs, path, sep = "\t", quote = FALSE)

### Read data from csv file ###
jobsII <- read.csv("C:/Users/wqemi/OneDrive - Florida State University/MyFSU_OneDrive/Documents/FSU/moderated mediation/sensitivity analysis/JOBS II data/JOBSII.csv")

### Gender ###
gender <- jobsII$sex
### Age ###
age <- jobsII$age
### Occupations ###
occupation <- jobsII$occp
### Dummy code of occupation ### 
occupation.1 <- ifelse(occupation == 2, 1, 0)
occupation.2 <- ifelse(occupation == 3, 1, 0)
occupation.3 <- ifelse(occupation == 4, 1, 0)
occupation.4 <- ifelse(occupation == 5, 1, 0)
occupation.5 <- ifelse(occupation == 6, 1, 0)
occupation.6 <- ifelse(occupation == 7, 1, 0)
### Level of economic hardship pre-treatment ###
economic.hardship <- jobsII$econ_hard
### Marital status ###
marital.status <- jobsII$marital
### Dummy code of marital.status ### 
marital.status.1 <- ifelse(marital.status == 2, 1, 0)
marital.status.2 <- ifelse(marital.status == 3, 1, 0)
marital.status.3 <- ifelse(marital.status == 4, 1, 0)
marital.status.4 <- ifelse(marital.status == 5, 1, 0)
### White or not white ###
white <- jobsII$nonwhite
white <- white - 1
### Eduction level ###
education <- jobsII$educ
### Dummy code of education ### 
education.1 <- ifelse(education == 2, 1, 0)
education.2 <- ifelse(education == 3, 1, 0)
education.3 <- ifelse(education == 4, 1, 0)
education.4 <- ifelse(education == 5, 1, 0)
### Condition ###
condition <- jobsII$treat
### Depression ###
depression <- jobsII$depress2
### Self-efficacy ###
self.efficacy <- jobsII$job_seek
### Sample size ###
N <- nrow(jobsII)
### Create a new dataset ###
jobs.data.new <- cbind(condition, gender, age, occupation, occupation.1, occupation.2, occupation.3, occupation.4, occupation.5, occupation.6, economic.hardship, marital.status, marital.status.1, marital.status.2, 
                 marital.status.3, marital.status.4, white, education, education.1, education.2, education.3, education.4, self.efficacy, depression)

### Factor loading ###
factor.loading <- 0.6
### number of items ###
item.number.M <- 6
item.number.Y <- 6
### cut point for M ###
cut.point.M <- c(.05, .15, .3, .5)
### cut point for Y ###
cut.point.Y <- c(.5, .7, .85, .95)
### Mediation analysis to extract coefficients ###
model <- '
### Regression model ###
self.efficacy ~ a * condition + gender + age + occupation.1 + occupation.2 + occupation.3 + occupation.4 + occupation.5 + occupation.6 + economic.hardship + white + marital.status.1 + marital.status.2 + marital.status.3 +
            marital.status.4 + education.1 + education.2 + education.3 + education.4
depression ~ b * self.efficacy + c * condition + gender + age + occupation.1 + occupation.2 + occupation.3 + occupation.4 + occupation.5 + occupation.6 + economic.hardship + white + marital.status.1 + marital.status.2 + 
            marital.status.3 + marital.status.4 + education.1 + education.2 + education.3 + education.4
'
fit <- sem(model, data = jobs.data.new)
coefficient.a <- parameterestimates(fit)[1, 5]
coefficient.b <- parameterestimates(fit)[20, 5]
coefficient.c <- parameterestimates(fit)[21, 5]
coefficient.gender.M <- parameterestimates(fit)[2, 5]
coefficient.age.M <- parameterestimates(fit)[3, 5]
coefficient.occupation.1.M <- parameterestimates(fit)[4, 5]
coefficient.occupation.2.M <- parameterestimates(fit)[5, 5]
coefficient.occupation.3.M <- parameterestimates(fit)[6, 5]
coefficient.occupation.4.M <- parameterestimates(fit)[7, 5]
coefficient.occupation.5.M <- parameterestimates(fit)[8, 5]
coefficient.occupation.6.M <- parameterestimates(fit)[9, 5]
coefficient.economic.hardship.M <- parameterestimates(fit)[10, 5]
coefficient.white.M <- parameterestimates(fit)[11, 5]
coefficient.marital.status.1.M <- parameterestimates(fit)[12, 5]
coefficient.marital.status.2.M <- parameterestimates(fit)[13, 5]
coefficient.marital.status.3.M <- parameterestimates(fit)[14, 5]
coefficient.marital.status.4.M <- parameterestimates(fit)[15, 5]
coefficient.education.1.M <- parameterestimates(fit)[16, 5]
coefficient.education.2.M <- parameterestimates(fit)[17, 5]
coefficient.education.3.M <- parameterestimates(fit)[18, 5]
coefficient.education.4.M <- parameterestimates(fit)[19, 5]
coefficient.gender.Y <- parameterestimates(fit)[22, 5]
coefficient.age.Y <- parameterestimates(fit)[23, 5]
coefficient.occupation.1.Y <- parameterestimates(fit)[24, 5]
coefficient.occupation.2.Y <- parameterestimates(fit)[25, 5]
coefficient.occupation.3.Y <- parameterestimates(fit)[26, 5]
coefficient.occupation.4.Y <- parameterestimates(fit)[27, 5]
coefficient.occupation.5.Y <- parameterestimates(fit)[28, 5]
coefficient.occupation.6.Y <- parameterestimates(fit)[29, 5]
coefficient.economic.hardship.Y <- parameterestimates(fit)[30, 5]
coefficient.white.Y <- parameterestimates(fit)[31, 5]
coefficient.marital.status.1.Y <- parameterestimates(fit)[32, 5]
coefficient.marital.status.2.Y <- parameterestimates(fit)[33, 5]
coefficient.marital.status.3.Y <- parameterestimates(fit)[34, 5]
coefficient.marital.status.4.Y <- parameterestimates(fit)[35, 5]
coefficient.education.1.Y <- parameterestimates(fit)[36, 5]
coefficient.education.2.Y <- parameterestimates(fit)[37, 5]
coefficient.education.3.Y <- parameterestimates(fit)[38, 5]
coefficient.education.4.Y <- parameterestimates(fit)[39, 5]

obs.M.total <- list()
obs.Y.total <- list()

for (k in 1:100) {
    
    ### Correlation between M and Y ###
    rho <- runif(1, -1, 1)
  
    ### Residual variance-covariance matrix of M and Y ###
    mu <- numeric(2)
    sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    ### Generate residuals for M and Y ###
    residual.MY <- mvrnorm(N, mu, sigma)
    ### Calculate latent variable M ###
    X <- condition
    M = coefficient.a * X + coefficient.gender.M * gender + coefficient.age.M * age + coefficient.occupation.1.M * occupation.1 + coefficient.occupation.2.M * occupation.2 + 
        coefficient.occupation.3.M * occupation.3 + coefficient.occupation.4.M * occupation.4 + coefficient.occupation.5.M * occupation.5 + coefficient.occupation.6.M * occupation.6 + 
        coefficient.economic.hardship.M * economic.hardship + coefficient.white.M * white + coefficient.marital.status.1.M * marital.status.1 + coefficient.marital.status.2.M * marital.status.2 +
        coefficient.marital.status.3.M * marital.status.3 + coefficient.marital.status.4.M * marital.status.4 + coefficient.education.1.M * education.1 +  coefficient.education.2.M * education.2 + 
        coefficient.education.3.M * education.3 + coefficient.education.4.M * education.4 + residual.MY[, 1]
    ### Calcualte latent variable Y ###
    Y = coefficient.b * M + coefficient.c * X + coefficient.gender.Y * gender + coefficient.age.Y * age + coefficient.occupation.1.Y * occupation.1 + coefficient.occupation.2.Y * occupation.2 +
        coefficient.occupation.3.Y * occupation.3 + coefficient.occupation.4.Y * occupation.4 + coefficient.occupation.5.Y * occupation.5 + coefficient.occupation.6.Y * occupation.6 + 
        coefficient.economic.hardship.Y * economic.hardship + coefficient.white.Y * white + coefficient.marital.status.1.Y * marital.status.1 + coefficient.marital.status.2.Y * marital.status.2 +
        coefficient.marital.status.3.Y * marital.status.3 + coefficient.marital.status.4.Y * marital.status.4 + coefficient.education.1.Y * education.1 + coefficient.education.2.Y * education.2 + 
        coefficient.education.3.Y * education.3 + coefficient.education.4.Y * education.4 + residual.MY[, 2]

    ### Job-search self-efficacy using a six-item index; Reliability: 0.87 ###
    ### Variance for residuals: 0.3227587 ###
    ### Mean matrix of meausrement error for M ###
    mu.measuremenet.error.M <- rep(0, item.number.M)
    ### Variance-covariance matrix of measurement error for M ###
    var.measuremenet.error.M <- diag(0.3227587, item.number.M)
    ### Generate residual data for observed variables of M
    obs.residual.M <- mvrnorm(N, mu.measuremenet.error.M, var.measuremenet.error.M)
    ### Factor loading matrix ###
    factor.loading.M <- rep(factor.loading, item.number.M)
    ### Generate continous observed data for M ###
    obs.M.total[[k]] <- M %*% t(factor.loading.M) + obs.residual.M
 
    ### Depressive symptoms level was measured with a subscale of 6 items; Reliability: 0.9 ###
    ### Variance for residuals: 0.24 ###
    ### Mean matrix of meausrement error for M ###
    mu.measuremenet.error.Y <- rep(0, item.number.Y)
    ### Variance-covariance matrix of measurement error for M ###
    var.measuremenet.error.Y <- diag(0.24, item.number.Y)
    ### Generate residual data for observed variables of M
    obs.residual.Y <- mvrnorm(N, mu.measuremenet.error.Y, var.measuremenet.error.Y)
    ### Factor loading matrix ###
    factor.loading.Y <- rep(factor.loading, item.number.Y)
    ### Generate continous observed data for M ###
    obs.Y.total[[k]] <- Y %*% t(factor.loading.Y) + obs.residual.Y
}

obs.M <- obs.Y <- matrix(NA, N, ncol(obs.M.total[[1]]))
### Average obs.M and obs.Y ###
for (i in 1:ncol(obs.M.total[[1]])) {
    result.M <- lapply(obs.M.total, "[", , i)
    result.Y <- lapply(obs.Y.total, "[", , i)
    obs.M[, i] <- rowMeans(do.call(cbind, result.M))
    obs.Y[, i] <- rowMeans(do.call(cbind, result.Y))
}

### Generate ordinal data for M ### 
ordinal.obs.M <- matrix(NA, N, item.number.M)
for (i in 1 : item.number.M) {
    cutoff1.M <- quantile(obs.M[, i], cut.point.M[1])
    cutoff2.M <- quantile(obs.M[, i], cut.point.M[2])
    cutoff3.M <- quantile(obs.M[, i], cut.point.M[3])
    cutoff4.M <- quantile(obs.M[, i], cut.point.M[4])
    ordinal.obs.M[, i] <- cut(obs.M[, i], breaks = c(min(obs.M[, i]) - 1, cutoff1.M, cutoff2.M, cutoff3.M, cutoff4.M, max(obs.M[, i]) + 1), labels = c(1, 2, 3, 4, 5))
}

### Generate ordinal data for Y ### 
ordinal.obs.Y <- matrix(NA, N, item.number.Y)
for (i in 1 : item.number.Y) {
    cutoff1.Y <- quantile(obs.Y[, i], cut.point.Y[1])
    cutoff2.Y <- quantile(obs.Y[, i], cut.point.Y[2])
    cutoff3.Y <- quantile(obs.Y[, i], cut.point.Y[3])
    cutoff4.Y <- quantile(obs.Y[, i], cut.point.Y[4])
    ordinal.obs.Y[, i] <- cut(obs.Y[, i], breaks = c(min(obs.Y[, i]) - 1, cutoff1.Y, cutoff2.Y, cutoff3.Y, cutoff4.Y, max(obs.Y[, i]) + 1), labels = c(1, 2, 3, 4, 5))
}

### Dataset with dummy variable
jobs.data.ordinal.analysis <- cbind(condition, gender, age, occupation.1, occupation.2, occupation.3, occupation.4, occupation.5, occupation.6, economic.hardship, marital.status.1, marital.status.2, marital.status.3, marital.status.4, white, education.1, education.2, education.3, education.4, ordinal.obs.M, ordinal.obs.Y)
filename <- "C:/Users/wqemi/OneDrive - Florida State University/MyFSU_OneDrive/Documents/FSU/moderated mediation/sensitivity analysis/JOBS II data/job_analysis_average.dat"
write.table(jobs.data.ordinal.analysis, filename, col.names = FALSE, row.names = FALSE, quote = FALSE)

filename <- "C:/Users/wqemi/OneDrive - Florida State University/MyFSU_OneDrive/Documents/FSU/moderated mediation/sensitivity analysis/JOBS II data/job_analysis_average.csv"
write.table(jobs.data.ordinal.analysis, filename, col.names = TRUE, row.names = FALSE, quote = FALSE)