####################### Generate data for bifactor example ##############################
library(magrittr)
library(MASS)
####### Choose the condition that transgressor value is low #####
n.apology <- 234 # low transgressor value with apology
n.noapology <- 240 # low transgressor value with no apology
N <- n.apology + n.noapology

### Factor loading ###
factor.loading <- .6
factor.loading.g <- .8
factor.loading.s <- .45
### Number of items ###
item.number.rv <- 6 # N of relationship value scale
item.number.rev <- 4 # N of subscale revenge of forgiveness
item.number.ben <- 6 # N of subscale benevolence of forgiveness
item.number.avo <- 7 # N of subscale avoidance of forgiveness
item.number.forgive <- sum(item.number.rev, item.number.ben, item.number.avo) # N of forgiveness
### Coefficients from  Forster et. al. (2021) figure 4. under the condition transgressor value: low ###
coef.a <- .502 # Coefficient from apology to relationship value
coef.b <- .777 # Coefficient from relationship value to forgiveness

rel.rv <- .95 # Reliability of relationship value: .95 
rel.forgive <- .901 # Reliability of forgiveness: .901 

### Cut point for relationship value ###
cut.point.M <- c(.01, .05, .17, .35, .65, .83, .95, .99)

### Cut point for forgiveness ###
cut.point.Y <- c(.1, .3, .7, .9)

obs.M.total <- list()
obs.Y.total <- list()

### Generate X ###
X <- rbinom(N, 1, (240 / N))


if (sum(X == 1) > n.apology) {

# If n.apology > 234, then randomly select a number and make it be 0.
    repeat {

        position <- sample(c(1:N), 1, replace = FALSE)

        if (X[position] == 1) {
            X[position] <- 0
        }

        if (sum(X == 1) == n.apology) {
            break
        }
        
    }

} else if (sum(X == 0) > n.noapology) {

# If n.noapology > 240, then randomly select a number and make it be 1.
    repeat {

        position <- sample(c(1:N), 1, replace = FALSE)

        if (X[position] == 0) {
            X[position] <- 1
        }

        if (sum(X == 0) == n.noapology) {
            break
        }
        
    }

}


### Write function to calculate residual variances for each latent model ###
### One factor model ###
omega.res <- function(rel, loading, nitem) {
    sum.res <- ((1 - rel) * (sum(loading))) ^ 2 / rel
    res <- sum.res / nitem
    return(res)
}
rv.loading <- c(factor.loading, item.number.rv)

rv.res <- omega.res(rel.rv, rv.loading, item.number.rv)

### Bifactor factor model ###
h.omega.res <- function(rel, gloading, sloading, nitem) {
    sum.general.loading <- sum(gloading)
    sum.specific.loading <- sum(sloading)
    sum.res <- (sum.general.loading ^ 2 - rel * (sum.general.loading ^ 2 + sum.specific.loading ^ 2)) / rel
    res <- sum.res / nitem
    return(res)
}

gloading <- rep(factor.loading.g, (item.number.rev + item.number.ben + item.number.avo))
sloading <- rep(factor.loading.s, (item.number.rev + item.number.ben))

forgive.res <- h.omega.res(rel.forgive, gloading, sloading, (item.number.rev + item.number.ben + item.number.avo))

for (k in 1:10) {
    
    ### Correlation between M and Y ###
    rho <- runif(1, -1, 1)
    ### Residual variance-covariance matrix of M and Y ###
    mu <- numeric(2)
    sigma <- matrix(c(1, rho, rho, 1), 2, 2)
    ### Generate residuals for M and Y ###
    residual.MY <- mvrnorm(N, mu, sigma)

    ### Calculate latent variable M ###
    M <- coef.a * X + residual.MY[, 1]

    ### Calcualte latent variable Y; also, this is the general factor of the bifactor model ###
    Y <- coef.b * M + residual.MY[, 2]
    
    ### Gnerate specific factors ###
    Y.ben <- rnorm(N, 0, 1)
    Y.rev <- rnorm(N, 0, 1)

    Y.bi <- c(Y, Y.ben, Y.rev) %>% matrix(N, 3) # Create a matrix including all factors of the bifactor model

    ### Reliability of relationship value: .95 ###
    ### Mean matrix of meausrement error for M ###
    mu.measuremenet.error.M <- rep(0, item.number.rv)
    ### Variance-covariance matrix of measurement error for M ###
    var.measuremenet.error.M <- diag(rv.res, item.number.rv)
    ### Generate residual data for observed variables of M
    obs.residual.M <- mvrnorm(N, mu.measuremenet.error.M, var.measuremenet.error.M)
    ### Factor loading matrix ###
    factor.loading.M <- rep(factor.loading, item.number.rv)
    ### Generate continous observed data for M ###
    obs.M.total[[k]] <- M %*% t(factor.loading.M) + obs.residual.M

    ### Reliability of forgiveness: .901 ###
    ### Mean matrix of meausrement error for Y ###
    mu.measuremenet.error.Y <- rep(0, item.number.forgive)
    ### Variance-covariance matrix of measurement error for Y ###
    var.measuremenet.error.Y <- diag(forgive.res, item.number.forgive)
    ### Generate residual data for observed variables of Y
    obs.residual.Y <- mvrnorm(N, mu.measuremenet.error.Y, var.measuremenet.error.Y)
    ### Factor loading matrix ###
    factor.loading.forgive <- rep(factor.loading.g, item.number.forgive)
    factor.loading.rev <- c(numeric(item.number.avo), rep(factor.loading.s, item.number.rev), numeric(item.number.ben))
    factor.loading.ben <- c(numeric(item.number.forgive - item.number.ben), rep(factor.loading.s, item.number.ben))
    factor.loading.Y <- c(factor.loading.forgive, factor.loading.rev, factor.loading.ben) %>% matrix(item.number.forgive, 3)
    ### Generate continous observed data for Y ###
    obs.Y.total[[k]] <- Y.bi %*% t(factor.loading.Y) + obs.residual.Y
}

obs.M <- matrix(NA, N, ncol(obs.M.total[[1]]))
obs.Y <- matrix(NA, N, ncol(obs.Y.total[[1]]))
### Average obs.M and obs.Y ###
for (i in 1:ncol(obs.M.total[[1]])) {
    result.M <- lapply(obs.M.total, "[", , i)
    obs.M[, i] <- rowMeans(do.call(cbind, result.M))
}

for (i in 1:ncol(obs.Y.total[[1]])) {
    result.Y <- lapply(obs.Y.total, "[", , i)
    obs.Y[, i] <- rowMeans(do.call(cbind, result.Y))
}

### Generate ordinal data for M ### 
ordinal.obs.M <- matrix(NA, N, item.number.rv)
for (i in 1:item.number.rv) {
    cutoff1.M <- quantile(obs.M[, i], cut.point.M[1])
    cutoff2.M <- quantile(obs.M[, i], cut.point.M[2])
    cutoff3.M <- quantile(obs.M[, i], cut.point.M[3])
    cutoff4.M <- quantile(obs.M[, i], cut.point.M[4])
    cutoff5.M <- quantile(obs.M[, i], cut.point.M[5])
    cutoff6.M <- quantile(obs.M[, i], cut.point.M[6])
    cutoff7.M <- quantile(obs.M[, i], cut.point.M[7])
    cutoff8.M <- quantile(obs.M[, i], cut.point.M[8])
    ordinal.obs.M[, i] <- cut(obs.M[, i], breaks = c(min(obs.M[, i]) - 1, 
    cutoff1.M, cutoff2.M, cutoff3.M, cutoff4.M, cutoff5.M, cutoff6.M, cutoff7.M, cutoff8.M,
    max(obs.M[, i]) + 1), labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9))
}

### Generate ordinal data for Y ### 
ordinal.obs.Y <- matrix(NA, N, item.number.forgive)
for (i in 1:item.number.forgive) {
    cutoff1.Y <- quantile(obs.Y[, i], cut.point.Y[1])
    cutoff2.Y <- quantile(obs.Y[, i], cut.point.Y[2])
    cutoff3.Y <- quantile(obs.Y[, i], cut.point.Y[3])
    cutoff4.Y <- quantile(obs.Y[, i], cut.point.Y[4])
    ordinal.obs.Y[, i] <- cut(obs.Y[, i], breaks = c(min(obs.Y[, i]) - 1, 
    cutoff1.Y, cutoff2.Y, cutoff3.Y, cutoff4.Y, max(obs.Y[, i]) + 1), labels = c(1, 2, 3, 4, 5))
}

### Dataset with dummy variable
apology <- cbind(X, ordinal.obs.M, ordinal.obs.Y)
filename <- "C:/Users/wqemi/OneDrive - Florida State University/MyFSU_OneDrive/Documents/FSU/moderated mediation/sensitivity analysis/apology/ap_analysis_average.dat"
write.table(apology, filename, col.names = FALSE, row.names = FALSE, quote = FALSE)

filename <- "C:/Users/wqemi/OneDrive - Florida State University/MyFSU_OneDrive/Documents/FSU/moderated mediation/sensitivity analysis/apology/ap_analysis_average.csv"
write.table(apology, filename, col.names = TRUE, row.names = FALSE, quote = FALSE)