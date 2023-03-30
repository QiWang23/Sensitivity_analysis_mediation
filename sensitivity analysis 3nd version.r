library(base)
library(data.table)
library(magrittr)
library(shiny)
library(shinyWidgets)
library(lavaan)
library(stringr)
library(ggplot2)
library(latex2exp)
library(rsconnect)

### Record warning message and error ###
myTryCatch <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
}

###### Function part for sensitivity analysis ######
###### Create dataset for plot ######
################ Delta method #################
###### Function includes X,M and Y ######
Analyze.XMY <- function(dataset, method) {
  
  # split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M and Y
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))
  
  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")

  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r * M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1 * M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1 * Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y and covariates for M######
Analyze.XMYCM <- function(dataset, method) {
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y, and CM
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CM <- sum(temp == "CM")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CM) {
      command <- paste0(command, " + ", "CM.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y and covariates for Y######
Analyze.XMYCY <- function(dataset, method) {
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y, and CY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CY <- sum(temp == "CY")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")

  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CY) {
      command <- paste0(command, " + ", "CY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y and covariates for M & Y######
Analyze.XMYCMY <- function(dataset, method) {
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y, and CY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CMY <- sum(temp == "CMY")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")

  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y, covariates for M and covariates for M & Y######
Analyze.XMYCMMY <- function(dataset, method) {
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y, CY and CMY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CM <- sum(temp == "CM")
  n.CMY <- sum(temp == "CMY")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))
  
  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")

  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CM) {
      command <- paste0(command, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y, covariates for Y and covariates for M & Y######
Analyze.XMYCYMY <- function(dataset, method) {
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y,CY and CMY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CY <- sum(temp == "CY")
  n.CMY <- sum(temp == "CMY")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")

  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CY) {
      command <- paste0(command, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y, covariates for M and covariates for Y######
Analyze.XMYCMCY <- function(dataset, method) {
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y,CM and CY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CM <- sum(temp == "CM")
  n.CY <- sum(temp == "CY")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")

  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CM) {
      command <- paste0(command, " + ", "CM.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CY) {
      command <- paste0(command, " + ", "CY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y, covariates for M, Y and M & Y######
Analyze.all <- function(dataset, method) {
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y, CM, CY and CMY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CM <- sum(temp == "CM")
  n.CY <- sum(temp == "CY")
  n.CMY <- sum(temp == "CMY")
  
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CM) {
      command <- paste0(command, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CY) {
      command <- paste0(command, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # } 

    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}

################ Bootstrap method #################
###### Function includes X,M and Y ######
Analyze.XMY.boot <- function(dataset, bootnum) {
  
  method <- "bootstrap"
  
  # split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M and Y
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {

      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y and covariates for M######
Analyze.XMYCM.boot <- function(dataset, bootnum) {
  
  method <- "bootstrap"
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y, and CM
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CM <- sum(temp == "CM")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CM) {
      command <- paste0(command, " + ", "CM.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {

      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y and covariates for Y######
Analyze.XMYCY.boot <- function(dataset, bootnum) {
  
  method <- "bootstrap"
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y, and CY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CY <- sum(temp == "CY")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CY) {
      command <- paste0(command, " + ", "CY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {

      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y and covariates for M & Y######
Analyze.XMYCMY.boot <- function(dataset, bootnum) {
  
  method <- "bootstrap"
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y, and CY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CMY <- sum(temp == "CMY")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {

      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y, covariates for M and covariates for M & Y######
Analyze.XMYCMMY.boot <- function(dataset, bootnum) {
  
  method <- "bootstrap"
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y, CY and CMY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CM <- sum(temp == "CM")
  n.CMY <- sum(temp == "CMY")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")

  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CM) {
      command <- paste0(command, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {

      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y, covariates for Y and covariates for M & Y######
Analyze.XMYCYMY.boot <- function(dataset, bootnum) {
  
  method <- "bootstrap"
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y,CY and CMY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CY <- sum(temp == "CY")
  n.CMY <- sum(temp == "CMY")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CY) {
      command <- paste0(command, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {

      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y, covariates for M and covariates for Y######
Analyze.XMYCMCY.boot <- function(dataset, bootnum) {
  
  method <- "bootstrap"
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y,CM and CY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CM <- sum(temp == "CM")
  n.CY <- sum(temp == "CY")
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CM) {
      command <- paste0(command, " + ", "CM.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CY) {
      command <- paste0(command, " + ", "CY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {

      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}
###########################################################

###### Function includes X, M, Y, covariates for M, Y and M & Y######
Analyze.all.boot <- function(dataset, bootnum) {
  
  method <- "bootstrap"
  
  # Split numbers and characters
  temp <- str_extract(names(dataset), "[aA-zZ]+")
  
  
  # Get N of M, Y, CM, CY and CMY
  n.M <- sum(temp == "M")
  n.Y <- sum(temp == "Y")
  n.CM <- sum(temp == "CM")
  n.CY <- sum(temp == "CY")
  n.CMY <- sum(temp == "CMY")
  
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    # Measurement model
    command <- character(0)
    command <- paste0(command, "M =~ NA * M.1")
    for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }
    
    # Regression model
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~ alpha*X")
    for (asd in 1:n.CM) {
      command <- paste0(command, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~ gamma*X + beta*M")
    for (asd in 1:n.CY) {
      command <- paste0(command, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
      command <- paste0(command, " + ", "CMY.", asd)
    }
    
    # Residual correlations
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ r*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "M ~~ 1*M")
    command <- paste0(command, "\n")
    command <- paste0(command, "Y ~~ 1*Y")
    
    # Mediation effect
    command <- paste0(command, "\n")
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {

      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}

##################################################################################### Bifactor model #################################################################################################
############################## M is the bifactor model, but Y is not #############################################
############################## Delta method ########################################
Analyze.bi.M <- function(dataset, mediator, method) {
  
  # split numbers and characters 
  temp <- str_extract(names(dataset), "[aA-zZ]+") # Get characters
  
  # Get N of X, M and Y
  n.X <- sum(temp == "X") ### Total number of X
  n.OM <- sum(temp == "M") ### Total number of M: N of the factors that only load on the general factor
  n.M <- sum(temp == "MS") ### Total number of MS

  ### Get N of specific factors of M ###
  str.MS <- str_split(names(dataset), pattern = " ") %>% unlist
  str.MS.2 <- strsplit(str.MS, "(?=[.])", perl = TRUE)
  str.MS.3 <- lapply(str.MS.2, "[", 1) %>% unlist
  
  str.MS.4 <- str.MS.3[grep("MS", str.MS.3)] ### Save specific factors of M (save variables which start with MS)
  MS.num <- gsub(".*?([0-9]+).*", "\\1", str.MS.4) %>% as.numeric() %>% unique ### Save numeric parts of string from MS

  n.Y <- sum(temp == "Y") ### Total number of Y
  n.CM <- sum(temp == "CM") ### Total number of CM
  n.CY <- sum(temp == "CY") ### Total number of CY
  n.CMY <- sum(temp == "CMY") ### Total number of CMY
  
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    if (n.OM == 0) {
      ### Perfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of M ###
      command <- character(0)

      command.MS.g.start <- list()
      command.MS.g.final <- list()
      n.MS <- list()

      ### Extract N of each MS ###
      for (n in 1:length(MS.num)) {
        n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
      }

      n.MS.num <- 1:length(n.MS) # Get length of MS

      ### Determine the factor of the bifactor model is the mediator in the model ###
      if (mediator %>% as.numeric() == 0) {
        mediator.focus <- "MG"
      } else if (mediator %>% as.numeric() != 0) {
        for (n in 1:length(n.MS.num)) {
          if (mediator %>% as.numeric() == n.MS.num[n]) {
            mediator.focus <- paste0("MS", n.MS.num[n])
          }
        }
      }

      ### Create the start syntax of MS: MS1., MS2., MS3. ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g.start[[n]] <- paste0("MS", n.MS.num[n], ".1")
      }

      ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g <- command.MS.g.start[[n]]
        for (asd in 1:(n.MS[[n]] - 1)) {
          command.MS.g <- paste0(command.MS.g, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
        }
        command.MS.g.final[[n]] <- command.MS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of M ###
      command <- paste0("MG =~ NA * ", command.MS.g.final[[1]])
      for (n in 2:length(command.MS.g.final)) {
        command <- paste0(command, " + ", command.MS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    } else if (n.OM != 0) {
      ### Imperfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of M ###
      command <- character(0)

      command.MS.g.start <- list()
      command.MS.g.final <- list()
      n.MS <- list()

      ### Extract N of each MS ###
      for (n in 1:length(MS.num)) {
        n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
      }

      n.MS.num <- 1:length(n.MS) # Get length of MS

      ### Determine the factor of the bifactor model is the mediator in the model ###
      if (mediator %>% as.numeric() == 0) {
        mediator.focus <- "MG"
      } else if (mediator %>% as.numeric() != 0) {
        for (n in 1:length(n.MS.num)) {
          if (mediator %>% as.numeric() == n.MS.num[n]) {
            mediator.focus <- paste0("MS", n.MS.num[n])
          }
        }
      }

      ### Create the start syntax of MS: MS1., MS2., MS3. ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g.start[[n]] <- paste0("MS", n.MS.num[n], ".1")
      }

      ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g <- command.MS.g.start[[n]]
        for (asd in 1:(n.MS[[n]] - 1)) {
          command.MS.g <- paste0(command.MS.g, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
        }
        command.MS.g.final[[n]] <- command.MS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of M ###
      command <- paste0("MG =~ NA * M.1")
      for (asd in 2:n.OM) {
        command <- paste0(command, " + ", "M.", asd)
      }

      for (n in 1:length(command.MS.g.final)) {
        command <- paste0(command, " + ", command.MS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    }
    ### Write syntax for the specific factors of M ###
    n.MS <- list()

    for (n in 1:length(MS.num)) {
      n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
    }

    command.MS.start <- list()
    command.MS <- character(0)
    command.MS.final <- list()

    ### Create the start syntax of MS: MS11, MS21, MS31 ###
    for (n in 1:length(n.MS.num)) {
      command.MS.start[[n]] <- paste0(command.MS, "MS", n.MS.num[n], " =~ NA * MS", n.MS.num[n], ".1")
    }

    ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
    for (n in 1:length(n.MS.num)) {
      command.MS <- command.MS.start[[n]]
      for (asd in 1:(n.MS[[n]] - 1)) {
        command.MS <- paste0(command.MS, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
      }
      command.MS.final[[n]] <- paste0(command.MS, "\n")
    }

    ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for specific factors of M ###
    command <- paste0(command, command.MS.final[[1]])
    for (n in 2:length(command.MS.final)) {
      command <- paste0(command, command.MS.final[[n]])
    }

    ### Write syntax for Y ###
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }

    command <- paste0(command, "\n")

    # Regression model
    if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY == 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY == 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      
      for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
      }
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY == 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
      }
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY != 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
      }
      for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
      }      
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY != 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
      }
      for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
      }
      for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
      }
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY != 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
      }
      for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
      }
      for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
      }
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY == 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
      }
      
      for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
      }
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY != 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
      }
      for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
      }
      command.m.new <- paste0(command.m.new, "\n")
      
      for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
      }
      for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
      }
      command.y.new <- paste0(command.y.new, "\n")
    }

    # Residual correlations
    command.r <- character(0) ### Create a new character to save the syntax for residuals
    command.r <- paste0(command.r, "Y ~~ r*M")
    command.r <- paste0(command.r, "\n")
    command.r <- paste0(command.r, "M ~~ 1*M") # For correlation between M and Y
    command.r <- paste0(command.r, "\n") 
    command.r <- paste0(command.r, "Y ~~ 1*Y", "\n") # For correlation between M and Y

    ### Replace mediator M according to researcher's choice ###
    command.r.new <- gsub("M", mediator.focus, command.r)
    command <- paste0(command, command.m.new, command.y.new, command.r.new)
   
    ### Unit variance identification (UVI) ###
    if (mediator == 0) {
      ### Mediator is the general factor ###
      for (asd in 1:length(n.MS.num)) {
        command <- paste0(command, "MS", asd, " ~~ 1 * MS", asd, "\n")
      }
    } else if (mediator != 0) {
      ### Mediator is the specific factor ###
      command <- paste0(command, "MG ~~ 1 * MG", "\n")
      n.MS.num.new <- n.MS.num[n.MS.num != mediator]
      for (asd in 1:length(n.MS.num.new)) {
        command <- paste0(command, "MS", n.MS.num.new[asd], " ~~ 1 * MS", n.MS.num.new[asd], "\n")
      }
    }

    ### Residual correlations of bifactor model ###
    ### Get correlations between the general factor and specific factors ###
    res.cor.MS <- function(x) {
      cor.MS <- paste0("MG ~~ 0 * MS", x)
      return(cor.MS)
    }

    cor.gs <- lapply(n.MS.num %>% as.list, res.cor.MS)

    command <- paste0(command, cor.gs[[1]], "\n")

    for (n in 2:length(cor.gs)) {
      command <- paste0(command, cor.gs[[n]], "\n")
    }

    ### Get correlations between sepcific factors ###

    com.ss <- combn(paste0("MS", 1:length(n.MS.num)), 2) # Get all combinations of specific factors

    cor.ss <- list()

    ### Get syntax of correlations between specific factors ###
    for (n in 1:ncol(com.ss)) {
      cor.ss[[n]] <- paste0(com.ss[1, n], " ~~ 0 * ", com.ss[2, n])
    }

    if (ncol(com.ss) == 1) {
      
      command <- paste0(command, cor.ss[[1]], "\n")
      
    } else if (ncol(com.ss) > 1) {
      
      command <- paste0(command, cor.ss[[1]], "\n")
      
      for (n in 2:length(cor.ss)) {
        command <- paste0(command, cor.ss[[n]], "\n")
      }
    }
      
    # Mediation effect
    command <- paste0(command, "ind := alpha*beta")
  
    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }

    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}
###########################################################

######################################################## Bootstrap #################################################################
Analyze.bi.M.boot <- function(dataset, mediator, bootnum) {
  
  method <- "bootstrap"

  # split numbers and characters 
  temp <- str_extract(names(dataset), "[aA-zZ]+") # Get characters
  
  # Get N of X, M and Y
  n.X <- sum(temp == "X") ### Total number of X
  n.OM <- sum(temp == "M") ### Total number of M: N of the factors that only load on the general factor
  n.M <- sum(temp == "MS") ### Total number of MS

  ### Get N of specific factors of M ###
  str.MS <- str_split(names(dataset), pattern = " ") %>% unlist
  str.MS.2 <- strsplit(str.MS, "(?=[.])", perl = TRUE)
  str.MS.3 <- lapply(str.MS.2, "[", 1) %>% unlist
  
  str.MS.4 <- str.MS.3[grep("MS", str.MS.3)] ### Save specific factors of M (save variables which start with MS)
  MS.num <- gsub(".*?([0-9]+).*", "\\1", str.MS.4) %>% as.numeric() %>% unique ### Save numeric parts of string from MS

  n.Y <- sum(temp == "Y") ### Total number of Y
  n.CM <- sum(temp == "CM") ### Total number of CM
  n.CY <- sum(temp == "CY") ### Total number of CY
  n.CMY <- sum(temp == "CMY") ### Total number of CMY
  
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    if (n.OM == 0) {
      ### Perfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of M ###
      command <- character(0)

      command.MS.g.start <- list()
      command.MS.g.final <- list()
      n.MS <- list()

      ### Extract N of each MS ###
      for (n in 1:length(MS.num)) {
        n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
      }

      n.MS.num <- 1:length(n.MS) # Get length of MS

      ### Determine the factor of the bifactor model is the mediator in the model ###
      if (mediator %>% as.numeric() == 0) {
        mediator.focus <- "MG"
      } else if (mediator %>% as.numeric() != 0) {
        for (n in 1:length(n.MS.num)) {
          if (mediator %>% as.numeric() == n.MS.num[n]) {
            mediator.focus <- paste0("MS", n.MS.num[n])
          }
        }
      }

      ### Create the start syntax of MS: MS1., MS2., MS3. ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g.start[[n]] <- paste0("MS", n.MS.num[n], ".1")
      }

      ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g <- command.MS.g.start[[n]]
        for (asd in 1:(n.MS[[n]] - 1)) {
          command.MS.g <- paste0(command.MS.g, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
        }
        command.MS.g.final[[n]] <- command.MS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of M ###
      command <- paste0("MG =~ 1 * ", command.MS.g.final[[1]])
      for (n in 2:length(command.MS.g.final)) {
        command <- paste0(command, " + ", command.MS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    } else if (n.OM != 0) {
      ### Imperfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of M ###
      command <- character(0)

      command.MS.g.start <- list()
      command.MS.g.final <- list()
      n.MS <- list()

      ### Extract N of each MS ###
      for (n in 1:length(MS.num)) {
        n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
      }

      n.MS.num <- 1:length(n.MS) # Get length of MS

      ### Determine the factor of the bifactor model is the mediator in the model ###
      if (mediator %>% as.numeric() == 0) {
        mediator.focus <- "MG"
      } else if (mediator %>% as.numeric() != 0) {
        for (n in 1:length(n.MS.num)) {
          if (mediator %>% as.numeric() == n.MS.num[n]) {
            mediator.focus <- paste0("MS", n.MS.num[n])
          }
        }
      }

      ### Create the start syntax of MS: MS1., MS2., MS3. ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g.start[[n]] <- paste0("MS", n.MS.num[n], ".1")
      }

      ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g <- command.MS.g.start[[n]]
        for (asd in 1:(n.MS[[n]] - 1)) {
          command.MS.g <- paste0(command.MS.g, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
        }
        command.MS.g.final[[n]] <- command.MS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of M ###
      command <- paste0("MG =~ NA * M.1")
      for (asd in 2:n.OM) {
        command <- paste0(command, " + ", "M.", asd)
      }

      for (n in 1:length(command.MS.g.final)) {
        command <- paste0(command, " + ", command.MS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    }
    ### Write syntax for the specific factors of M ###
    n.MS <- list()

    for (n in 1:length(MS.num)) {
      n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
    }

    command.MS.start <- list()
    command.MS <- character(0)
    command.MS.final <- list()

    ### Create the start syntax of MS: MS11, MS21, MS31 ###
    for (n in 1:length(n.MS.num)) {
      command.MS.start[[n]] <- paste0(command.MS, "MS", n.MS.num[n], " =~ NA * MS", n.MS.num[n], ".1")
    }

    ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
    for (n in 1:length(n.MS.num)) {
      command.MS <- command.MS.start[[n]]
      for (asd in 1:(n.MS[[n]] - 1)) {
        command.MS <- paste0(command.MS, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
      }
      command.MS.final[[n]] <- paste0(command.MS, "\n")
    }

    ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for specific factors of M ###
    command <- paste0(command, command.MS.final[[1]])
    for (n in 2:length(command.MS.final)) {
      command <- paste0(command, command.MS.final[[n]])
    }

    ### Write syntax for Y ###
    command <- paste0(command, "Y =~ NA * Y.1")
    for (asd in 1:(n.Y - 1)) {
      command <- paste0(command, " + ", "Y.", asd + 1)
    }

    command <- paste0(command, "\n")

    # Regression model
    if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY == 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY == 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      
      for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
      }
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY == 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
      }
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY != 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
      }
      for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
      }      
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY != 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
      }
      for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
      }
      for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
      }
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY != 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
      }
      for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
      }
      for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
      }
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY == 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
      }
      
      for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
      }
      
      command.m.new <- paste0(command.m.new, "\n")
      command.y.new <- paste0(command.y.new, "\n")
      
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY != 0) {
      
      command.m <- character(0) ### Create a new character to save the syntax for M
      command.y <- character(0) ### Create a new character to save the syntax for Y
      command.m <- paste0(command.m, "M ~ alpha*X")
      command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
      
      ### Replace mediator M according to researcher's choice ###
      command.m.new <- gsub("M", mediator.focus, command.m)
      command.y.new <- gsub("M", mediator.focus, command.y)
      
      for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
      }
      for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
      }
      command.m.new <- paste0(command.m.new, "\n")
      
      for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
      }
      for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
      }
      command.y.new <- paste0(command.y.new, "\n")
    }

    # Residual correlations
    command.r <- character(0) ### Create a new character to save the syntax for residuals
    command.r <- paste0(command.r, "Y ~~ r*M")
    command.r <- paste0(command.r, "\n")
    command.r <- paste0(command.r, "M ~~ 1*M") # For correlation between M and Y
    command.r <- paste0(command.r, "\n") 
    command.r <- paste0(command.r, "Y ~~ 1*Y", "\n") # For correlation between M and Y

    ### Replace mediator M according to researcher's choice ###
    command.r.new <- gsub("M", mediator.focus, command.r)
    command <- paste0(command, command.m.new, command.y.new, command.r.new)

    ### Unit variance identification (UVI) ###
    if (mediator == 0) {
      ### Mediator is the general factor ###
      for (asd in 1:length(n.MS.num)) {
        command <- paste0(command, "MS", asd, " ~~ 1 * MS", asd, "\n")
      }
    } else if (mediator != 0) {
      ### Mediator is the specific factor ###
      command <- paste0(command, "MG ~~ 1 * MG", "\n")
      n.MS.num.new <- n.MS.num[n.MS.num != mediator]
      for (asd in 1:length(n.MS.num.new)) {
        command <- paste0(command, "MS", n.MS.num.new[asd], " ~~ 1 * MS", n.MS.num.new[asd], "\n")
      }
    }

    ### Residual correlations of bifactor model ###
    ### Get correlations between the general factor and specific factors ###
    res.cor.MS <- function(x) {
      cor.MS <- paste0("MG ~~ 0 * MS", x)
      return(cor.MS)
    }

    cor.gs <- lapply(n.MS.num %>% as.list, res.cor.MS)

    command <- paste0(command, cor.gs[[1]], "\n")

    for (n in 2:length(cor.gs)) {
      command <- paste0(command, cor.gs[[n]], "\n")
    }

    ### Get correlations between sepcific factors ###

    com.ss <- combn(paste0("MS", 1:length(n.MS.num)), 2) # Get all combinations of specific factors

    cor.ss <- list()

    ### Get syntax of correlations between specific factors ###
    for (n in 1:ncol(com.ss)) {
      cor.ss[[n]] <- paste0(com.ss[1, n], " ~~ 0 * ", com.ss[2, n])
    }

    if (ncol(com.ss) == 1) {
      
      command <- paste0(command, cor.ss[[1]], "\n")
      
    } else if (ncol(com.ss) > 1) {
      
      command <- paste0(command, cor.ss[[1]], "\n")
      
      for (n in 2:length(cor.ss)) {
        command <- paste0(command, cor.ss[[n]], "\n")
      }
    }
      
    # Mediation effect
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {
      ### If there is only one warning ###
      ### If warnings include "bootstrap runs resulted in nonadmissible solutions." ###
      ### or "bootstrap runs failed or did not converge." ###
      ### Then 1 will be given to indicate these warnings ###
      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      ### If there are more than 1 warning message ###
      ### Check warning message one by one ###
      ### war.temp is created to store values that warnings include ###
      ### "bootstrap runs resulted in nonadmissible solutions." ###
      ### or "bootstrap runs failed or did not converge." ###
      ### If all warning messages are "bootstrap runs resulted in nonadmissible ###
      ### solutions." or "bootstrap runs failed or did not converge. ###
      ### Then, 1 is given to indicate this situation. ###
      ### If there are some other warnings that are not belonged to previous warnings, ###
      ### 2 is given to indicate this. ###
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}
###########################################################

############################################## Y is the bifactor model, but M is not #######################################################
############################################## Delta method #######################################################
Analyze.bi.Y <- function(dataset, response, method) {
  
  # split numbers and characters 
  temp <- str_extract(names(dataset), "[aA-zZ]+") # Get characters
  
  # Get N of X, M and Y
  n.X <- sum(temp == "X") ### Total number of X
  n.M <- sum(temp == "M") ### Total number of MS; Also, this is the number of the general factor of M

  ### Get N of specific factors of Y ###
  str.YS <- str_split(names(dataset), pattern = " ") %>% unlist
  str.YS.2 <- strsplit(str.YS, "(?=[.])", perl = TRUE)
  str.YS.3 <- lapply(str.YS.2, "[", 1) %>% unlist
  
  str.YS.4 <- str.YS.3[grep("YS", str.YS.3)] ### Save specific factors of M (save variables which start with YS)
  YS.num <- gsub(".*?([0-9]+).*", "\\1", str.YS.4) %>% as.numeric() %>% unique ### Save numeric parts of string from YS

  n.Y <- sum(temp == "YS") ### Total number of YS
  n.OY <- sum(temp == "Y") ### Total number of Y: N of the factors that only load on the general factor
  n.CM <- sum(temp == "CM") ### Total number of CM
  n.CY <- sum(temp == "CY") ### Total number of CY
  n.CMY <- sum(temp == "CMY") ### Total number of CMY
  
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")

  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    if (n.OY == 0) {
      ### Perfect bifactor structure ###
      # Measurement model
      ### Write syntax for M ###
      command <- character(0)
      command <- paste0(command, "M =~ NA * M.1")
      for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
      }

      command <- paste0(command, "\n")

      ### Write syntax for the general factor of Y ###
      command.YS.g.start <- list()
      command.YS.g.final <- list()
      n.YS <- list()

      ### Extract N of each YS ###
      for (n in 1:length(YS.num)) {
      n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
      }

      n.YS.num <- 1:length(n.YS) # Get length of YS

      ### Determine the factor of the bifactor model is the response in the model ###
      if (response %>% as.numeric() == 0) {
      response.focus <- "YG"
      } else if (response %>% as.numeric() != 0) {
      for (n in 1:length(n.YS.num)) {
          if (response %>% as.numeric() == n.YS.num[n]) {
          response.focus <- paste0("YS", n.YS.num[n])
          }
      }
      }

      ### Create the start syntax of YS: YS1., YS2., YS3. ###
      for (n in 1:length(n.YS.num)) {
      command.YS.g.start[[n]] <- paste0("YS", n.YS.num[n], ".1")
      }

      ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
      for (n in 1:length(n.YS.num)) {
      command.YS.g <- command.YS.g.start[[n]]
      for (asd in 1:(n.YS[[n]] - 1)) {
          command.YS.g <- paste0(command.YS.g, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
      }
      command.YS.g.final[[n]] <- command.YS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of Y ###
      command <- paste0(command, "YG =~ 1 * ", command.YS.g.final[[1]])
      for (n in 2:length(command.YS.g.final)) {
      command <- paste0(command, " + ", command.YS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    } else if (n.OY != 0) {
      ### ImPerfect bifactor structure ###
      # Measurement model
      ### Write syntax for M ###
      command <- character(0)
      command <- paste0(command, "M =~ NA * M.1")
      for (asd in 1:(n.M - 1)) {
        command <- paste0(command, " + ", "M.", asd + 1)
      }

      command <- paste0(command, "\n")

      ### Write syntax for the general factor of Y ###
      command.YS.g.start <- list()
      command.YS.g.final <- list()
      n.YS <- list()

      ### Extract N of each YS ###
      for (n in 1:length(YS.num)) {
        n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
      }

      n.YS.num <- 1:length(n.YS) # Get length of YS

      ### Determine the factor of the bifactor model is the response in the model ###
      if (response %>% as.numeric() == 0) {
        response.focus <- "YG"
      } else if (response %>% as.numeric() != 0) {
        for (n in 1:length(n.YS.num)) {
          if (response %>% as.numeric() == n.YS.num[n]) {
            response.focus <- paste0("YS", n.YS.num[n])
          }
        }
      }

      ### Create the start syntax of YS: YS1., YS2., YS3. ###
      for (n in 1:length(n.YS.num)) {
        command.YS.g.start[[n]] <- paste0("YS", n.YS.num[n], ".1")
      }

      ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
      for (n in 1:length(n.YS.num)) {
        command.YS.g <- command.YS.g.start[[n]]
        for (asd in 1:(n.YS[[n]] - 1)) {
          command.YS.g <- paste0(command.YS.g, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
        }
        command.YS.g.final[[n]] <- command.YS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of Y ###
      command <- paste0(command, "YG =~ NA * Y.1")
      for (asd in 2:n.OY) {
        command <- paste0(command, " + ", "Y.", asd)
      }
      for (n in 1:length(command.YS.g.final)) {
        command <- paste0(command, " + ", command.YS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    }

    ### Write syntax for the specific factors of Y ###
    n.YS <- list()

    for (n in 1:length(YS.num)) {
    n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
    }

    command.YS.start <- list()
    command.YS <- character(0)
    command.YS.final <- list()

    ### Create the start syntax of YS: YS11, YS21, YS31 ###
    for (n in 1:length(n.YS.num)) {
    command.YS.start[[n]] <- paste0(command.YS, "YS", n.YS.num[n], " =~ NA * YS", n.YS.num[n], ".1")
    }

    ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
    for (n in 1:length(n.YS.num)) {
    command.YS <- command.YS.start[[n]]
    for (asd in 1:(n.YS[[n]] - 1)) {
        command.YS <- paste0(command.YS, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
    }
    command.YS.final[[n]] <- paste0(command.YS, "\n")
    }

    ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for specific factors of Y ###
    command <- paste0(command, command.YS.final[[1]])
    for (n in 2:length(command.YS.final)) {
    command <- paste0(command, command.YS.final[[n]])
    }

    # Regression model
    if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    command.m.new <- command.m
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    ### Replace response Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y)
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }      
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    command.m.new <- paste0(command.m.new, "\n")
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    command.y.new <- paste0(command.y.new, "\n")
    }

    # Residual correlations
    command.r <- character(0) ### Create a new character to save the syntax for residuals
    command.r <- paste0(command.r, "Y ~~ r*M")
    command.r <- paste0(command.r, "\n")
    command.r <- paste0(command.r, "M ~~ 1*M") # For correlation between M and Y
    command.r <- paste0(command.r, "\n") 
    command.r <- paste0(command.r, "Y ~~ 1*Y", "\n") # For correlation between M and Y

    ### Replace response Y according to researcher's choice ###
    command.r.new <- gsub("Y", response.focus, command.r)
    command <- paste0(command, command.m.new, command.y.new, command.r.new)

    ### Unit variance identification (UVI) ###
    if (response == 0) {
      ### Response is the general factor ###
      for (asd in 1:length(n.YS.num)) {
        command <- paste0(command, "YS", asd, " ~~ 1 * YS", asd, "\n")
      }
    } else if (response != 0) {
      ### Response is the specific factor ###
      command <- paste0(command, "YG ~~ 1 * YG", "\n")
      n.YS.num.new <- n.YS.num[n.YS.num != response]
      for (asd in 1:length(n.YS.num.new)) {
        command <- paste0(command, "YS", n.YS.num.new[asd], " ~~ 1 * YS", n.YS.num.new[asd], "\n")
      }
    }

    ### Residual correlations of bifactor model ###
    ### Get correlations between the general factor and specific factors ###
    res.cor.YS <- function(x) {
    cor.YS <- paste0("YG ~~ 0 * YS", x)
    return(cor.YS)
    }

    cor.gs <- lapply(n.YS.num %>% as.list, res.cor.YS)

    command <- paste0(command, cor.gs[[1]], "\n")

    for (n in 2:length(cor.gs)) {
    command <- paste0(command, cor.gs[[n]], "\n")
    }

    ### Get correlations between sepcific factors ###

    com.ss <- combn(paste0("YS", 1:length(n.YS.num)), 2) # Get all combinations of specific factors

    cor.ss <- list()

    ### Get syntax of correlations between specific factors ###
    for (n in 1:ncol(com.ss)) {
    cor.ss[[n]] <- paste0(com.ss[1, n], " ~~ 0 * ", com.ss[2, n])
    }

    if (ncol(com.ss) == 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    } else if (ncol(com.ss) > 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    for (n in 2:length(cor.ss)) {
        command <- paste0(command, cor.ss[[n]], "\n")
    }
    }

    # Mediation effect
    command <- paste0(command, "ind := alpha*beta")

    # # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }

    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}
###########################################################

############################################## Bootstrap #######################################################
Analyze.bi.Y.boot <- function(dataset, response, bootnum) {

  method <- "bootstrap"
  
  # split numbers and characters 
  temp <- str_extract(names(dataset), "[aA-zZ]+") # Get characters
  
  # Get N of X, M and Y
  n.X <- sum(temp == "X") ### Total number of X
  n.M <- sum(temp == "M") ### Total number of MS; Also, this is the number of the general factor of M

  ### Get N of specific factors of Y ###
  str.YS <- str_split(names(dataset), pattern = " ") %>% unlist
  str.YS.2 <- strsplit(str.YS, "(?=[.])", perl = TRUE)
  str.YS.3 <- lapply(str.YS.2, "[", 1) %>% unlist
  
  str.YS.4 <- str.YS.3[grep("YS", str.YS.3)] ### Save specific factors of M (save variables which start with YS)
  YS.num <- gsub(".*?([0-9]+).*", "\\1", str.YS.4) %>% as.numeric() %>% unique ### Save numeric parts of string from YS

  n.Y <- sum(temp == "YS") ### Total number of YS
  n.OY <- sum(temp == "Y") ### Total number of Y: N of the factors that only load on the general factor
  n.CM <- sum(temp == "CM") ### Total number of CM
  n.CY <- sum(temp == "CY") ### Total number of CY
  n.CMY <- sum(temp == "CMY") ### Total number of CMY
  
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    if (n.OY == 0) {
      ### Perfect bifactor structure ###
      # Measurement model
      ### Write syntax for M ###
      command <- character(0)
      command <- paste0(command, "M =~ NA * M.1")
      for (asd in 1:(n.M - 1)) {
      command <- paste0(command, " + ", "M.", asd + 1)
      }

      command <- paste0(command, "\n")

      ### Write syntax for the general factor of Y ###
      command.YS.g.start <- list()
      command.YS.g.final <- list()
      n.YS <- list()

      ### Extract N of each YS ###
      for (n in 1:length(YS.num)) {
      n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
      }

      n.YS.num <- 1:length(n.YS) # Get length of YS

      ### Determine the factor of the bifactor model is the response in the model ###
      if (response %>% as.numeric() == 0) {
      response.focus <- "YG"
      } else if (response %>% as.numeric() != 0) {
      for (n in 1:length(n.YS.num)) {
          if (response %>% as.numeric() == n.YS.num[n]) {
          response.focus <- paste0("YS", n.YS.num[n])
          }
      }
      }

      ### Create the start syntax of YS: YS1., YS2., YS3. ###
      for (n in 1:length(n.YS.num)) {
      command.YS.g.start[[n]] <- paste0("YS", n.YS.num[n], ".1")
      }

      ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
      for (n in 1:length(n.YS.num)) {
      command.YS.g <- command.YS.g.start[[n]]
      for (asd in 1:(n.YS[[n]] - 1)) {
          command.YS.g <- paste0(command.YS.g, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
      }
      command.YS.g.final[[n]] <- command.YS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of Y ###
      command <- paste0(command, "YG =~ NA * ", command.YS.g.final[[1]])
      for (n in 2:length(command.YS.g.final)) {
      command <- paste0(command, " + ", command.YS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    } else if (n.OY != 0) {
      ### ImPerfect bifactor structure ###
      # Measurement model
      ### Write syntax for M ###
      command <- character(0)
      command <- paste0(command, "M =~ NA * M.1")
      for (asd in 1:(n.M - 1)) {
        command <- paste0(command, " + ", "M.", asd + 1)
      }

      command <- paste0(command, "\n")

      ### Write syntax for the general factor of Y ###
      command.YS.g.start <- list()
      command.YS.g.final <- list()
      n.YS <- list()

      ### Extract N of each YS ###
      for (n in 1:length(YS.num)) {
        n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
      }

      n.YS.num <- 1:length(n.YS) # Get length of YS

      ### Determine the factor of the bifactor model is the response in the model ###
      if (response %>% as.numeric() == 0) {
        response.focus <- "YG"
      } else if (response %>% as.numeric() != 0) {
        for (n in 1:length(n.YS.num)) {
          if (response %>% as.numeric() == n.YS.num[n]) {
            response.focus <- paste0("YS", n.YS.num[n])
          }
        }
      }

      ### Create the start syntax of YS: YS1., YS2., YS3. ###
      for (n in 1:length(n.YS.num)) {
        command.YS.g.start[[n]] <- paste0("YS", n.YS.num[n], ".1")
      }

      ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
      for (n in 1:length(n.YS.num)) {
        command.YS.g <- command.YS.g.start[[n]]
        for (asd in 1:(n.YS[[n]] - 1)) {
          command.YS.g <- paste0(command.YS.g, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
        }
        command.YS.g.final[[n]] <- command.YS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of Y ###
      command <- paste0(command, "YG =~ NA * Y.1")
      for (asd in 2:n.OY) {
        command <- paste0(command, " + ", "Y.", asd)
      }
      for (n in 1:length(command.YS.g.final)) {
        command <- paste0(command, " + ", command.YS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    }

    ### Write syntax for the specific factors of Y ###
    n.YS <- list()

    for (n in 1:length(YS.num)) {
    n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
    }

    command.YS.start <- list()
    command.YS <- character(0)
    command.YS.final <- list()

    ### Create the start syntax of YS: YS11, YS21, YS31 ###
    for (n in 1:length(n.YS.num)) {
    command.YS.start[[n]] <- paste0(command.YS, "YS", n.YS.num[n], " =~ NA * YS", n.YS.num[n], ".1")
    }

    ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
    for (n in 1:length(n.YS.num)) {
    command.YS <- command.YS.start[[n]]
    for (asd in 1:(n.YS[[n]] - 1)) {
        command.YS <- paste0(command.YS, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
    }
    command.YS.final[[n]] <- paste0(command.YS, "\n")
    }

    ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for specific factors of Y ###
    command <- paste0(command, command.YS.final[[1]])
    for (n in 2:length(command.YS.final)) {
    command <- paste0(command, command.YS.final[[n]])
    }

    # Regression model
    if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    command.m.new <- command.m
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    ### Replace response Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y)
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }      
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace response Y according to researcher's choice ###
    command.m.new <- command.m
    command.y.new <- gsub("Y", response.focus, command.y)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    command.m.new <- paste0(command.m.new, "\n")
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    command.y.new <- paste0(command.y.new, "\n")
    }

    # Residual correlations
    command.r <- character(0) ### Create a new character to save the syntax for residuals
    command.r <- paste0(command.r, "Y ~~ r*M")
    command.r <- paste0(command.r, "\n")
    command.r <- paste0(command.r, "M ~~ 1*M") # For correlation between M and Y
    command.r <- paste0(command.r, "\n") 
    command.r <- paste0(command.r, "Y ~~ 1*Y", "\n") # For correlation between M and Y

    ### Replace response Y according to researcher's choice ###
    command.r.new <- gsub("Y", response.focus, command.r)
    command <- paste0(command, command.m.new, command.y.new, command.r.new)

    ### Unit variance identification (UVI) ###
    if (response == 0) {
      ### Response is the general factor ###
      for (asd in 1:length(n.YS.num)) {
        command <- paste0(command, "YS", asd, " ~~ 1 * YS", asd, "\n")
      }
    } else if (response != 0) {
      ### Response is the specific factor ###
      command <- paste0(command, "YG ~~ 1 * YG", "\n")
      n.YS.num.new <- n.YS.num[n.YS.num != response]
      for (asd in 1:length(n.YS.num.new)) {
        command <- paste0(command, "YS", n.YS.num.new[asd], " ~~ 1 * YS", n.YS.num.new[asd], "\n")
      }
    }

    ### Residual correlations of bifactor model ###
    ### Get correlations between the general factor and specific factors ###
    res.cor.YS <- function(x) {
    cor.YS <- paste0("YG ~~ 0 * YS", x)
    return(cor.YS)
    }

    cor.gs <- lapply(n.YS.num %>% as.list, res.cor.YS)

    command <- paste0(command, cor.gs[[1]], "\n")

    for (n in 2:length(cor.gs)) {
    command <- paste0(command, cor.gs[[n]], "\n")
    }

    ### Get correlations between sepcific factors ###

    com.ss <- combn(paste0("YS", 1:length(n.YS.num)), 2) # Get all combinations of specific factors

    cor.ss <- list()

    ### Get syntax of correlations between specific factors ###
    for (n in 1:ncol(com.ss)) {
    cor.ss[[n]] <- paste0(com.ss[1, n], " ~~ 0 * ", com.ss[2, n])
    }

    if (ncol(com.ss) == 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    } else if (ncol(com.ss) > 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    for (n in 2:length(cor.ss)) {
        command <- paste0(command, cor.ss[[n]], "\n")
    }
    }

    # Mediation effect
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {

      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}
###########################################################

####################################################### Both M and Y are bifactor models ##########################################################
############################## Delta method ####################################
Analyze.bi.MY <- function(dataset, mediator, response, method) {
  
  # split numbers and characters 
  temp <- str_extract(names(dataset), "[aA-zZ]+") # Get characters
  
  # Get N of X, M and Y
  n.X <- sum(temp == "X") ### Total number of X
  n.M <- sum(temp == "MS") ### Total number of MS
  n.OM <- sum(temp == "M") ### Total number of M: N of the factors that only load on the general factor

  ### Get N of specific factors of M ###
  str.MS <- str_split(names(dataset), pattern = " ") %>% unlist
  str.MS.2 <- strsplit(str.MS, "(?=[.])", perl = TRUE)
  str.MS.3 <- lapply(str.MS.2, "[", 1) %>% unlist
  
  str.MS.4 <- str.MS.3[grep("MS", str.MS.3)] ### Save specific factors of M (save variables which start with MS)
  MS.num <- gsub(".*?([0-9]+).*", "\\1", str.MS.4) %>% as.numeric() %>% unique ### Save numeric parts of string from MS

  ### Get N of specific factors of Y ###
  str.YS <- str_split(names(dataset), pattern = " ") %>% unlist
  str.YS.2 <- strsplit(str.YS, "(?=[.])", perl = TRUE)
  str.YS.3 <- lapply(str.YS.2, "[", 1) %>% unlist
  
  str.YS.4 <- str.YS.3[grep("YS", str.YS.3)] ### Save specific factors of M (save variables which start with YS)
  YS.num <- gsub(".*?([0-9]+).*", "\\1", str.YS.4) %>% as.numeric() %>% unique ### Save numeric parts of string from YS

  n.Y <- sum(temp == "YS") ### Total number of Y
  n.OY <- sum(temp == "Y") ### Total number of Y: N of the factors that only load on the general factor
  n.CM <- sum(temp == "CM") ### Total number of CM
  n.CY <- sum(temp == "CY") ### Total number of CY
  n.CMY <- sum(temp == "CMY") ### Total number of CMY
  
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]

    if (n.OM == 0) {
      ### Perfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of M ###
      command <- character(0)

      command.MS.g.start <- list()
      command.MS.g.final <- list()
      n.MS <- list()

      ### Extract N of each MS ###
      for (n in 1:length(MS.num)) {
        n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
      }

      n.MS.num <- 1:length(n.MS) # Get length of MS

      ### Determine the factor of the bifactor model is the mediator in the model ###
      if (mediator %>% as.numeric() == 0) {
        mediator.focus <- "MG"
      } else if (mediator %>% as.numeric() != 0) {
        for (n in 1:length(n.MS.num)) {
          if (mediator %>% as.numeric() == n.MS.num[n]) {
            mediator.focus <- paste0("MS", n.MS.num[n])
          }
        }
      }

      ### Create the start syntax of MS: MS1., MS2., MS3. ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g.start[[n]] <- paste0("MS", n.MS.num[n], ".1")
      }

      ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g <- command.MS.g.start[[n]]
        for (asd in 1:(n.MS[[n]] - 1)) {
          command.MS.g <- paste0(command.MS.g, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
        }
        command.MS.g.final[[n]] <- command.MS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of M ###
      command <- paste0("MG =~ NA * ", command.MS.g.final[[1]])
      for (n in 2:length(command.MS.g.final)) {
        command <- paste0(command, " + ", command.MS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    } else if (n.OM != 0) {
      ### Imperfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of M ###
      command <- character(0)

      command.MS.g.start <- list()
      command.MS.g.final <- list()
      n.MS <- list()

      ### Extract N of each MS ###
      for (n in 1:length(MS.num)) {
        n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
      }

      n.MS.num <- 1:length(n.MS) # Get length of MS

      ### Determine the factor of the bifactor model is the mediator in the model ###
      if (mediator %>% as.numeric() == 0) {
        mediator.focus <- "MG"
      } else if (mediator %>% as.numeric() != 0) {
        for (n in 1:length(n.MS.num)) {
          if (mediator %>% as.numeric() == n.MS.num[n]) {
            mediator.focus <- paste0("MS", n.MS.num[n])
          }
        }
      }

      ### Create the start syntax of MS: MS1., MS2., MS3. ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g.start[[n]] <- paste0("MS", n.MS.num[n], ".1")
      }

      ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g <- command.MS.g.start[[n]]
        for (asd in 1:(n.MS[[n]] - 1)) {
          command.MS.g <- paste0(command.MS.g, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
        }
        command.MS.g.final[[n]] <- command.MS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of M ###
      command <- paste0("MG =~ NA * M.1")
      for (asd in 2:n.OM) {
        command <- paste0(command, " + ", "M.", asd)
      }

      for (n in 1:length(command.MS.g.final)) {
        command <- paste0(command, " + ", command.MS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    }

    ### Write syntax for the specific factors of M ###
    n.MS <- list()

    for (n in 1:length(MS.num)) {
    n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
    }

    command.MS.start <- list()
    command.MS <- character(0)
    command.MS.final <- list()

    ### Create the start syntax of MS: MS11, MS21, MS31 ###
    for (n in 1:length(n.MS.num)) {
    command.MS.start[[n]] <- paste0(command.MS, "MS", n.MS.num[n], " =~ NA * MS", n.MS.num[n], ".1")
    }

    ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
    for (n in 1:length(n.MS.num)) {
    command.MS <- command.MS.start[[n]]
    for (asd in 1:(n.MS[[n]] - 1)) {
        command.MS <- paste0(command.MS, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
    }
    command.MS.final[[n]] <- paste0(command.MS, "\n")
    }

    ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for specific factors of M ###
    command <- paste0(command, command.MS.final[[1]])
    for (n in 2:length(command.MS.final)) {
    command <- paste0(command, command.MS.final[[n]])
    }
    
    if (n.OY == 0) {
      ### Perfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of Y ###
      command.YS.g.start <- list()
      command.YS.g.final <- list()
      n.YS <- list()

      ### Extract N of each YS ###
      for (n in 1:length(YS.num)) {
      n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
      }

      n.YS.num <- 1:length(n.YS) # Get length of YS

      ### Determine the factor of the bifactor model is the response in the model ###
      if (response %>% as.numeric() == 0) {
      response.focus <- "YG"
      } else if (response %>% as.numeric() != 0) {
      for (n in 1:length(n.YS.num)) {
          if (response %>% as.numeric() == n.YS.num[n]) {
          response.focus <- paste0("YS", n.YS.num[n])
          }
      }
      }

      ### Create the start syntax of YS: YS1., YS2., YS3. ###
      for (n in 1:length(n.YS.num)) {
      command.YS.g.start[[n]] <- paste0("YS", n.YS.num[n], ".1")
      }

      ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
      for (n in 1:length(n.YS.num)) {
      command.YS.g <- command.YS.g.start[[n]]
      for (asd in 1:(n.YS[[n]] - 1)) {
          command.YS.g <- paste0(command.YS.g, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
      }
      command.YS.g.final[[n]] <- command.YS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of Y ###
      command <- paste0(command, "YG =~ NA * ", command.YS.g.final[[1]])
      for (n in 2:length(command.YS.g.final)) {
      command <- paste0(command, " + ", command.YS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    } else if (n.OY != 0) {
      ### Imperfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of Y ###
      command.YS.g.start <- list()
      command.YS.g.final <- list()
      n.YS <- list()

      ### Extract N of each YS ###
      for (n in 1:length(YS.num)) {
        n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
      }

      n.YS.num <- 1:length(n.YS) # Get length of YS

      ### Determine the factor of the bifactor model is the response in the model ###
      if (response %>% as.numeric() == 0) {
        response.focus <- "YG"
      } else if (response %>% as.numeric() != 0) {
        for (n in 1:length(n.YS.num)) {
          if (response %>% as.numeric() == n.YS.num[n]) {
            response.focus <- paste0("YS", n.YS.num[n])
          }
        }
      }

      ### Create the start syntax of YS: YS1., YS2., YS3. ###
      for (n in 1:length(n.YS.num)) {
        command.YS.g.start[[n]] <- paste0("YS", n.YS.num[n], ".1")
      }

      ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
      for (n in 1:length(n.YS.num)) {
        command.YS.g <- command.YS.g.start[[n]]
        for (asd in 1:(n.YS[[n]] - 1)) {
          command.YS.g <- paste0(command.YS.g, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
        }
        command.YS.g.final[[n]] <- command.YS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of Y ###
      command <- paste0(command, "YG =~ NA * Y.1")
      for (asd in 2:n.OY) {
        command <- paste0(command, " + ", "Y.", asd)
      }
      for (n in 1:length(command.YS.g.final)) {
        command <- paste0(command, " + ", command.YS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    }

    ### Write syntax for the specific factors of Y ###
    n.YS <- list()

    for (n in 1:length(YS.num)) {
    n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
    }

    command.YS.start <- list()
    command.YS <- character(0)
    command.YS.final <- list()

    ### Create the start syntax of YS: YS11, YS21, YS31 ###
    for (n in 1:length(n.YS.num)) {
    command.YS.start[[n]] <- paste0(command.YS, "YS", n.YS.num[n], " =~ NA * YS", n.YS.num[n], ".1")
    }

    ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
    for (n in 1:length(n.YS.num)) {
    command.YS <- command.YS.start[[n]]
    for (asd in 1:(n.YS[[n]] - 1)) {
        command.YS <- paste0(command.YS, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
    }
    command.YS.final[[n]] <- paste0(command.YS, "\n")
    }

    ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for specific factors of Y ###
    command <- paste0(command, command.YS.final[[1]])
    for (n in 2:length(command.YS.final)) {
    command <- paste0(command, command.YS.final[[n]])
    }

    # Regression model
    if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    
    ### Replace mediator M according to researcher's choice ###
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }      
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    command.m.new <- paste0(command.m.new, "\n")
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    command.y.new <- paste0(command.y.new, "\n")
    }

    # Residual correlations
    command.r <- character(0) ### Create a new character to save the syntax for residuals
    command.r <- paste0(command.r, "Y ~~ r*M")
    command.r <- paste0(command.r, "\n")
    command.r <- paste0(command.r, "M ~~ 1*M") # For correlation between M and Y
    command.r <- paste0(command.r, "\n") 
    command.r <- paste0(command.r, "Y ~~ 1*Y", "\n") # For correlation between M and Y

    ### Replace mediator M according to researcher's choice ###
    command.r.new.1 <- gsub("M", mediator.focus, command.r)
    command.r.new <- gsub("Y", response.focus, command.r.new.1)

    command <- paste0(command, command.m.new, command.y.new, command.r.new)

    ### Unit variance identification (UVI) ###
    if (mediator == 0) {
      ### Mediator is the general factor ###
      for (asd in 1:length(n.MS.num)) {
        command <- paste0(command, "MS", asd, " ~~ 1 * MS", asd, "\n")
      }
    } else if (mediator != 0) {
      ### Mediator is the specific factor ###
      command <- paste0(command, "MG ~~ 1 * MG", "\n")
      n.MS.num.new <- n.MS.num[n.MS.num != mediator]
      for (asd in 1:length(n.MS.num.new)) {
        command <- paste0(command, "MS", n.MS.num.new[asd], " ~~ 1 * MS", n.MS.num.new[asd], "\n")
      }
    }

    if (response == 0) {
      ### Response is the general factor ###
      for (asd in 1:length(n.YS.num)) {
        command <- paste0(command, "YS", asd, " ~~ 1 * YS", asd, "\n")
      }
    } else if (response != 0) {
      ### Response is the specific factor ###
      command <- paste0(command, "YG ~~ 1 * YG", "\n")
      n.YS.num.new <- n.YS.num[n.YS.num != response]
      for (asd in 1:length(n.YS.num.new)) {
        command <- paste0(command, "YS", n.YS.num.new[asd], " ~~ 1 * YS", n.YS.num.new[asd], "\n")
      }
    }

    ### Residual correlations of bifactor model ###
    ### Get correlations between the general factor and specific factors of M ###
    res.cor.MS <- function(x) {
    cor.MS <- paste0("MG ~~ 0 * MS", x)
    return(cor.MS)
    }

    cor.gs <- lapply(n.MS.num %>% as.list, res.cor.MS)

    command <- paste0(command, cor.gs[[1]], "\n")

    for (n in 2:length(cor.gs)) {
    command <- paste0(command, cor.gs[[n]], "\n")
    }

    ### Get correlations between sepcific factors ###

    com.ss <- combn(paste0("MS", 1:length(n.MS.num)), 2) # Get all combinations of specific factors

    cor.ss <- list()

    ### Get syntax of correlations between specific factors ###
    for (n in 1:ncol(com.ss)) {
    cor.ss[[n]] <- paste0(com.ss[1, n], " ~~ 0 * ", com.ss[2, n])
    }

    if (ncol(com.ss) == 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    } else if (ncol(com.ss) > 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    for (n in 2:length(cor.ss)) {
        command <- paste0(command, cor.ss[[n]], "\n")
    }
    }

    ### Residual correlations of bifactor model ###
    ### Get correlations between the general factor and specific factors of Y ###
    res.cor.YS <- function(x) {
    cor.YS <- paste0("YG ~~ 0 * YS", x)
    return(cor.YS)
    }

    cor.gs <- lapply(n.YS.num %>% as.list, res.cor.YS)

    command <- paste0(command, cor.gs[[1]], "\n")

    for (n in 2:length(cor.gs)) {
    command <- paste0(command, cor.gs[[n]], "\n")
    }

    ### Get correlations between sepcific factors ###

    com.ss <- combn(paste0("YS", 1:length(n.YS.num)), 2) # Get all combinations of specific factors

    cor.ss <- list()

    ### Get syntax of correlations between specific factors ###
    for (n in 1:ncol(com.ss)) {
    cor.ss[[n]] <- paste0(com.ss[1, n], " ~~ 0 * ", com.ss[2, n])
    }

    if (ncol(com.ss) == 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    } else if (ncol(com.ss) > 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
      for (n in 2:length(cor.ss)) {
          command <- paste0(command, cor.ss[[n]], "\n")
      }
    }

    # Mediation effect
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, estimator = method))
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
        war.res[i, 1] <- 0
    } else {
        war.res[i, 1] <- 1
    }

    if (is.null(war$error)) {
        war.res[i, 2] <- 0
    } else {
        war.res[i, 2] <- 1
    }

    meffect.est[i] <- tail(parameterEstimates(output),1)$est
    meffect.upper[i] <- tail(parameterEstimates(output),1)$ci.upper
    meffect.lower[i] <- tail(parameterEstimates(output),1)$ci.lower
    
  }
  result <- data.frame(meffect.est, meffect.upper, meffect.lower, rho, war.res)
  return(result)
}
###########################################################

################################################################# Bootstrap ######################################################################
Analyze.bi.MY.boot <- function(dataset, mediator, response, bootnum) {

  method <- "bootstrap"
  
  # split numbers and characters 
  temp <- str_extract(names(dataset), "[aA-zZ]+") # Get characters
  
  # Get N of X, M and Y
  n.X <- sum(temp == "X") ### Total number of X
  n.M <- sum(temp == "MS") ### Total number of MS
  n.OM <- sum(temp == "M") ### Total number of M: N of the factors that only load on the general factor

  ### Get N of specific factors of M ###
  str.MS <- str_split(names(dataset), pattern = " ") %>% unlist
  str.MS.2 <- strsplit(str.MS, "(?=[.])", perl = TRUE)
  str.MS.3 <- lapply(str.MS.2, "[", 1) %>% unlist
  
  str.MS.4 <- str.MS.3[grep("MS", str.MS.3)] ### Save specific factors of M (save variables which start with MS)
  MS.num <- gsub(".*?([0-9]+).*", "\\1", str.MS.4) %>% as.numeric() %>% unique ### Save numeric parts of string from MS

  ### Get N of specific factors of Y ###
  str.YS <- str_split(names(dataset), pattern = " ") %>% unlist
  str.YS.2 <- strsplit(str.YS, "(?=[.])", perl = TRUE)
  str.YS.3 <- lapply(str.YS.2, "[", 1) %>% unlist
  
  str.YS.4 <- str.YS.3[grep("YS", str.YS.3)] ### Save specific factors of M (save variables which start with YS)
  YS.num <- gsub(".*?([0-9]+).*", "\\1", str.YS.4) %>% as.numeric() %>% unique ### Save numeric parts of string from YS

  n.Y <- sum(temp == "YS") ### Total number of Y
  n.OY <- sum(temp == "Y") ### Total number of Y: N of the factors that only load on the general factor
  n.CM <- sum(temp == "CM") ### Total number of CM
  n.CY <- sum(temp == "CY") ### Total number of CY
  n.CMY <- sum(temp == "CMY") ### Total number of CMY
  
  
  rho <- seq(-1, 1, 0.1)
  meffect.est <- numeric(length(rho))
  meffect.upper <- numeric(length(rho))
  meffect.lower <- numeric(length(rho))

  ### Create a matrix to save warnings and errors ###
  war.res <- matrix(NA, length(rho), 2) %>% as.data.frame()
  names(war.res) <- c("warning", "error")
  
  for (i in 1:length(rho)) {
    
    rho1 <- rho[i]
    
    if (n.OM == 0) {
      ### Perfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of M ###
      command <- character(0)

      command.MS.g.start <- list()
      command.MS.g.final <- list()
      n.MS <- list()

      ### Extract N of each MS ###
      for (n in 1:length(MS.num)) {
        n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
      }

      n.MS.num <- 1:length(n.MS) # Get length of MS

      ### Determine the factor of the bifactor model is the mediator in the model ###
      if (mediator %>% as.numeric() == 0) {
        mediator.focus <- "MG"
      } else if (mediator %>% as.numeric() != 0) {
        for (n in 1:length(n.MS.num)) {
          if (mediator %>% as.numeric() == n.MS.num[n]) {
            mediator.focus <- paste0("MS", n.MS.num[n])
          }
        }
      }

      ### Create the start syntax of MS: MS1., MS2., MS3. ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g.start[[n]] <- paste0("MS", n.MS.num[n], ".1")
      }

      ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g <- command.MS.g.start[[n]]
        for (asd in 1:(n.MS[[n]] - 1)) {
          command.MS.g <- paste0(command.MS.g, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
        }
        command.MS.g.final[[n]] <- command.MS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of M ###
      command <- paste0("MG =~ NA * ", command.MS.g.final[[1]])
      for (n in 2:length(command.MS.g.final)) {
        command <- paste0(command, " + ", command.MS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    } else if (n.OM != 0) {
      ### Imperfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of M ###
      command <- character(0)

      command.MS.g.start <- list()
      command.MS.g.final <- list()
      n.MS <- list()

      ### Extract N of each MS ###
      for (n in 1:length(MS.num)) {
        n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
      }

      n.MS.num <- 1:length(n.MS) # Get length of MS

      ### Determine the factor of the bifactor model is the mediator in the model ###
      if (mediator %>% as.numeric() == 0) {
        mediator.focus <- "MG"
      } else if (mediator %>% as.numeric() != 0) {
        for (n in 1:length(n.MS.num)) {
          if (mediator %>% as.numeric() == n.MS.num[n]) {
            mediator.focus <- paste0("MS", n.MS.num[n])
          }
        }
      }

      ### Create the start syntax of MS: MS1., MS2., MS3. ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g.start[[n]] <- paste0("MS", n.MS.num[n], ".1")
      }

      ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
      for (n in 1:length(n.MS.num)) {
        command.MS.g <- command.MS.g.start[[n]]
        for (asd in 1:(n.MS[[n]] - 1)) {
          command.MS.g <- paste0(command.MS.g, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
        }
        command.MS.g.final[[n]] <- command.MS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of M ###
      command <- paste0("MG =~ NA * M.1")
      for (asd in 2:n.OM) {
        command <- paste0(command, " + ", "M.", asd)
      }

      for (n in 1:length(command.MS.g.final)) {
        command <- paste0(command, " + ", command.MS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    }

    ### Write syntax for the specific factors of M ###
    n.MS <- list()

    for (n in 1:length(MS.num)) {
    n.MS[[n]] <- sum(str.MS.4 == paste0("MS", MS.num[n]))
    }

    command.MS.start <- list()
    command.MS <- character(0)
    command.MS.final <- list()

    ### Create the start syntax of MS: MS11, MS21, MS31 ###
    for (n in 1:length(n.MS.num)) {
    command.MS.start[[n]] <- paste0(command.MS, "MS", n.MS.num[n], " =~ NA * MS", n.MS.num[n], ".1")
    }

    ### Create the rest syntax of MS (e.g., MS11, MS12, MS21, MS22, MS31, MS32 etc.) ###
    for (n in 1:length(n.MS.num)) {
    command.MS <- command.MS.start[[n]]
    for (asd in 1:(n.MS[[n]] - 1)) {
        command.MS <- paste0(command.MS, " + ", "MS", n.MS.num[[n]], ".", asd + 1)
    }
    command.MS.final[[n]] <- paste0(command.MS, "\n")
    }

    ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for specific factors of M ###
    command <- paste0(command, command.MS.final[[1]])
    for (n in 2:length(command.MS.final)) {
    command <- paste0(command, command.MS.final[[n]])
    }
    
    if (n.OY == 0) {
      ### Perfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of Y ###
      command.YS.g.start <- list()
      command.YS.g.final <- list()
      n.YS <- list()

      ### Extract N of each YS ###
      for (n in 1:length(YS.num)) {
      n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
      }

      n.YS.num <- 1:length(n.YS) # Get length of YS

      ### Determine the factor of the bifactor model is the response in the model ###
      if (response %>% as.numeric() == 0) {
      response.focus <- "YG"
      } else if (response %>% as.numeric() != 0) {
      for (n in 1:length(n.YS.num)) {
          if (response %>% as.numeric() == n.YS.num[n]) {
          response.focus <- paste0("YS", n.YS.num[n])
          }
      }
      }

      ### Create the start syntax of YS: YS1., YS2., YS3. ###
      for (n in 1:length(n.YS.num)) {
      command.YS.g.start[[n]] <- paste0("YS", n.YS.num[n], ".1")
      }

      ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
      for (n in 1:length(n.YS.num)) {
      command.YS.g <- command.YS.g.start[[n]]
      for (asd in 1:(n.YS[[n]] - 1)) {
          command.YS.g <- paste0(command.YS.g, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
      }
      command.YS.g.final[[n]] <- command.YS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of Y ###
      command <- paste0(command, "YG =~ NA * ", command.YS.g.final[[1]])
      for (n in 2:length(command.YS.g.final)) {
      command <- paste0(command, " + ", command.YS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    } else if (n.OY != 0) {
      ### Imperfect bifactor structure ###
      # Measurement model
      ### Write syntax for the general factor of Y ###
      command.YS.g.start <- list()
      command.YS.g.final <- list()
      n.YS <- list()

      ### Extract N of each YS ###
      for (n in 1:length(YS.num)) {
        n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
      }

      n.YS.num <- 1:length(n.YS) # Get length of YS

      ### Determine the factor of the bifactor model is the response in the model ###
      if (response %>% as.numeric() == 0) {
        response.focus <- "YG"
      } else if (response %>% as.numeric() != 0) {
        for (n in 1:length(n.YS.num)) {
          if (response %>% as.numeric() == n.YS.num[n]) {
            response.focus <- paste0("YS", n.YS.num[n])
          }
        }
      }

      ### Create the start syntax of YS: YS1., YS2., YS3. ###
      for (n in 1:length(n.YS.num)) {
        command.YS.g.start[[n]] <- paste0("YS", n.YS.num[n], ".1")
      }

      ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
      for (n in 1:length(n.YS.num)) {
        command.YS.g <- command.YS.g.start[[n]]
        for (asd in 1:(n.YS[[n]] - 1)) {
          command.YS.g <- paste0(command.YS.g, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
        }
        command.YS.g.final[[n]] <- command.YS.g
      }

      ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for the general factor of Y ###
      command <- paste0(command, "YG =~ NA * Y.1")
      for (asd in 2:n.OY) {
        command <- paste0(command, " + ", "Y.", asd)
      }
      for (n in 1:length(command.YS.g.final)) {
        command <- paste0(command, " + ", command.YS.g.final[[n]])
      }

      command <- paste0(command, "\n")
    }

    ### Write syntax for the specific factors of Y ###
    n.YS <- list()

    for (n in 1:length(YS.num)) {
    n.YS[[n]] <- sum(str.YS.4 == paste0("YS", YS.num[n]))
    }

    command.YS.start <- list()
    command.YS <- character(0)
    command.YS.final <- list()

    ### Create the start syntax of YS: YS11, YS21, YS31 ###
    for (n in 1:length(n.YS.num)) {
    command.YS.start[[n]] <- paste0(command.YS, "YS", n.YS.num[n], " =~ NA * YS", n.YS.num[n], ".1")
    }

    ### Create the rest syntax of YS (e.g., YS11, YS12, YS21, YS22, YS31, YS32 etc.) ###
    for (n in 1:length(n.YS.num)) {
    command.YS <- command.YS.start[[n]]
    for (asd in 1:(n.YS[[n]] - 1)) {
        command.YS <- paste0(command.YS, " + ", "YS", n.YS.num[[n]], ".", asd + 1)
    }
    command.YS.final[[n]] <- paste0(command.YS, "\n")
    }

    ### Combine the start sytanx and the rest of syntax to create the complete syntaxs for specific factors of Y ###
    command <- paste0(command, command.YS.final[[1]])
    for (n in 2:length(command.YS.final)) {
    command <- paste0(command, command.YS.final[[n]])
    }

    # Regression model
    if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    
    ### Replace mediator M according to researcher's choice ###
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY == 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }      
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY == 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM == 0 & n.CY != 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY == 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    
    command.m.new <- paste0(command.m.new, "\n")
    command.y.new <- paste0(command.y.new, "\n")
    
    } else if (n.X != 0 & n.Y != 0 & n.M != 0 & n.CM != 0 & n.CY != 0 & n.CMY != 0) {
    
    command.m <- character(0) ### Create a new character to save the syntax for M
    command.y <- character(0) ### Create a new character to save the syntax for Y
    command.m <- paste0(command.m, "M ~ alpha*X")
    command.y <- paste0(command.y, "Y ~ gamma*X + beta*M")
    
    ### Replace mediator M according to researcher's choice ###
    command.m.new <- gsub("M", mediator.focus, command.m)
    command.y.new.1 <- gsub("M", mediator.focus, command.y)
    
    ### Replace mediator Y according to researcher's choice ###
    command.y.new <- gsub("Y", response.focus, command.y.new.1)
    
    for (asd in 1:n.CM) {
        command.m.new <- paste0(command.m.new, " + ", "CM.", asd)
    }
    for (asd in 1:n.CMY) {
        command.m.new <- paste0(command.m.new, " + ", "CMY.", asd)
    }
    command.m.new <- paste0(command.m.new, "\n")
    
    for (asd in 1:n.CY) {
        command.y.new <- paste0(command.y.new, " + ", "CY.", asd)
    }
    for (asd in 1:n.CMY) {
        command.y.new <- paste0(command.y.new, " + ", "CMY.", asd)
    }
    command.y.new <- paste0(command.y.new, "\n")
    }

    # Residual correlations
    command.r <- character(0) ### Create a new character to save the syntax for residuals
    command.r <- paste0(command.r, "Y ~~ r*M")
    command.r <- paste0(command.r, "\n")
    command.r <- paste0(command.r, "M ~~ 1*M") # For correlation between M and Y
    command.r <- paste0(command.r, "\n") 
    command.r <- paste0(command.r, "Y ~~ 1*Y", "\n") # For correlation between M and Y

    ### Replace mediator M according to researcher's choice ###
    command.r.new.1 <- gsub("M", mediator.focus, command.r)
    command.r.new <- gsub("Y", response.focus, command.r.new.1)

    command <- paste0(command, command.m.new, command.y.new, command.r.new)

    ### Unit variance identification (UVI) ###
    if (mediator == 0) {
      ### Mediator is the general factor ###
      for (asd in 1:length(n.MS.num)) {
        command <- paste0(command, "MS", asd, " ~~ 1 * MS", asd, "\n")
      }
    } else if (mediator != 0) {
      ### Mediator is the specific factor ###
      command <- paste0(command, "MG ~~ 1 * MG", "\n")
      n.MS.num.new <- n.MS.num[n.MS.num != mediator]
      for (asd in 1:length(n.MS.num.new)) {
        command <- paste0(command, "MS", n.MS.num.new[asd], " ~~ 1 * MS", n.MS.num.new[asd], "\n")
      }
    }

    if (response == 0) {
      ### Response is the general factor ###
      for (asd in 1:length(n.YS.num)) {
        command <- paste0(command, "YS", asd, " ~~ 1 * YS", asd, "\n")
      }
    } else if (response != 0) {
      ### Response is the specific factor ###
      command <- paste0(command, "YG ~~ 1 * YG", "\n")
      n.YS.num.new <- n.YS.num[n.YS.num != response]
      for (asd in 1:length(n.YS.num.new)) {
        command <- paste0(command, "YS", n.YS.num.new[asd], " ~~ 1 * YS", n.YS.num.new[asd], "\n")
      }
    }

    ### Residual correlations of bifactor model ###
    ### Get correlations between the general factor and specific factors of M ###
    res.cor.MS <- function(x) {
    cor.MS <- paste0("MG ~~ 0 * MS", x)
    return(cor.MS)
    }

    cor.gs <- lapply(n.MS.num %>% as.list, res.cor.MS)

    command <- paste0(command, cor.gs[[1]], "\n")

    for (n in 2:length(cor.gs)) {
    command <- paste0(command, cor.gs[[n]], "\n")
    }

    ### Get correlations between sepcific factors ###

    com.ss <- combn(paste0("MS", 1:length(n.MS.num)), 2) # Get all combinations of specific factors

    cor.ss <- list()

    ### Get syntax of correlations between specific factors ###
    for (n in 1:ncol(com.ss)) {
    cor.ss[[n]] <- paste0(com.ss[1, n], " ~~ 0 * ", com.ss[2, n])
    }

    if (ncol(com.ss) == 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    } else if (ncol(com.ss) > 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    for (n in 2:length(cor.ss)) {
        command <- paste0(command, cor.ss[[n]], "\n")
    }
    }

    ### Residual correlations of bifactor model ###
    ### Get correlations between the general factor and specific factors of Y ###
    res.cor.YS <- function(x) {
    cor.YS <- paste0("YG ~~ 0 * YS", x)
    return(cor.YS)
    }

    cor.gs <- lapply(n.YS.num %>% as.list, res.cor.YS)

    command <- paste0(command, cor.gs[[1]], "\n")

    for (n in 2:length(cor.gs)) {
    command <- paste0(command, cor.gs[[n]], "\n")
    }

    ### Get correlations between sepcific factors ###

    com.ss <- combn(paste0("YS", 1:length(n.YS.num)), 2) # Get all combinations of specific factors

    cor.ss <- list()

    ### Get syntax of correlations between specific factors ###
    for (n in 1:ncol(com.ss)) {
    cor.ss[[n]] <- paste0(com.ss[1, n], " ~~ 0 * ", com.ss[2, n])
    }

    if (ncol(com.ss) == 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    } else if (ncol(com.ss) > 1) {
    
    command <- paste0(command, cor.ss[[1]], "\n")
    
    for (n in 2:length(cor.ss)) {
        command <- paste0(command, cor.ss[[n]], "\n")
    }
    }

    # Mediation effect
    command <- paste0(command, "ind := alpha*beta")

    # No direct effect: no c prime 
    # if (effect == 0) {
    #   command <- gsub("gamma", 0, command)
    # }
    
    command <- gsub("r",rho1,command)
    
    ### Check warnings and errors ###
    war <- myTryCatch(output <- sem(command, data = dataset, se = method, bootstrap = bootnum))
    
    ### Save warnings and error ###
    if (is.null(war$warnings)) {
      ### No warnings ###
      war.res[i, 1] <- 0
      
    } else if (length(war$warnings) == 1) {

      if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings)
          | grepl('bootstrap runs failed or did not converge.', war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
      
    } else if (length(war$warnings) > 1) {
      
      war.temp <- length(war$warnings)
      
      for (n in 1:length(war$warnings)) {
        
        if (grepl('bootstrap runs resulted in nonadmissible solutions.', war$warnings[n])
            | grepl('bootstrap runs failed or did not converge.', war$warnings[n])) {
          war.temp[n] <- 1
        } else {
          war.temp[n] <- 0
        }
      }
      
      if (sum(war.temp) == length(war$warnings)) {
        war.res[i, 1] <- 0
      } else {
        war.res[i, 1] <- 1
      }
    }
    
    if (is.null(war$error)) {
      war.res[i, 2] <- 0
    } else {
      war.res[i, 2] <- 1
    }
    
    if (war.res$warning[i] == 0 & war.res$error[i] == 0) {
      meffect.est[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$est
      meffect.upper[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.upper
      meffect.lower[i] <- tail(parameterEstimates(output, boot.ci.type = "perc", level = 0.95, ci = TRUE), 1)$ci.lower
    }
    
  }
  result <- data.frame(meffect.est, meffect.upper,meffect.lower,rho, war.res)
  return(result)
}
###########################################################

###### Draw plot ######
DrawPlot <- function(result) {
  
  # if (range[1] <= range[2]) {
    
  #   rho.range <- seq(range[1], range[2], .1)
    
  # } else if (range[1] > range[2]) {
    
  #   rho.range <- seq(range[2], range[1], .1)
    
  # }
  
  ### Extract result to analyze based on the rho range ###
  # result.new <- result[result$rho >= rho.range[1] & result$rho < (rho.range[length(rho.range)] + .1), ]
  
  ### If there is no warnings and all results can be used based on the rho range ###
  if (sum(result$warning) == 0 & sum(result$error == 0)) {
    ggplot(result, aes(x = rho, y = meffect.est)) +
    # geom_point() +
    geom_line() +
    # geom_smooth(method = "loess", colour = "black") +
    # Remove background 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    # Add a line that indicates estimated mediation effect
    geom_hline(yintercept = 0, color = "black") +
    geom_vline(xintercept = 0, color = "black") +
    xlab(TeX('$\\rho$')) +
    ylab(TeX('Mediation Effect Estimate: $\\alpha\\beta$')) +
    geom_ribbon(aes(ymin = meffect.lower, ymax = meffect.upper), linetype = 2, alpha = 0.1) +
    # scale_x_continuous(breaks = seq(-1.0, 1.0, by = 0.1)) +
    # Change font of xlabel and ylabel
    theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
    
  } else if (sum(result$warning) != 0 | sum(result$error) != 0) {
    
    result.new <- result.new[result$warning == 0 & result$error == 0, ]
    
    ggplot(result.new, aes(x = rho, y = meffect.est)) +
    # geom_point() +
    geom_line() +
    # geom_smooth(method = "loess", colour = "black") +
    # Remove background 
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    # Add a line that indicates estimated mediation effect
    geom_hline(yintercept = 0, color = "black") +
    geom_vline(xintercept = 0, color = "black") +
    xlab(TeX('$\\rho$')) +
    ylab(TeX('Mediation Effect Estimate: $\\alpha\\beta$')) +
    geom_ribbon(aes(ymin = meffect.lower, ymax = meffect.upper), linetype = 2, alpha = 0.1) +
    # scale_x_continuous(breaks = seq(-1.0, 1.0, by = 0.1)) +
    # Change font of xlabel and ylabel
    theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
  }
}


# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Sensitivity analysis"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("file1", "Choose File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Horizontal line ----
      tags$hr(),

      # Input: Select variables ----
      textInput("predictor", "Which columns are X?", ""),     
     
      ### Choose bifactor model or not ###
      selectInput(
        "modelM",
        "Is M a bifactor model?",
        c(` ` = "0", `No` = "No", `Yes` = "Yes")
      ),

      conditionalPanel(
        condition = "input.modelM == 'No'",
        textInput("mediator", "Which columns are indicators of M?", "")
      ),

      # conditionalPanel(
      #   condition = "input.modelM == 'Yes'",
      #   selectInput(
      #   "modelstructM",
      #   "Is the bifactor structure perfect? (Note: a perfect bifactor structure means all of items will load both on the general factor and one specific factor; an imperfect bifactor structure means some items will load only on the general factor, and some items will load both on the general factor and one specific factor)",
      #   c(` ` = "0", `No` = "No", `Yes` = "Yes")
      #   ),
      # ),

      conditionalPanel(
        condition = "input.modelM == 'Yes'",
        textInput("mediatorog", "Which columns are indicators that are exclusively loaded on the general factor of M? (Note: Please leave this box empty if there are no indicators that are exclusively loaded on the general factor)", "")
      ),

      conditionalPanel(
        condition = "input.modelM == 'Yes'",
        textInput("mediatornum", "How many specific facotors of M? (At least 2 specific factors)", "")
      ),

      conditionalPanel(
        condition = "input.modelM == 'Yes' && input.mediatornum == '2'",
        textInput("mediators1", "Which columns are indicators of the specific factor 1 of M?", ""),
        textInput("mediators2", "Which columns are indicators of the specific factor 2 of M?", "")
      ),

      conditionalPanel(
        condition = "input.modelM == 'Yes' && input.mediatornum == '3'",
        textInput("mediators1", "Which columns are indicators of the specific factor 1 of M?", ""),
        textInput("mediators2", "Which columns are indicators of the specific factor 2 of M?", ""),
        textInput("mediators3", "Which columns are indicators of the specific factor 3 of M?", "")
      ),

       conditionalPanel(
        condition = "input.modelM == 'Yes' && input.mediatornum == '4'",
        textInput("mediators1", "Which columns are indicators of the specific factor 1 of M?", ""),
        textInput("mediators2", "Which columns are indicators of the specific factor 2 of M?", ""),
        textInput("mediators3", "Which columns are indicators of the specific factor 3 of M?", ""),
        textInput("mediators4", "Which columns are indicators of the specific factor 4 of M?", ""),
      ),

      conditionalPanel(
        condition = "input.modelM == 'Yes' && input.mediatornum == '5'",
        textInput("mediators1", "Which columns are indicators of the specific factor 1 of M?", ""),
        textInput("mediators2", "Which columns are indicators of the specific factor 2 of M?", ""),
        textInput("mediators3", "Which columns are indicators of the specific factor 3 of M?", ""),
        textInput("mediators4", "Which columns are indicators of the specific factor 4 of M?", ""),
        textInput("mediators5", "Which columns are indicators of the specific factor 5 of M?", "")
      ),

      conditionalPanel(
        condition = "input.modelM == 'Yes' && input.mediatornum == '6'",
        textInput("mediators1", "Which columns are indicators of the specific factor 1 of M?", ""),
        textInput("mediators2", "Which columns are indicators of the specific factor 2 of M?", ""),
        textInput("mediators3", "Which columns are indicators of the specific factor 3 of M?", ""),
        textInput("mediators4", "Which columns are indicators of the specific factor 4 of M?", ""),
        textInput("mediators5", "Which columns are indicators of the specific factor 5 of M?", ""),
        textInput("mediators6", "Which columns are indicators of the specific factor 6 of M?", "")
      ),


      conditionalPanel(
        condition = "input.modelM == 'Yes' && input.mediatornum == '7'",
        textInput("mediators1", "Which columns are indicators of the specific factor 1 of M?", ""),
        textInput("mediators2", "Which columns are indicators of the specific factor 2 of M?", ""),
        textInput("mediators3", "Which columns are indicators of the specific factor 3 of M?", ""),
        textInput("mediators4", "Which columns are indicators of the specific factor 4 of M?", ""),
        textInput("mediators5", "Which columns are indicators of the specific factor 5 of M?", ""),
        textInput("mediators6", "Which columns are indicators of the specific factor 6 of M?", ""),
        textInput("mediators7", "Which columns are indicators of the specific factor 7 of M?", "")
      ),


      conditionalPanel(
        condition = "input.modelM == 'Yes' && input.mediatornum == '8'",
        textInput("mediators1", "Which columns are indicators of the specific factor 1 of M?", ""),
        textInput("mediators2", "Which columns are indicators of the specific factor 2 of M?", ""),
        textInput("mediators3", "Which columns are indicators of the specific factor 3 of M?", ""),
        textInput("mediators4", "Which columns are indicators of the specific factor 4 of M?", ""),
        textInput("mediators5", "Which columns are indicators of the specific factor 5 of M?", ""),
        textInput("mediators6", "Which columns are indicators of the specific factor 6 of M?", ""),
        textInput("mediators7", "Which columns are indicators of the specific factor 7 of M?", ""),
        textInput("mediators8", "Which columns are indicators of the specific factor 8 of M?", "")
      ),


      conditionalPanel(
        condition = "input.modelM == 'Yes' && input.mediatornum == '9'",
        textInput("mediators1", "Which columns are indicators of the specific factor 1 of M?", ""),
        textInput("mediators2", "Which columns are indicators of the specific factor 2 of M?", ""),
        textInput("mediators3", "Which columns are indicators of the specific factor 3 of M?", ""),
        textInput("mediators4", "Which columns are indicators of the specific factor 4 of M?", ""),
        textInput("mediators5", "Which columns are indicators of the specific factor 5 of M?", ""),
        textInput("mediators6", "Which columns are indicators of the specific factor 6 of M?", ""),
        textInput("mediators7", "Which columns are indicators of the specific factor 7 of M?", ""),
        textInput("mediators8", "Which columns are indicators of the specific factor 8 of M?", ""),
        textInput("mediators9", "Which columns are indicators of the specific factor 9 of M?", "")
      ),


      conditionalPanel(
        condition = "input.modelM == 'Yes' && input.mediatornum == '10'",
        textInput("mediators1", "Which columns are indicators of the specific factor 1 of M?", ""),
        textInput("mediators2", "Which columns are indicators of the specific factor 2 of M?", ""),
        textInput("mediators3", "Which columns are indicators of the specific factor 3 of M?", ""),
        textInput("mediators4", "Which columns are indicators of the specific factor 4 of M?", ""),
        textInput("mediators5", "Which columns are indicators of the specific factor 5 of M?", ""),
        textInput("mediators6", "Which columns are indicators of the specific factor 6 of M?", ""),
        textInput("mediators7", "Which columns are indicators of the specific factor 7 of M?", ""),
        textInput("mediators8", "Which columns are indicators of the specific factor 8 of M?", ""),
        textInput("mediators9", "Which columns are indicators of the specific factor 9 of M?", ""),
        textInput("mediators10", "Which columns are indicators of the specific factor 10 of M?", "")
      ),

      # ### Let applied researchers name factors of the bifactor model M ###
      #   conditionalPanel(
      #   condition = "input.modelM == 'Yes' && input.mediatornum == '2'",
      #   textInput("namemg", "Please input the name for the general factor", ""),
      #   textInput("namems1", "Please input the name for the specific factor 1", ""),
      #   textInput("namems2", "Please input the name for the specific factor 2", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelM == 'Yes' && input.mediatornum == '3'",
      #   textInput("namemg", "Please input the name for the general factor", ""),
      #   textInput("namems1", "Please input the name for the specific factor 1", ""),
      #   textInput("namems2", "Please input the name for the specific factor 2", ""),
      #   textInput("namems3", "Please input the name for the specific factor 3", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelM == 'Yes' && input.mediatornum == '4'",
      #   textInput("namemg", "Please input the name for the general factor", ""),
      #   textInput("namems1", "Please input the name for the specific factor 1", ""),
      #   textInput("namems2", "Please input the name for the specific factor 2", ""),
      #   textInput("namems3", "Please input the name for the specific factor 3", ""),
      #   textInput("namems4", "Please input the name for the specific factor 4", "")
      # ),

      #  conditionalPanel(
      #   condition = "input.modelM == 'Yes' && input.mediatornum == '5'",
      #   textInput("namemg", "Please input the name for the general factor", ""),
      #   textInput("namems1", "Please input the name for the specific factor 1", ""),
      #   textInput("namems2", "Please input the name for the specific factor 2", ""),
      #   textInput("namems3", "Please input the name for the specific factor 3", ""),
      #   textInput("namems4", "Please input the name for the specific factor 4", ""),
      #   textInput("namems5", "Please input the name for the specific factor 5", "")
      # ),

      #  conditionalPanel(
      #   condition = "input.modelM == 'Yes' && input.mediatornum == '6'",
      #   textInput("namemg", "Please input the name for the general factor", ""),
      #   textInput("namems1", "Please input the name for the specific factor 1", ""),
      #   textInput("namems2", "Please input the name for the specific factor 2", ""),
      #   textInput("namems3", "Please input the name for the specific factor 3", ""),
      #   textInput("namems4", "Please input the name for the specific factor 4", ""),
      #   textInput("namems5", "Please input the name for the specific factor 5", ""),
      #   textInput("namems6", "Please input the name for the specific factor 6", "")
      # ),

      #  conditionalPanel(
      #   condition = "input.modelM == 'Yes' && input.mediatornum == '7'",
      #   textInput("namemg", "Please input the name for the general factor", ""),
      #   textInput("namems1", "Please input the name for the specific factor 1", ""),
      #   textInput("namems2", "Please input the name for the specific factor 2", ""),
      #   textInput("namems3", "Please input the name for the specific factor 3", ""),
      #   textInput("namems4", "Please input the name for the specific factor 4", ""),
      #   textInput("namems5", "Please input the name for the specific factor 5", ""),
      #   textInput("namems6", "Please input the name for the specific factor 6", ""),
      #   textInput("namems7", "Please input the name for the specific factor 7", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelM == 'Yes' && input.mediatornum == '8'",
      #   textInput("namemg", "Please input the name for the general factor?", ""),
      #   textInput("namems1", "Please input the name for the specific factor 1", ""),
      #   textInput("namems2", "Please input the name for the specific factor 2", ""),
      #   textInput("namems3", "Please input the name for the specific factor 3", ""),
      #   textInput("namems4", "Please input the name for the specific factor 4", ""),
      #   textInput("namems5", "Please input the name for the specific factor 5", ""),
      #   textInput("namems6", "Please input the name for the specific factor 6", ""),
      #   textInput("namems7", "Please input the name for the specific factor 7", ""),
      #   textInput("namems8", "Please input the name for the specific factor 8", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelM == 'Yes' && input.mediatornum == '9'",
      #   textInput("namemg", "Please input the name for the general factor", ""),
      #   textInput("namems1", "Please input the name for the specific factor 1", ""),
      #   textInput("namems2", "Please input the name for the specific factor 2", ""),
      #   textInput("namems3", "Please input the name for the specific factor 3", ""),
      #   textInput("namems4", "Please input the name for the specific factor 4", ""),
      #   textInput("namems5", "Please input the name for the specific factor 5", ""),
      #   textInput("namems6", "Please input the name for the specific factor 6", ""),
      #   textInput("namems7", "Please input the name for the specific factor 7", ""),
      #   textInput("namems8", "Please input the name for the specific factor 8", ""),
      #   textInput("namems9", "Please input the name for the specific factor 9", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelM == 'Yes' && input.mediatornum == '10'",
      #   textInput("namemg", "Please input the name for the general factor", ""),
      #   textInput("namems1", "Please input the name for the specific factor 1", ""),
      #   textInput("namems2", "Please input the name for the specific factor 2", ""),
      #   textInput("namems3", "Please input the name for the specific factor 3", ""),
      #   textInput("namems4", "Please input the name for the specific factor 4", ""),
      #   textInput("namems5", "Please input the name for the specific factor 5", ""),
      #   textInput("namems6", "Please input the name for the specific factor 6", ""),
      #   textInput("namems7", "Please input the name for the specific factor 7", ""),
      #   textInput("namems8", "Please input the name for the specific factor 8", ""),
      #   textInput("namems9", "Please input the name for the specific factor 9", ""),
      #   textInput("namems10", "Please input the name for the specific factor 10", "")
      # ),

      ### Let applied researchers choose the mediator of the bifactor model ###
      conditionalPanel(
        condition = "input.modelM == 'Yes' && input.mediatornum > 1",
        textInput(
        "choicem",
        "Which factor of the bifactor model is the mediator in the model? (Note: Please input numbers. 0 indicates the general factor; 1 indicates the specific factor 1; 2 indicates the specific factor 2; 3 indicates the specific factor 3; etc.)", 
        "")
      ),


      selectInput(
        "modelY",
        "Is Y a bifactor model?",
        c(` ` = "0", `No` = "No", `Yes` = "Yes")
      ),

      conditionalPanel(
        condition = "input.modelY == 'No'",
        textInput("response", "Which columns are indicators of Y?", "")
      ),

      # conditionalPanel(
      #   condition = "input.modelY == 'Yes'",
      #   selectInput(
      #   "modelstructY",
      #   "Is the bifactor structure perfect? (Note: a perfect bifactor structure means all of items will load both on the general factor and one specific factor; an imperfect bifactor structure means some items will load only on the general factor, and some items will load both on the general factor and one specific factor)",
      #   c(` ` = "0", `No` = "No", `Yes` = "Yes")
      #   ),
      # ),

      conditionalPanel(
        condition = "input.modelY == 'Yes'",
        textInput("responseog", "Which columns are indicators that are exclusively loaded on the general factor of Y? (Note: Please leave this box empty if there are no indicators that are exclusively loaded on the general factor)", "")
      ),

      conditionalPanel(
        condition = "input.modelY == 'Yes'",
        textInput("responsenum", "How many specific facotors of Y? (At least 2 specific factors)", "")
      ),

      conditionalPanel(
        condition = "input.modelY == 'Yes' && input.responsenum == '2'",
        textInput("responses1", "Which columns are indicators of the specific factor 1 of Y?", ""),
        textInput("responses2", "Which columns are indicators of the specific factor 2 of Y?", "")
      ),


      conditionalPanel(
        condition = "input.modelY == 'Yes' && input.responsenum == '3'",
        textInput("responses1", "Which columns are indicators of the specific factor 1 of Y?", ""),
        textInput("responses2", "Which columns are indicators of the specific factor 2 of Y?", ""),
        textInput("responses3", "Which columns are indicators of the specific factor 3 of Y?", "")
      ),

      conditionalPanel(
        condition = "input.modelY == 'Yes' && input.responsenum == '4'",
        textInput("responses1", "Which columns are indicators of the specific factor 1 of Y?", ""),
        textInput("responses2", "Which columns are indicators of the specific factor 2 of Y?", ""),
        textInput("responses3", "Which columns are indicators of the specific factor 3 of Y?", ""),
        textInput("responses4", "Which columns are indicators of the specific factor 4 of Y?", "")
      ),


      conditionalPanel(
        condition = "input.modelY == 'Yes' && input.responsenum == '5'",
        textInput("responses1", "Which columns are indicators of the specific factor 1 of Y?", ""),
        textInput("responses2", "Which columns are indicators of the specific factor 2 of Y?", ""),
        textInput("responses3", "Which columns are indicators of the specific factor 3 of Y?", ""),
        textInput("responses4", "Which columns are indicators of the specific factor 4 of Y?", ""),
        textInput("responses5", "Which columns are indicators of the specific factor 5 of Y?", "")
      ),

      conditionalPanel(
        condition = "input.modelY == 'Yes' && input.responsenum == '6'",
        textInput("responses1", "Which columns are indicators of the specific factor 1 of Y?", ""),
        textInput("responses2", "Which columns are indicators of the specific factor 2 of Y?", ""),
        textInput("responses3", "Which columns are indicators of the specific factor 3 of Y?", ""),
        textInput("responses4", "Which columns are indicators of the specific factor 4 of Y?", ""),
        textInput("responses5", "Which columns are indicators of the specific factor 5 of Y?", ""),
        textInput("responses6", "Which columns are indicators of the specific factor 6 of Y?", "")
      ),

      conditionalPanel(
        condition = "input.modelY == 'Yes' && input.responsenum == '7'",
        textInput("responses1", "Which columns are indicators of the specific factor 1 of Y?", ""),
        textInput("responses2", "Which columns are indicators of the specific factor 2 of Y?", ""),
        textInput("responses3", "Which columns are indicators of the specific factor 3 of Y?", ""),
        textInput("responses4", "Which columns are indicators of the specific factor 4 of Y?", ""),
        textInput("responses5", "Which columns are indicators of the specific factor 5 of Y?", ""),
        textInput("responses6", "Which columns are indicators of the specific factor 6 of Y?", ""),
        textInput("responses7", "Which columns are indicators of the specific factor 7 of Y?", "")
      ),

       conditionalPanel(
        condition = "input.modelY == 'Yes' && input.responsenum == '8'",
        textInput("responses1", "Which columns are indicators of the specific factor 1 of Y?", ""),
        textInput("responses2", "Which columns are indicators of the specific factor 2 of Y?", ""),
        textInput("responses3", "Which columns are indicators of the specific factor 3 of Y?", ""),
        textInput("responses4", "Which columns are indicators of the specific factor 4 of Y?", ""),
        textInput("responses5", "Which columns are indicators of the specific factor 5 of Y?", ""),
        textInput("responses6", "Which columns are indicators of the specific factor 6 of Y?", ""),
        textInput("responses7", "Which columns are indicators of the specific factor 7 of Y?", ""),
        textInput("responses8", "Which columns are indicators of the specific factor 8 of Y?", "")
      ),

      conditionalPanel(
        condition = "input.modelY == 'Yes' && input.responsenum == '9'",
        textInput("responses1", "Which columns are indicators of the specific factor 1 of Y?", ""),
        textInput("responses2", "Which columns are indicators of the specific factor 2 of Y?", ""),
        textInput("responses3", "Which columns are indicators of the specific factor 3 of Y?", ""),
        textInput("responses4", "Which columns are indicators of the specific factor 4 of Y?", ""),
        textInput("responses5", "Which columns are indicators of the specific factor 5 of Y?", ""),
        textInput("responses6", "Which columns are indicators of the specific factor 6 of Y?", ""),
        textInput("responses7", "Which columns are indicators of the specific factor 7 of Y?", ""),
        textInput("responses8", "Which columns are indicators of the specific factor 8 of Y?", ""),
        textInput("responses9", "Which columns are indicators of the specific factor 9 of Y?", "")
      ),

      conditionalPanel(
        condition = "input.modelY == 'Yes' && input.responsenum == '10'",
        textInput("responses1", "Which columns are indicators of the specific factor 1 of Y?", ""),
        textInput("responses2", "Which columns are indicators of the specific factor 2 of Y?", ""),
        textInput("responses3", "Which columns are indicators of the specific factor 3 of Y?", ""),
        textInput("responses4", "Which columns are indicators of the specific factor 4 of Y?", ""),
        textInput("responses5", "Which columns are indicators of the specific factor 5 of Y?", ""),
        textInput("responses6", "Which columns are indicators of the specific factor 6 of Y?", ""),
        textInput("responses7", "Which columns are indicators of the specific factor 7 of Y?", ""),
        textInput("responses8", "Which columns are indicators of the specific factor 8 of Y?", ""),
        textInput("responses9", "Which columns are indicators of the specific factor 9 of Y?", ""),
        textInput("responses10", "Which columns are indicators of the specific factor 10 of Y?", "")
      ),

      # ### Let applied researchers name factors of the bifactor model Y ###
      # conditionalPanel(
      #   condition = "input.modelY == 'Yes' && input.responsenum == '2'",
      #   textInput("nameyg", "Please input the name for the general factor", ""),
      #   textInput("nameys1", "Please input the name for the specific factor 1", ""),
      #   textInput("nameys2", "Please input the name for the specific factor 2", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelY == 'Yes' && input.responsenum == '3'",
      #   textInput("nameyg", "Please input the name for the general factor", ""),
      #   textInput("nameys1", "Please input the name for the specific factor 1", ""),
      #   textInput("nameys2", "Please input the name for the specific factor 2", ""),
      #   textInput("nameys3", "Please input the name for the specific factor 3", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelY == 'Yes' && input.responsenum == '4'",
      #   textInput("nameyg", "Please input the name for the general factor", ""),
      #   textInput("nameys1", "Please input the name for the specific factor 1", ""),
      #   textInput("nameys2", "Please input the name for the specific factor 2", ""),
      #   textInput("nameys3", "Please input the name for the specific factor 3", ""),
      #   textInput("nameys4", "Please input the name for the specific factor 4", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelY == 'Yes' && input.responsenum == '5'",
      #   textInput("nameyg", "Please input the name for the general factor", ""),
      #   textInput("nameys1", "Please input the name for the specific factor 1", ""),
      #   textInput("nameys2", "Please input the name for the specific factor 2", ""),
      #   textInput("nameys3", "Please input the name for the specific factor 3", ""),
      #   textInput("nameys4", "Please input the name for the specific factor 4", ""),
      #   textInput("nameys5", "Please input the name for the specific factor 5", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelY == 'Yes' && input.responsenum == '6'",
      #   textInput("nameyg", "Please input the name for the general factor", ""),
      #   textInput("nameys1", "Please input the name for the specific factor 1", ""),
      #   textInput("nameys2", "Please input the name for the specific factor 2", ""),
      #   textInput("nameys3", "Please input the name for the specific factor 3", ""),
      #   textInput("nameys4", "Please input the name for the specific factor 4", ""),
      #   textInput("nameys5", "Please input the name for the specific factor 5", ""),
      #   textInput("nameys6", "Please input the name for the specific factor 6", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelY == 'Yes' && input.responsenum == '7'",
      #   textInput("nameyg", "Please input the name for the general factor", ""),
      #   textInput("nameys1", "Please input the name for the specific factor 1", ""),
      #   textInput("nameys2", "Please input the name for the specific factor 2", ""),
      #   textInput("nameys3", "Please input the name for the specific factor 3", ""),
      #   textInput("nameys4", "Please input the name for the specific factor 4", ""),
      #   textInput("nameys5", "Please input the name for the specific factor 5", ""),
      #   textInput("nameys6", "Please input the name for the specific factor 6", ""),
      #   textInput("nameys7", "Please input the name for the specific factor 7", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelY == 'Yes' && input.responsenum == '8'",
      #   textInput("nameyg", "Please input the name for the general factor", ""),
      #   textInput("nameys1", "Please input the name for the specific factor 1", ""),
      #   textInput("nameys2", "Please input the name for the specific factor 2", ""),
      #   textInput("nameys3", "Please input the name for the specific factor 3", ""),
      #   textInput("nameys4", "Please input the name for the specific factor 4", ""),
      #   textInput("nameys5", "Please input the name for the specific factor 5", ""),
      #   textInput("nameys6", "Please input the name for the specific factor 6", ""),
      #   textInput("nameys7", "Please input the name for the specific factor 7", ""),
      #   textInput("nameys8", "Please input the name for the specific factor 8", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelY == 'Yes' && input.responsenum == '9'",
      #   textInput("nameyg", "Please input the name for the general factor", ""),
      #   textInput("nameys1", "Please input the name for the specific factor 1", ""),
      #   textInput("nameys2", "Please input the name for the specific factor 2", ""),
      #   textInput("nameys3", "Please input the name for the specific factor 3", ""),
      #   textInput("nameys4", "Please input the name for the specific factor 4", ""),
      #   textInput("nameys5", "Please input the name for the specific factor 5", ""),
      #   textInput("nameys6", "Please input the name for the specific factor 6", ""),
      #   textInput("nameys7", "Please input the name for the specific factor 7", ""),
      #   textInput("nameys8", "Please input the name for the specific factor 8", ""),
      #   textInput("nameys9", "Please input the name for the specific factor 9", "")
      # ),

      # conditionalPanel(
      #   condition = "input.modelY == 'Yes' && input.responsenum == '10'",
      #   textInput("nameyg", "Please input the name for the general factor", ""),
      #   textInput("nameys1", "Please input the name for the specific factor 1", ""),
      #   textInput("nameys2", "Please input the name for the specific factor 2", ""),
      #   textInput("nameys3", "Please input the name for the specific factor 3", ""),
      #   textInput("nameys4", "Please input the name for the specific factor 4", ""),
      #   textInput("nameys5", "Please input the name for the specific factor 5", ""),
      #   textInput("nameys6", "Please input the name for the specific factor 6", ""),
      #   textInput("nameys7", "Please input the name for the specific factor 7", ""),
      #   textInput("nameys8", "Please input the name for the specific factor 8", ""),
      #   textInput("nameys9", "Please input the name for the specific factor 9", ""),
      #   textInput("nameys10", "Please input the name for the specific factor 10", "")
      # ),


      ### Let applied researchers choose the response of the bifactor model ###
      conditionalPanel(
        condition = "input.modelY == 'Yes' && input.responsenum > 1",
        textInput(
        "choicey",
        "Which factor of the bifactor model is the dependent variable in the model? (Note: Please input numbers. 0 indicates the general factor; 1 indicates the specific factor 1; 2 indicates the specific factor 2; 3 indicates the specific factor 3; etc.)", 
        "")
      ),

      tags$hr(),
      
      # Input: Checkbox if file has header ----
      #checkboxInput("header", "Header", TRUE),
      
      
      # Inout: Select covariates ---
      selectInput(
        "covariate",
        "Do you have covariates?",
        c(`No` = "0", `Yes` = "Yes")
      ),
      # Only show this panel if covariate == Yes
      conditionalPanel(
        condition = "input.covariate == 'Yes'",
        textInput("conM", "Choose your covariates for M (Please leave this box blank if you do not have any covariates for M)", ""),
        textInput("conY", "Choose your covariates for Y (Please leave this box blank if you do not have any covariates for Y)", ""),
        textInput("conMY", "Choose your covariates for both M & Y (Please leave this box blank if you do not have any covariates for both M & Y)", "")
      ),
      
      # Horizontal line ----
      # tags$hr(),

      # selectInput(
      #   "mediationeffect",
      #   "Is this complete mediation or partial mediation? (Note: complete mediation indicates that the direct effect is zero; partial mediation indicates that the direct effect is not zero. 'Yes' means complete mediation; 'No' means partial mediation.)",
      #   c(`No` = "1", `Yes` = "0")
      # ),

      # tags$hr(),

      # numericRangeInput(
      # inputId = "corrange", label = "Correlation range input:",
      # value = c(-1, 1)
      # ),

      tags$hr(),

      selectInput(
        "method",
        "Which is your estimation method? (ML = Maximum Likelihood Estimation; MLR = Robust Maximum Likelihood Estimation.)",
        c(`ML` = "ML", `MLR` = "MLR", `Bootstrap` = "Bootstrap")
      ),

      ### Let applied researchers to choose number of bootstrap by themselves ###
      conditionalPanel(
        condition = "input.method == 'Bootstrap'",
        textInput(
        "bootstrapnum",
        "Input the number of bootstrap replicates", 
        "")
      ),

      # Horizontal line ----
      tags$hr(),
      
      
      actionButton("analysis","Analyze"),
      
      # Download button
      tags$hr(),
      
      downloadButton('export')
      
      # Input: Select number of rows to display ----
      #radioButtons("disp", "Display",
      #choices = c(Head = "head",
      #All = "all"),
      #selected = "head")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      tableOutput("contents"),
      tableOutput("updated"),
      h6("Note: CM = covariates for M; CY = covariates for Y; CMY = covariates for both M and Y; MS = specific factors of M; YS = specific factors of Y"),
      plotOutput("plot")
      # verbatimTextOutput("print1"),
      # verbatimTextOutput("print2"),
      # verbatimTextOutput("print3"),
      # verbatimTextOutput("print4")
      
    )
    
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    
    df <- fread(input$file1$datapath)
    
    return(head(df))
    
  })
  
  output$updated <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    
    df.new <- fread(input$file1$datapath)
    
    index.x <- as.numeric(input$predictor)
    names(df.new)[index.x] <- "X"
    
    ######## Read textinput string from the front interface #######
    textinp <- c(input$mediator, input$mediatorog, input$mediators1, input$mediators2, input$mediators3, input$mediators4, input$mediators5, input$mediators6, input$mediators7, input$mediators8, input$mediators9, 
    input$mediators10, input$response, input$responseog, input$responses1, input$responses2, input$responses3, input$responses4, input$responses5, input$responses6, input$responses7, input$responses8, 
    input$responses9, input$responses10, input$conM, input$conY, input$conMY) %>% as.list() ### Make textinput string as a list
    
    ### Create a list to save variable names ###
    var.name <- c("M", "M", paste0("M", paste0("S", 1:10)), "Y", "Y", paste0("Y", paste0("S", 1:10)), paste0("C", c("M", "Y", "MY"))) %>% as.list()

    ###### A function which can read numbers only and numbers with hyphen ######
    read.num <- function(input.var) {
      
      index <- input.var %>% str_split(., pattern = ";") %>% unlist() %>% as.numeric() ### Transfer textinput string to numeric
      if (sum(is.na(index)) == 0) {
      ### If applied researchers input numbers without hyphen in the textinput box ###
        index <- index[!is.na(index)]
        # names(df.new)[index] <- paste0("M", 1:length(index))
      } else if (sum(is.na(index)) != 0) {
      ### If applied researchers input numbers with hyphen in the textinput box ###
      
      ### Get all numbers without hyphen ###
        number <- input.var %>% str_split(., pattern = ";") %>% unlist() %>% as.numeric()
        number <- number[!is.na(number)]

      ### Numbers with hyphen ###
        index.str <- input.var %>% str_split(., pattern = ";") %>% unlist()
        index.hyphen <- strsplit(index.str, "(?=[-])", perl = TRUE)
        index.leng <- lapply(index.hyphen, length) %>% unlist ### Get length for list elements
        hyphen.num <- list() ### Create a list to save numbers which are indicated by the hyphen
        hyphen.num.2 <- list() ### Create a list to save numbers before and after the hyphen
        for (n in 1:length(index.hyphen)) {
          if (index.leng[[n]] == 3) {
            ### Calculate the count of numbers that hyphen indicates ### 
            diff <- tail(index.hyphen[[n]], 1) %>% as.numeric() - head(index.hyphen[[n]], 1) %>% as.numeric()
            ### Save numbers before and after the hyphen ###
            hyphen.ab <- c(head(index.hyphen[[n]], 1) %>% as.numeric(), tail(index.hyphen[[n]], 1) %>% as.numeric())
            ### Create a vector to save numbers for hyphen ###
            hyphen.diff <- numeric(diff - 1)
            ### Create numbers for the hyphen ###
            for (i in 1:(diff - 1)) {
              hyphen.diff[i] <- head(index.hyphen[[n]], 1) %>% as.numeric() + i
            }
            hyphen.num.2[[n]] <- hyphen.ab
            hyphen.num[[n]] <- hyphen.diff
          }
          ### Combine all numbers and remove duplicates and sort numbers by descreasing order ###
          index <- c(number, hyphen.num %>% unlist, hyphen.num.2 %>% unlist) %>% unique() %>% sort(., decreasing = FALSE)
          # names(df.new)[index] <- paste0(var.name, 1:length(index))
        }
      }
      return(index)
    }

    ### Save all variable numbers in a list ###
    index.num <- lapply(textinp, read.num)

    ### Give variable names according to their roles ###
    for (i in 1:length(index.num)) {
      names(df.new)[index.num[[i]]] <- paste0(var.name[[i]], ".", 1:length(index.num[[i]]))
    }

    return(head(df.new))
    
  })

  ### Draw plots ###  
  plot.sen <- eventReactive(input$analysis, {
    
   ### Refine variables for the plot which is the same the table ###

    df.new <- fread(input$file1$datapath)
    
    index.x <- as.numeric(input$predictor)
    names(df.new)[index.x] <- "X"
    
    ######## Read textinput string from the front interface #######
    textinp <- c(input$mediator, input$mediatorog, input$mediators1, input$mediators2, input$mediators3, input$mediators4, input$mediators5, input$mediators6, input$mediators7, input$mediators8, input$mediators9, 
    input$mediators10, input$response, input$responseog, input$responses1, input$responses2, input$responses3, input$responses4, input$responses5, input$responses6, input$responses7, input$responses8, 
    input$responses9, input$responses10, input$conM, input$conY, input$conMY) 
    
    ### Create a list to save variable names ###
    var.name <- c("M", "M", paste0("M", paste0("S", 1:10)), "Y", "Y", paste0("Y", paste0("S", 1:10)), paste0("C", c("M", "Y", "MY"))) 
    ###### A function which can read numbers only and numbers with hyphen ######
    read.num <- function(input.var) {
      
      index <- input.var %>% str_split(., pattern = ";") %>% unlist() %>% as.numeric() ### Transfer textinput string to numeric
      if (sum(is.na(index)) == 0) {
      ### If applied researchers input numbers without hyphen in the textinput box ###
        index <- index[!is.na(index)]
        # names(df.new)[index] <- paste0("M", 1:length(index))
      } else if (sum(is.na(index)) != 0) {
      ### If applied researchers input numbers with hyphen in the textinput box ###
      
      ### Get all numbers without hyphen ###
        number <- input.var %>% str_split(., pattern = ";") %>% unlist() %>% as.numeric()
        number <- number[!is.na(number)]

      ### Numbers with hyphen ###
        index.str <- input.var %>% str_split(., pattern = ";") %>% unlist()
        index.hyphen <- strsplit(index.str, "(?=[-])", perl = TRUE)
        index.leng <- lapply(index.hyphen, length) %>% unlist ### Get length for list elements
        hyphen.num <- list() ### Create a list to save numbers which are indicated by the hyphen
        hyphen.num.2 <- list() ### Create a list to save numbers before and after the hyphen
        for (n in 1:length(index.hyphen)) {
          if (index.leng[[n]] == 3) {
            ### Calculate the count of numbers that hyphen indicates ### 
            diff <- tail(index.hyphen[[n]], 1) %>% as.numeric() - head(index.hyphen[[n]], 1) %>% as.numeric()
            ### Save numbers before and after the hyphen ###
            hyphen.ab <- c(head(index.hyphen[[n]], 1) %>% as.numeric(), tail(index.hyphen[[n]], 1) %>% as.numeric())
            ### Create a vector to save numbers for hyphen ###
            hyphen.diff <- numeric(diff - 1)
            ### Create numbers for the hyphen ###
            for (i in 1:(diff - 1)) {
              hyphen.diff[i] <- head(index.hyphen[[n]], 1) %>% as.numeric() + i
            }
            hyphen.num.2[[n]] <- hyphen.ab
            hyphen.num[[n]] <- hyphen.diff
          }
          ### Combine all numbers and remove duplicates and sort numbers by descreasing order ###
          index <- c(number, hyphen.num %>% unlist, hyphen.num.2 %>% unlist) %>% unique() %>% sort(., decreasing = FALSE)
          # names(df.new)[index] <- paste0(var.name, 1:length(index))
        }
      }
      return(index)
    }

    ### Save all variable numbers in a list ###
    index.num <- lapply(textinp, read.num)

    ### Give variable names according to their roles ###
    for (i in 1:length(index.num)) {
      names(df.new)[index.num[[i]]] <- paste0(var.name[[i]], ".", 1:length(index.num[[i]]))
    }

    ### Create a list to save variable names for the index.num list ###
    var.name.2 <- c("m", "m1", paste0("m", paste0("S", 1:10)), "y", "y1", paste0("y", paste0("s", 1:10)), paste0("c", c("m", "y", "my"))) 
    ### Name the list: index.num ###
    names(index.num) <- var.name.2
    
    index.m.num <- input$mediatornum %>% as.numeric
    index.y.num <- input$responsenum %>% as.numeric
    index.method <- input$method
    index.modelM <- input$modelM
    index.modelY <- input$modelY
    index.m.ch <- input$choicem
    index.y.ch <- input$choicey
    index.boot.num <- input$bootstrapnum %>% as.numeric
    cor.range <- input$corrange
    # index.effect <- input$mediationeffect %>% as.numeric
    
    if (index.method != "Bootstrap") {

        if (index.modelM == "No" & index.modelY == "No") {

            if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cm) == 0 & length(index.num$cy) == 0 & length(index.num$cmy) == 0) {
            result <- Analyze.XMY(df.new, index.method)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cm) != 0 & length(index.num$cy) == 0 & length(index.num$cmy) == 0) {
            result <- Analyze.XMYCM(df.new, index.method)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cy) != 0 & length(index.num$cm) == 0 & length(index.num$cmy) == 0) {
            result <- Analyze.XMYCY(df.new, index.method)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cmy) != 0 & length(index.num$cm) == 0 & length(index.num$cy) == 0) {
            result <- Analyze.XMYCMY(df.new, index.method)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cmy) != 0 & length(index.num$cm) != 0 & length(index.num$cy) == 0) {
            result <- Analyze.XMYCMMY(df.new, index.method)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cmy) != 0 & length(index.num$cm) == 0 & length(index.num$cy) != 0) {
            result <- Analyze.XMYCYMY(df.new, index.method)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cmy) == 0 & length(index.num$cm) != 0 & length(index.num$cy) != 0) {
            result <- Analyze.XMYCMCY(df.new, index.method)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cmy) != 0 & length(index.num$cm) != 0 & length(index.num$cy) != 0) {
            result <- Analyze.all(df.new, index.method)
            DrawPlot(result)
            }

        } else if (index.modelM == "Yes" & index.modelY == "No") {

          result <- Analyze.bi.M(df.new, index.m.ch, index.method)
          DrawPlot(result)

        } else if (index.modelM == "No" & index.modelY == "Yes") {
          
          result <- Analyze.bi.Y(df.new, index.y.ch, index.method)
          DrawPlot(result)

        } else if (index.modelM == "Yes" & index.modelY == "Yes") {

          result <- Analyze.bi.MY(df.new, index.m.ch, index.y.ch, index.method)
          DrawPlot(result)

        }

    } else if (index.method == "Bootstrap") {

        if (index.modelM == "No" & index.modelY == "No") {

            if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cm) == 0 & length(index.num$cy) == 0 & length(index.num$cmy) == 0) {
            result <- Analyze.XMY.boot(df.new, index.boot.num)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cm) != 0 & length(index.num$cy) == 0 & length(index.num$cmy) == 0) {
            result <- Analyze.XMYCM.boot(df.new, index.boot.num)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cy) != 0 & length(index.num$cm) == 0 & length(index.num$cmy) == 0) {
            result <- Analyze.XMYCY.boot(df.new, index.boot.num)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cmy) != 0 & length(index.num$cm) == 0 & length(index.num$cy) == 0) {
            result <- Analyze.XMYCMY.boot(df.new, index.boot.num)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cmy) != 0 & length(index.num$cm) != 0 & length(index.num$cy) == 0) {
            result <- Analyze.XMYCMMY.boot(df.new, index.boot.num)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(iindex.num$cmy) != 0 & length(index.num$cm) == 0 & length(index.num$cy) != 0) {
            result <- Analyze.XMYCYMY.boot(df.new, index.boot.num)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cmy) == 0 & length(index.num$cm) != 0 & length(index.num$cy) != 0) {
            result <- Analyze.XMYCMCY.boot(df.new, index.boot.num)
            DrawPlot(result)
            } else if (length(index.x) != 0 & length(index.num$y) != 0 & length(index.num$m) != 0 & length(index.num$cmy) != 0 & length(index.num$cm) != 0 & length(index.num$cy) != 0) {
            result <- Analyze.all.boot(df.new, index.boot.num)
            DrawPlot(result)
            }
  
        } else if (index.modelM == "Yes" & index.modelY == "No") {

          result <- Analyze.bi.M.boot(df.new, index.m.ch, index.boot.num)
          DrawPlot(result)

        } else if (index.modelM == "No" & index.modelY == "Yes") {

          result <- Analyze.bi.Y.boot(df.new, index.y.ch, index.boot.num)
          DrawPlot(result)

        } else if (index.modelM == "Yes" & index.modelY == "Yes") {

          result <- Analyze.bi.MY.boot(df.new, index.m.ch, index.y.ch, index.boot.num)
          DrawPlot(result)

        }
    }
    
    
  })


  output$plot <- renderPlot({
    
    plot.sen()
    
  })

  # output$print1 <- renderPrint({ 
    
  #   index.method <- input$method
  #   index.method
    
  # })
  
  # output$print2 <- renderPrint({ 
    
  #   index.modelM <- input$modelM
  #   index.modelM
    
  # })
  
  # output$print3 <- renderPrint({ 
    
  #   index.modelY <- input$modelY
  #   index.modelY
    
  # })

  # output$print4 <- renderPrint({ 
    
  #   input$mediator
    
  # })

  
  # Download pdf file with A4 format
  output$export <- downloadHandler(
    filename = function() {"plot.pdf"},
    content = function(file) {
      ggsave(file, plot = plot.sen(), device = "pdf", width = 297, height = 210,, units = "mm")
    }
  )
}
# Run the app ----
shinyApp(ui, server)