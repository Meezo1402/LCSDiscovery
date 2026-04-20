# ==============================================================================
# 01_HMM_fitting.R
#
# Fits a K-state Hidden Markov Model with custom negative binomial emission
# distributions to longitudinal psychological symptom data (AUDIT, PHQ-8,
# PCL-5, PSQI) across three study phases.
#
# Key steps:
#   1. Multiple imputation of missing data via predictive mean matching
#   2. Define custom negative binomial response class for depmixS4
#   3. Fit 2- to 5-state HMMs with 100 random restarts per model
#   4. Select optimal model via BIC and state membership viability
#   5. Extract Viterbi-decoded state assignments and posterior probabilities
#
# Dependencies:
#   R 4.1.3 (arm64)
#   mice 3.18
#   depmixS4 1.5-0
#   gamlss 5.4-12
#   gamlss.dist 6.0-5
#
# Author: Mazen Istanbouli
# Date: 02/20/2026
# ==============================================================================

# --- Dependencies ---
library("mice")
library("depmixS4")
library("gamlss")
library("gamlss.dist")

# ==============================================================================
# 1. Multiple Imputation
# ==============================================================================

# Impute missing values using predictive mean matching (pmm) with 5 imputations
tempData <- mice(psych, m = 5, maxit = 50, 
                 meth = 'pmm', seed = 500)
summary(tempData)
completedData <- complete(tempData, 1)

# ==============================================================================
# 2. Custom Negative Binomial Response Class for depmixS4
#    Extends the response class to use the NBI (Type I) parameterization
#    from gamlss. Parameters mu and sigma are stored on the log scale to
#    ensure positivity during EM optimization.
# ==============================================================================

setClass("negbinom", contains = "response")

setGeneric("negbinom", function(y, pstart = NULL, fixed = NULL, ...) 
  standardGeneric("negbinom"))

# Constructor: initializes the negbinom response object with starting values
setMethod("negbinom",
          signature(y = "ANY"),
          function(y, pstart = NULL, fixed = NULL, ...) {
            y <- matrix(y, length(y))
            x <- matrix(1)
            parameters <- list()
            npar <- 2
            if(is.null(fixed)) fixed <- as.logical(rep(0, npar))
            if(!is.null(pstart)) {
              if(length(pstart) != npar) stop("length of 'pstart' must be ", npar)
              parameters$mu <- log(pstart[1])       # log-link for mean
              parameters$sigma <- log(pstart[2])     # log-link for dispersion
            }
            mod <- new("negbinom", parameters = parameters, fixed = fixed,
                       x = x, y = y, npar = npar)
            mod
          }
)

# Density: computes NBI density, exponentiating back to natural scale
setMethod("dens", "negbinom",
          function(object, log = FALSE) {
            dNBI(object@y, mu = exp(object@parameters$mu),
                 sigma = exp(object@parameters$sigma),
                 log = log)
          }
)

# Parameter extraction (getpars)
setMethod("getpars", "response",
          function(object, which = "pars", ...) {
            switch(which,
                   "pars" = {
                     parameters <- numeric()
                     parameters <- unlist(object@parameters)
                     pars <- parameters
                   },
                   "fixed" = {
                     pars <- object@fixed
                   }
            )
            return(pars)
          }
)

# Parameter setting (setpars)
setMethod("setpars", "negbinom",
          function(object, values, which = "pars", ...) {
            npar <- npar(object)
            if(length(values) != npar) stop("length of 'values' must be ", npar)
            nms <- names(object@parameters)
            switch(which,
                   "pars" = {
                     object@parameters$mu <- values[1]
                     object@parameters$sigma <- values[2]
                   },
                   "fixed" = {
                     object@fixed <- as.logical(values)
                   }
            )
            names(object@parameters) <- nms
            return(object)
          }
)

# Predict: returns the log-scale mean parameter
setMethod("predict", "negbinom",
          function(object) {
            ret <- object@parameters$mu
            return(ret)
          }
)

# Fit: called by the EM algorithm's M-step with posterior weights (gamma_t)
# Weights are clamped to avoid numerical instability in gamlss
setMethod("fit", "negbinom",
          function(object, w) {
            if(missing(w)) w <- NULL
            y <- object@y
            w[w < 1e-10] <- 0
            w[w > 1e10] <- 1e10
            fit <- gamlss(y ~ 1, weights = w, family = NBI(),
                          control = gamlss.control(n.cyc = 100, trace = FALSE),
                          mu.start = exp(object@parameters$mu),
                          sigma.start = exp(object@parameters$sigma))
            pars <- c(fit$mu.coefficients, fit$sigma.coefficients)
            object <- setpars(object, pars)
            object
          }
)

# ==============================================================================
# 3. Model Selection with Multiple Random Restarts
#    For each candidate number of states (2-5), the model is fit 100 times
#    with random starting values. The solution with the highest log-likelihood
#    is retained, and BIC is used for model comparison.
# ==============================================================================

set.seed(123)
n_restarts <- 100
state_range <- 2:5

results <- data.frame(nstates = integer(), 
                      best_ll = numeric(), 
                      bic = numeric())

best_models <- list()

start_time <- Sys.time()

for (nstate in state_range) {
  
  cat("\n=== Fitting", nstate, "state model ===\n")
  
  best_ll <- -Inf
  best_fit <- NULL
  
  for (r in 1:n_restarts) {
    
    rModels <- list()
    transition <- list()
    
    for (i in seq_len(nstate)) {
      rModels[[i]] <- list()
      rModels[[i]][[1]] <- negbinom(completedData$AUDIT_score, 
                                    pstart = c(round(mean(completedData$AUDIT_score), 0), 
                                               round(sd(completedData$AUDIT_score), 1)))
      rModels[[i]][[2]] <- negbinom(completedData$PSQI_score, 
                                    pstart = c(round(mean(completedData$PSQI_score), 0), 
                                               round(sd(completedData$PSQI_score), 1)))
      rModels[[i]][[3]] <- negbinom(completedData$PCL_score, 
                                    pstart = c(round(mean(completedData$PCL_score), 0), 
                                               round(sd(completedData$PCL_score), 1)))
      rModels[[i]][[4]] <- negbinom(completedData$PHQ_score, 
                                    pstart = c(round(mean(completedData$PHQ_score), 0), 
                                               round(sd(completedData$PHQ_score), 1)))
      transition[[i]] <- transInit(~1, nst = nstate, data = completedData)
    }
    
    inMod <- transInit(~1, ns = nstate, 
                       data = data.frame(matrix(1, nrow = nrow(completedData)/3, ncol = 1), 
                                         stringsAsFactors = FALSE), 
                       family = multinomial("identity"))
    
    mod <- makeDepmix(response = rModels, prior = inMod, transition = transition, 
                      ntimes = rep(3, 107), homogeneous = FALSE)
    
    fm <- tryCatch({
      fit(mod, verbose = FALSE, emc = em.control(random.start = TRUE))
    }, error = function(e) {
      cat("  Restart", r, "failed:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(fm)) {
      current_ll <- logLik(fm)
      if (current_ll > best_ll) {
        best_ll <- current_ll
        best_fit <- fm
        cat("  Restart", r, "- New best LL:", current_ll, "\n")
      }
    }
  }
  
  if (!is.null(best_fit)) {
    model_bic <- BIC(best_fit)
    results <- rbind(results, data.frame(nstates = nstate, 
                                         best_ll = best_ll, 
                                         bic = model_bic))
    best_models[[nstate]] <- best_fit
    cat("  Final best LL:", best_ll, "| BIC:", model_bic, "\n")
  }
}

end_time <- Sys.time()
tot_time <- round(end_time - start_time, 2)
print(tot_time)
print(results)

# ==============================================================================
# 4. Final 4-State Model (Selected)
#    The 4-state model was selected despite the 3-state model yielding a 
#    lower BIC, because the 3-state solution produced a degenerate state 
#    containing only 2 of 321 observations. The 4-state model provides 
#    balanced state membership (70, 36, 107, 108) suitable for downstream 
#    multi-omics analysis.
# ==============================================================================

nstate <- 4
rModels <- list()
transition <- list()
for(i in seq_len(nstate)){
  rModels[[i]] <- list()
  rModels[[i]][[1]] <- negbinom(completedData$AUDIT_score, pstart = c(round(mean(completedData$AUDIT_score),0), round(sd(completedData$AUDIT_score),1)))
  rModels[[i]][[2]] <- negbinom(completedData$PSQI_score, pstart = c(round(mean(completedData$PSQI_score),0), round(sd(completedData$PSQI_score),1)))
  rModels[[i]][[3]] <- negbinom(completedData$PCL_score, pstart = c(round(mean(completedData$PCL_score),0), round(sd(completedData$PCL_score),1)))
  rModels[[i]][[4]] <- negbinom(completedData$PHQ_score, pstart = c(round(mean(completedData$PHQ_score),0), round(sd(completedData$PHQ_score),1)))
  transition[[i]] <- transInit(~1, nst = nstate, data = completedData)
}
inMod <- transInit(~1, ns = nstate, data = data.frame(matrix(1, nrow = nrow(completedData)/3, ncol = 1), stringsAsFactors = FALSE), family = multinomial("identity"))
set.seed(123)
mod <- makeDepmix(response = rModels, prior = inMod, transition = transition, ntimes = rep(3,107), homogeneous = FALSE)
fm <- fit(mod, verbose = TRUE, emc = em.control(random.start = TRUE))
