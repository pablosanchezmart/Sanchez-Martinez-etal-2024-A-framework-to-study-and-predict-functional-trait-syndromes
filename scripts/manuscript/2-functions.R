########## FUNCTIONS ###############################################################################################################################

# Feb-2022
# Sanchez-Martinez, Pablo
# Ackerly, David


### Auxiliar functions --------------------------------------------------------- ####

# Complete data (inluded in package)
completePhyloData <- function(tree = tr, dataset = df, desiredCols = "SLA") {
  phyData <- list()
  completeVec <- complete.cases(dataset[, desiredCols])
  completeData <- dataset[completeVec, ]
  
  phylo <- keep.tip(tree, tree$tip.label[tree$tip.label %in% completeData$animal])
  phyData$dta <- as.data.frame(completeData[completeData$animal %in% phylo$tip.label, ])
  rownames(phyData$dta) <- phyData$dta$animal
  phyData$dta <- phyData$dta[phylo$tip.label, ]
  phyData$phylo <- phylo
  
  return(phyData)
}


# Calculate mean for numeric variables and unique values for descriptor factors
meanOrMode <- function(x){
  if(is.numeric(x) | is.double(x) | is.integer(x)){
    meanX <- mean(x)
    return(meanX)
  }
  if(is.character(x) | is.factor(x)){
    modeX <- unique(x)
    return(modeX)
  }
}


# optimize prior importance (specially for small datasets, where prior can have a strong influence on results)
optimizeMultiPriorFun <- function(pearsonCorrelation = NA,
                                  minCorDistance = 0.1, maxNumberIterations = 20, model.specifications = modelSpecifications,
                                  optiNiter = 1000, formula = fix.frml, dta = modellingData, variable1 = variable1, variable2 = variable2){
  print("optimizing priors")
  corDist <- abs(pearsonCorrelation)
  Niter <- 1
  toTheTop <- T
  toTheBottom <-  T
  
  print(paste0("corDist: ", round(corDist, 3), ", nu:", round(model.specifications$multiresponse_prior$G$G1$nu, 3)))
  
  while(corDist > minCorDistance  && Niter < maxNumberIterations){
    Niter <- Niter + 1
    
    if(toTheTop){
      
      # sum
      model.specifications$multiresponse_prior$G$G1$nu <- model.specifications$multiresponse_prior$G$G1$nu * 1.5
      
      mdl <- MCMCglmm(fixed = as.formula(formula), 
                      random = ~ us(trait):animal, rcov = ~us(trait):units, 
                      data= dta$dta, pedigree = dta$phylo, 
                      family = c("gaussian", "gaussian"), 
                      prior = model.specifications$multiresponse_prior,
                      nitt = optiNiter,
                      burnin = 10,
                      thin = 2,
                      verbose = F)
      
      CVphylo <- mdl$VCV[, paste0("trait", variable1, ":trait", variable2, ".animal")]
      CVres <- mdl$VCV[, paste0("trait", variable1, ":trait", variable2, ".units")]
      
      Vphylo1 <- mdl$VCV[, paste0("trait", variable1, ":trait", variable1, ".animal")]
      Vres1 <- mdl$VCV[, paste0("trait", variable1, ":trait", variable1, ".units")]
      
      Vphylo2 <- mdl$VCV[, paste0("trait", variable2, ":trait", variable2, ".animal")]
      Vres2 <- mdl$VCV[, paste0("trait", variable2, ":trait", variable2, ".units")]
      
      totalCorrlation_sum  <- (CVphylo + CVres) /
        sqrt( (Vphylo1 + Vres1) * (Vphylo2 + Vres2) )
      
      corDist_sum <-  abs(abs(pearsonCorrelation) - abs(mean(totalCorrlation_sum)))
    }
    
    if(toTheBottom){
      
      # substraction
      model.specifications$multiresponse_prior$G$G1$nu <- model.specifications$multiresponse_prior$G$G1$nu * 0.75
      
      mdl <- MCMCglmm(fixed = as.formula(formula), 
                      random = ~ us(trait):animal, rcov = ~us(trait):units, 
                      data= dta$dta, pedigree = dta$phylo, 
                      family = c("gaussian", "gaussian"), 
                      prior = model.specifications$multiresponse_prior,
                      nitt = model.specifications$number_interations,
                      burnin = model.specifications$burning_iterations,
                      thin = model.specifications$thinning_iterations,
                      verbose = F)
      
      CVphylo <- mdl$VCV[, paste0("trait", variable1, ":trait", variable2, ".animal")]
      CVres <- mdl$VCV[, paste0("trait", variable1, ":trait", variable2, ".units")]
      
      Vphylo1 <- mdl$VCV[, paste0("trait", variable1, ":trait", variable1, ".animal")]
      Vres1 <- mdl$VCV[, paste0("trait", variable1, ":trait", variable1, ".units")]
      
      Vphylo2 <- mdl$VCV[, paste0("trait", variable2, ":trait", variable2, ".animal")]
      Vres2 <- mdl$VCV[, paste0("trait", variable2, ":trait", variable2, ".units")]
      
      totalCorrlation_subs  <- (CVphylo + CVres) /
        sqrt( (Vphylo1 + Vres1) * (Vphylo2 + Vres2) )
      
      corDist_subs <-  abs(abs(pearsonCorrelation) - abs(mean(totalCorrlation_subs)))
    }
    
    if(corDist_sum < corDist_subs && corDist_sum < corDist){
      toTheTop <- T
      toTheBottom <- F
      corDist <- corDist_sum
    }
    if(corDist_subs < corDist_sum && corDist_subs < corDist){
      toTheTop <- F
      toTheBottom <- T
      corDist <- corDist_subs
    }
    
    print(paste0("corDist: ", round(corDist, 3), ", nu:", round(model.specifications$multiresponse_prior$G$G1$nu, 3)))
  }
  print(paste0("number of iterations: ", Niter))
  return(model.specifications)
}


### Simulate test data (included in package) ----------------------------------- ####

simulateDataSet <- function(tree = tr){
  
  require(castor)
  traitNames <- c("BM_HC_1", "BM_HC_2", "BM_HC_predictor", "BM_LC_1",  "BM_LC_2",   "BM_LC_predictor")
  diffMat <- matrix(c(1, 0.9, 0.8, 0, 0.1, 0.2,
                     0.9, 1, 0.8, 0, 0.1, 0.2,
                     0.8, 0.8, 1, 0, 0.1, 0.2,
                     0, 0, 0, 1, 0.9, 0.8,
                     0.1, 0.1, 0.1, 0.9, 1, 0.8,
                     0.2, 0.2, 0.2, 0.8, 0.8, 1
                    ), ncol = 6, dimnames = list(traitNames, traitNames)
                    )
  
  BM.df = as.data.frame(simulate_bm_model(tr, diffusivity=diffMat)$tip_states)
  colnames(BM.df) <- traitNames
  BM.df <- cbind("animal" = tr$tip.label, BM.df)
  
  rownames(BM.df) <- BM.df$animal
  
  nonBM.df <- rnorm_multi(n = length(tr$tip.label),
                         mu = c(0, 0, 0, 0, 0, 0),
                         sd = c(1, 1, 1, 1, 1, 1),
                         r = diffMat, 
                         varnames = c("nonBM_HC_1", "nonBM_HC_2", "nonBM_HC_predictor", "nonBM_LC_1",  "nonBM_LC_2",   "nonBM_LC_predictor"),
                         empirical = FALSE)
  
  
  nonBM.df <- cbind("animal" = tr$tip.label, nonBM.df)
  
  rownames(nonBM.df) <- nonBM.df$animal
  
  df <- merge(BM.df, nonBM.df, by = "animal")
  cor(df[, -1])
  return(df)
}


### Bayesian models diagnostics (included in the package) ---------------------- ####

## Model diagnostics (included in package)
diagnoseModels <- function(model = mdl){
  
  require(coda)
  
  model$autocFix <- autocorr.diag(model$Sol)
  model$autocRan <- autocorr.diag(model$VCV)
  model$heidelFix <- heidel.diag(model$Sol)
  model$heidelRan <- heidel.diag(model$VCV)
  model$effSizeFix <- effectiveSize(model$VCV)
  model$effSizeRan <- effectiveSize(model$Sol)
  name <- model$name
  
  mdlDiagn <- data.frame("model" = name, "AutcorrFix" = "T", "AutcorrRan" = "T", "HeidFix" = "T", "heidRan" = "T", "effSizeFix" = "T", "effSizeRan" = "T")
  
  # Autocorrelation
  if(any(abs(model$autocFix[-1, ]) > 0.1)){
    cat(name, " autocFix Failed. Increase thinning. \n")
    mdlDiagn$AutcorrFix <- "F"
  }
  if(any(abs(model$autocRan[-1, ]) > 0.1)){
    cat(name, " autocRan Failed. Increase thinning. \n")
    mdlDiagn$AutcorrRan <- "F"
  }
  # Convergence (stationary)
  i <- length(summary(model)$solutions)/5
  heidelFix <- model$heidelFix[1:i, 1:3]
  if(any(heidelFix == 0)){
    cat(name, " heidelFix Failed. Increase iterations/burning. \n")
    mdlDiagn$HeidFix <- "F"
  }
  i <- (length(summary(model)$Gcovariances)/4) + 1
  heidelRan <- model$heidelRan[1:i, 1:3]
  if(any(heidelRan == 0)){
    cat(name, " heidelRan Failed. Increase iterations/burning. \n")
    mdlDiagn$heidRan <- "F"
  }
  # Effect size
  if(any(model$effSizeFix < 1000)){
    cat(name, " EffSizeFix Failed. Increase iterations. \n")
    mdlDiagn$effSizeFix <- "F"
  }
  if(any(model$effSizeRan < 1000)){
    cat(name, " EffSizeRan Failed. Increase iterations.  \n")
    mdlDiagn$effSizeRan <- "F"
  }
  return(mdlDiagn)
}


### Define model specifications (included in package)--------------------------- ####

# Function to set mcmcglmm specifications
# 
# Arguments:
# 
# nitter: number of iterations (integer)
# burni: number of iterations burned (to ensure convergence) (integer)
# thinni: iterations not sampled (e.g., thinni = 10 will sample 1 of each 10 contiguous iterations, to avoid mcmc autocorrelation) (integer)
# Prior2: non-informative prior for phylogenetic effect and residual variance. Uniresponse models. (numeric matrix)
# priorMulti2: non-informative prior for phylogenetic effect and residual variance. Biresponse models. (numeric matrix)
# 
# Output:
# 
# models.specifications.list: list with model specifications
# 
# Pablo Sanchez-Martinez, 2022
# p.sanchez@creaf.uab.cat

defineModelsSpecifications <- function(nitter = 100, burni = 10, thinni = 2,
                                          # One response variable, phylogeny and residuals
                                          prior2 = list(
                                            R = list(V = 1, nu = 0.002), 
                                            G = list(G1 = list(V = 1, nu = 0.002))
                                          ),
                                          # Two response variables, phylogeny and residuals
                                          priorMulti2 = list(
                                            R=list(V=diag(2)/2,nu=2), 
                                            G=list(G1=list(V=diag(2)/2,nu=2))
                                          )
) {
  
  models.specifications.list <- list("number_interations" = nitter, 
                                     "burning_iterations" = burni, 
                                     "thinning_iterations" = thinni, 
                                     "uniresponse_prior" = prior2, 
                                     "multiresponse_prior" = priorMulti2)
  return(models.specifications.list)
}
