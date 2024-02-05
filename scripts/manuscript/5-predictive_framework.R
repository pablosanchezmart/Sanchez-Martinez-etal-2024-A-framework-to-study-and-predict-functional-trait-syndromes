#### PREDICTIONS --------------------------------------------------------------- ####

# Pablo Sanchez Martinez

# 22-11-2021

# dataGroup <- "traits"

for(propNA in c(0.2, 0.5)){ #, 0.7
  for(dataGroup in c("traits")){  #"simulations"
    
    source("scripts/manuscript/1-init.R")
    
    print(paste0("Extracting results for: ", dataGroup))
    
    #### Variance-covariance structures ------------------------------------------ ####
    
    for(f in list.files(path = paste0("outputs/outputs_", dataGroup, "/models_outputs/"),  full.names = T)){
      load(file = f)
    }
    
    #### DATA -------------------------------------------------------------------- ####
    
    df_pred <- df %>% filter(animal %in% tr$tip.label)
    df_pred <- df_pred %>% dplyr::rename(taxon = animal) %>% dplyr::select(taxon, traits, all_of(predictors)) %>% as.data.frame()
    
    # delete observations with no trait at all
    df_pred$allTraitsNA <- 0
    for(j in 1:nrow(df_pred)) {
      df_pred$allTraitsNA[j] <- ifelse(all(is.na(df[j, traits])), 1, 0)
    }
    df_pred <- df_pred %>% filter(allTraitsNA == 0) %>% select(-allTraitsNA)
    
    # only complete cases for traits
    # df_pred <- df_pred[complete.cases(df_pred), ]
    
    # delete observations with no environmental variables
    for(prd in predictors){
      df_pred <- df_pred[which(!is.na(df_pred[, prd])), ]
    }
    
    tr_pred <- keep.tip(tr, tr$tip.label[tr$tip.label %in% df_pred$taxon])
    
    
    #### IMPUTATION METHODS COMPARISON ------------------------------------------- ####
    
    performance_results <- data.frame("Type" = character(), "Variable" = character(), "NRMSE" = numeric(), "R2" = numeric())
    
    for(i in 1:numberIterations){
      
      # NA production
      df_miss <- df_pred
      df_miss[, traits] <- missForest::prodNA(as.data.frame(df_pred[, traits]), propNA)
      df_imp_mean <- df_miss
      df_imp_mice <- df_miss
      
      #### MEAN METHODOLOGY ------------------------------------------------------ ####
      
      for(trait in traits){
        # imputation
        df_imp_mean[, trait][which(is.na(df_imp_mean[, trait]))] <- mean(df_miss[, trait], na.rm = TRUE)
      }
      
      #### MICE METHODOLOGY ------------------------------------------------------ ####
      
      df_imp_mice$taxon <- NULL
      tempData <- mice(df_imp_mice, m=5, maxit=50, meth='pmm', seed=500)
      df_imp_mice <- complete(tempData, 1)
      
      
      #### RPHYLOPARS METHODOLOGY ------------------------------------------------ ####
      
      df_miss_phylopars <- cbind("species" = df_miss$taxon, df_miss[, traits])
      phyloPars_temp <- phylopars(trait_data = df_miss_phylopars, tree = tr_pred)
      
      df_imp_phylopars <- as.data.frame(phyloPars_temp$anc_recon[df_miss$taxon, ])
      df_imp_phylopars$taxon <- df_miss$taxon
      
      
      #### TrEvol METHODOLOGY ---------------------------------------------------- ####
      
      df_imp_TrEvol <- imputeTraits(imputationVariables = traits, predictors = predictors,
                                    dataset = df_miss, phylogeny = tr_pred, 
                                    # numberOfPhyloCoordinates = nPhyloCoord,
                                    correlationsTraitsResults = traitsCovariancePartitionResults,
                                    varianceResults = traitsVariancePartitionResults,
                                    # orderCriterium = "Total_coordination",
                                    numberOfPhyloCoordinates = 5,
                                    prodNAs = 0, IterationsNumber = 5, clustersNumber = clNumber, forceRun = forceRunImputation, 
                                    numberOfImputationRounds = numberOfImpRounds)
      
      
      #### PERFORMANCE RESULTS --------------------------------------------------- ####
      
      for(trait in traits){
        
        performance_results <- rbind(performance_results, 
                                     cbind(modelPerformanceCalculation(type = "Mean", 
                                                                 miss = df_miss,
                                                                 imp = df_imp_mean,
                                                                 true = df_pred,
                                                                 trait = trait),
                                           "formula" = paste0("mean(", trait, ")")
                                     )
        )
        
        performance_results <- rbind(performance_results, 
                                    cbind( modelPerformanceCalculation(type = "MICE",
                                                                 miss = df_miss,
                                                                 imp = df_imp_mice,
                                                                 true = df_pred, 
                                                                 trait = trait),
                                           "formula" = paste0("MICE(", paste0(c(traits, predictors),  collapse = ", "), ")")
                                    )
        )
        
        performance_results <- rbind(performance_results, 
                                     cbind(modelPerformanceCalculation(type = "phylopars", 
                                                                 miss = df_miss, 
                                                                 imp = df_imp_phylopars, 
                                                                 true = df_pred, 
                                                                 trait = trait),
                                           "formula" = paste0("phylopars(", paste0(traits, collapse = ", "), ")")
                                     )
        )
        
        for(nround in 1:numberOfImpRounds){
          
          performance_results <- rbind(performance_results,
                                       cbind(modelPerformanceCalculation(type = paste0("TrEvol Imputation ", nround), 
                                                                   miss = df_miss, 
                                                                   imp = df_imp_TrEvol[[paste0("round", nround)]]$ximp, 
                                                                   true = df_pred, 
                                                                   trait = trait),
                                             "formula" = df_imp_TrEvol[[paste0("round", nround)]]$modelFormula
                                       )
          )
          
        }
        
      }
      
    }
    
    performance_results <- plyr::arrange(performance_results, Type)
    
    save.image(file = paste0(outputs.dir, "/predictions/imputation_comparisons_nIter_",numberIterations, "NAprop_", propNA, ".RData"))
    print(paste0(outputs.dir, "/predictions/imputation_comparisons_nIter_",numberIterations, "NAprop_", propNA, ".RData"))
    
  }
}


#### PLOT RESULTS MODEL PERFORMANCES ------------------------------------------- ####

# dataGroup <- "simulations"

for(propNA in c(0.2, 0.5)){
  for(dataGroup in c("simulations", "traits")){
    
    source("scripts/manuscript/1-init.R")
    
    print(paste0("Extracting results for: ", dataGroup))
    
    #### TrEvol METHODOLOGY ------------------------------------------------------ ####
    
    #### Predictions results ####
    
    load(file = paste0(outputs.dir, "/predictions/imputation_comparisons_nIter_",numberIterations, "NAprop_", propNA, ".RData"))
    print(paste0(outputs.dir, "/predictions/imputation_comparisons_nIter_",numberIterations, "NAprop_", propNA, ".RData"))
    
    source("scripts/manuscript/1-init.R")
    
    w <- 3
    h <- 3
    
    performance_results$Variable <- str_replace_all(performance_results$Variable, "_", " ")
    
    performance_results_trEvol <- performance_results %>% filter(!is.na(str_extract(Type, "TrEvol")))
    performance_results_trEvol
    
    performance_results_plot <- performance_results %>% filter(Type %in% c("Mean", "MICE", "phylopars", paste0("TrEvol Imputation ", numberOfImpRounds)))
    
    # for(vari in unique(performance_results_trEvol$Variable)){
    #   
    #   rnd <- performance_results_trEvol %>% filter(Variable == vari)
    #   rnd <-  rnd %>% filter(R2 == max(rnd$R2)) %>% pull(Type)
    #   
    #   rnd_rslts <- performance_results_trEvol %>% filter(Type == rnd, Variable == vari)
    #   rnd_rslts$Type <- "TrEvol Imputation"
    #   performance_results_plot <- rbind(performance_results_plot, rnd_rslts)
    # }
    
    
    ### TrEvol imputation rounds performance ####
    
    ## NRMSE
    
    pimp1_vs2_nrmse_trEvol <- ggplot2::ggplot(performance_results_trEvol, ggplot2::aes(x = Variable, y = NRMSE, color = Type)) +
      ggplot2::geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
      ggplot2::geom_boxplot() + theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
    pimp1_vs2_nrmse_trEvol
    
    pdf(file = paste0(results.dir, "/figures/predictions/imputation_TrEvol_rounds_nrmse_nIter_", numberIterations, "NAprop_", propNA,  ".pdf"), height = h, width = w)
    print(pimp1_vs2_nrmse_trEvol)
    dev.off()
    print(paste0(results.dir, "/figures/predictions/imputation_TrEvol_rounds_nrmse_nIter_", numberIterations, "NAprop_", propNA,  ".pdf"))
    
    ## R2
    
    pimp1_vs2_r2_trEvol <- ggplot2::ggplot(performance_results_trEvol, ggplot2::aes(x = Variable, y = R2, color = Type)) +
      ggplot2::geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
      ggplot2::geom_boxplot() + ylim(0, 1)
    pimp1_vs2_r2_trEvol
    
    legend1 <- get_legend(pimp1_vs2_r2_trEvol)
    
    pimp1_vs2_r2_trEvol <- pimp1_vs2_r2_trEvol + theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
    
    pdf(file = paste0(results.dir, "/figures/predictions/imputation_rounds_r2_nIter_", numberIterations, "NAprop_", propNA,  ".pdf"), height = h, width = w)
    print(pimp1_vs2_r2_trEvol)
    dev.off()
    print(paste0(results.dir, "/figures/predictions/imputation_rounds_r2_nIter_", numberIterations, "NAprop_", propNA,  ".pdf"))
    
    ## legend
    
    pdf(file = paste0(results.dir, "/figures/predictions/legend_rounds.pdf"))
    plot(legend1)
    dev.off()
    print(paste0(results.dir, "/figures/predictions/legend_rounds.pdf"))
    
    
    ### Mean, MICE, TrEvol imp1 and TrEvol impo2 comparation ####
    
    ## NRMSE
    
    pimp1_vs2_nrmse <- ggplot2::ggplot(performance_results_plot, ggplot2::aes(x = Variable, y = NRMSE, color = Type)) +
      ggplot2::geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
      ggplot2::geom_boxplot() + theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank()) + ylim(0, 1.5)
    pimp1_vs2_nrmse
    
    pdf(file = paste0(results.dir, "/figures/predictions/imputation_comparison_nrmse_nIter_", numberIterations, "NAprop_", propNA,  ".pdf"), height = h, width = w)
    print(pimp1_vs2_nrmse)
    dev.off()
    print(paste0(results.dir, "/figures/predictions/imputation_comparison_nrmse_nIter_", numberIterations, "NAprop_", propNA,  ".pdf"))
    
    ## R2
    
    pimp1_vs2_r2 <- ggplot2::ggplot(performance_results_plot, ggplot2::aes(x = Variable, y = R2, color = Type)) +
      ggplot2::geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
      ggplot2::geom_boxplot() + ylim(0, 1)
    pimp1_vs2_r2
    
    legend1 <- get_legend(pimp1_vs2_r2)
    
    pimp1_vs2_r2 <- pimp1_vs2_r2 + theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
    
    pdf(file = paste0(results.dir, "/figures/predictions/imputation_comparison_r2_nIter_", numberIterations, "NAprop_", propNA,  ".pdf"), height = h, width = w)
    print(pimp1_vs2_r2)
    dev.off()
    print(paste0(results.dir, "/figures/predictions/imputation_comparison_r2_nIter_", numberIterations, "NAprop_", propNA,  ".pdf"))
    
    ## legend
    
    pdf(file = paste0(results.dir, "/figures/predictions/legend_comparison.pdf"))
    plot(legend1)
    dev.off()
    print(paste0(results.dir, "/figures/predictions/legend_comparison.pdf"))
  }
} #, 0.7
