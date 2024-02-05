########## COMPARISON BETWEEN BI AND UNIRESPONSE MCMCGLMM MODELS ###############

# Pablo Sanchez Martinez

# 22-11-2021

# dataGroup <- "traits"

for(dataGroup in c("simulations", "traits")){
  
  print(paste0("Running models for ", dataGroup))
  
  source("scripts/manuscript/1-init.R")

  #### PHYLOGENETIC SIGNAL ------------------------------------------------------- ####

  # variancePartitionResults <- computeVariancePartition(traits = vars, dataset = df, phylogeny = tr,
  #                                                       model.specifications = modelSpecifications, force.run = forceRunModels)
  # 
  # #### CORRELATIONS -------------------------------------------------------------- ####
  # 
  # covariancePartitionResults <- computeVariancePartition(traits = vars, dataset = df, phylogeny = tr,
  #                                                   model.specifications = modelSpecifications, force.run = forceRunModels)
  
  
  covariancePartitionResults <- computeVarianceCovariancePartition(traits = vars, dataset = df, phylogeny = tr,
                                                                   model.specifications = modelSpecifications, force.run = forceRunModels)

  for(predictor in predictors){
    
    #### PARTIAL RESIDUAL PHYLOGENETIC SIGNAL | PREDICTOR ------------------------ ####
    
    # PartialVariancePartitionResults <- computeVariancePartition(traits = traits, environmental.variables = predictor, dataset = df, phylogeny = tr, 
    #                                                              model.specifications = modelSpecifications, force.run = forceRunModels)
    # 
    # write.csv(PartialVariancePartitionResults$varianceResults, 
    #           paste0(outputs.dir, "/phylogenetic_variance_covariance/PartialVariancePartitionResults_", predictor, ".csv"), row.names = F)
    # print(paste0(outputs.dir, "/phylogenetic_variance_covariance/PartialVariancePartitionResults_", predictor, ".csv"))
    # 
    # #### RESIDUAL CORRELATIONS | PREDICTOR---------------------------------------- ####
    # 
    # PartialCovariancePartitionResults <- computeCovariancePartition(traits = traits, environmental.variables = predictor, dataset = df, phylogeny = tr, 
    #                                                             model.specifications = modelSpecifications, force.run = forceRunModels)
    # 
    # write.csv(PartialCovariancePartitionResults$covarianceResults, paste0(outputs.dir, "/phylogenetic_variance_covariance/PartialCovariancePartitionResults_", predictor, ".csv"), row.names = F)
    # print(paste0(outputs.dir, "/phylogenetic_variance_covariance/PartialCovariancePartitionResults_", predictor, ".csv"))
    
    
    #### ALTOGETHER ####
    
    
    PartialVarianceCovariancePartitionResults <- computeVarianceCovariancePartition(traits = traits, environmental.variables = predictor, dataset = df, phylogeny = tr, 
                                                                    model.specifications = modelSpecifications, force.run = forceRunModels)
    
    write.csv(PartialVarianceCovariancePartitionResults$varianceResults, 
              paste0(outputs.dir, "/phylogenetic_variance_covariance/PartialVariancePartitionResults_", predictor, ".csv"), row.names = F)
    print(paste0(outputs.dir, "/phylogenetic_variance_covariance/PartialVariancePartitionResults_", predictor, ".csv"))
    
    write.csv(PartialVarianceCovariancePartitionResults$covarianceResults, paste0(outputs.dir, "/phylogenetic_variance_covariance/PartialCovariancePartitionResults_", predictor, ".csv"), row.names = F)
    print(paste0(outputs.dir, "/phylogenetic_variance_covariance/PartialCovariancePartitionResults_", predictor, ".csv"))
}

}