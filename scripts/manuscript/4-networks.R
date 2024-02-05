########## RESULTS #############################################################

# Pablo Sanchez Martinez

# 22-11-2021

remove(list = ls())

# dataGroup <-  "traits"
# dataGroup <-  "simulations"

for(dataGroup in c("simulations", "traits")){
  
  if(dataGroup == "simulations"){
    onlySig <- T
    edgeLab <- T
    threshold <- 0.2
  } else{
    onlySig <- T
    edgeLab <- T
    threshold <- 0
  }
  
  source("scripts/manuscript/1-init.R")
  
  print(paste0("Extracting results for: ", dataGroup))
  
  #### LOAD MODELS OUTPUTS ----------------------------------------------------- ####
  
  variancePartition.list <- list()
  covariancePartition.list <- list()
  
  filenames <- list.files(path = paste0("outputs/outputs_", dataGroup, "/models_outputs/"), pattern = ".RData",  full.names = F)
  
  for(i in 1:length(filenames)){
    
    fileName <- paste0("outputs/outputs_", dataGroup, "/models_outputs/", filenames[i])
    
    name <- str_remove(filenames[i], ".RData")
    
    if(str_detect(fileName, "traitsVariancePartition")){
      variancePartition.list[name] <- mget(load(file = fileName))
      
      write.csv(variancePartition.list[[name]]$varianceResults, paste0(results.dir, "/tables/", name, ".csv"), row.names = F)
      print(paste0(results.dir, "/tables/", name, ".csv"))
    }
    
    if(str_detect(fileName, "traitsCovariancePartition")){
      covariancePartition.list[name] <- mget(load(file = fileName))
      
      write.csv(covariancePartition.list[[name]]$covarianceResults, paste0(results.dir, "/tables/", name, ".csv"), row.names = F)
      print(paste0(results.dir, "/tables/", name, ".csv"))
    }
  }
  
  
  #### VCV NETWORK PLOTS ------------------------------------------------------- ####
  
  grouping_vars <- data.frame("variable" = c(traits, predictors), "group" = c(rep("trait", length(traits)), rep("predictor", length(predictors))))
  
  ### Total correlation
  
  pdf(paste0(results.dir, "/figures/VCV_networks/total_coordination.pdf"), width = w, height = h)
  TC <- plotNetwork(covariance.results = covariancePartition.list$traitsCovariancePartitionResults,
                    covariance.type = "Total_coordination", group.variables = grouping_vars, edge.label = edgeLab, layout = qgraphLayout,
                    only.significant = onlySig, threshold = threshold, node.label = allNodeLabels, networkMetrixTextSize = txtSize)
  dev.off()
  print(paste0(results.dir, "/figures/VCV_networks/total_coordination.pdf"))
  
  
  ### Total correlation of traits
  
  pdf(paste0(results.dir, "/figures/VCV_networks/total_coordination_traits.pdf"), width = w, height = h)
  TC_traits <- plotNetwork(covariance.results = covariancePartition.list$traitsCovariancePartitionResults,
                    covariance.type = "Total_coordination", group.variables = grouping_vars, edge.label = edgeLab, layout = qgraphLayout,
                    only.significant = onlySig, threshold = threshold, node.label = allNodeLabels, networkMetrixTextSize = txtSize,
                    not.show.variables = predictors)
  dev.off()
  print(paste0(results.dir, "/figures/VCV_networks/total_coordination_traits.pdf"))
  
  
  
  ### Phylogenetic correlation (joined phylogenetic signal)
  # 
  # pdf(paste0(results.dir, "/figures/VCV_networks/total_phylogenetic_conservatism.pdf"), width = w, height = h)
  # PC <- plotNetwork(covariance.results = covariancePartition.list$traitsCovariancePartitionResults$covarianceResults,
  #                             variance.results = variancePartition.list$traitsVariancePartitionResults$varianceResults,
  #                             variance.type = "Total_phylogenetic_conservatism", covariance.type = "Total_coordinated_phylogenetic_conservatism",
  #                             group.variables = grouping_vars, edge.label = F, layout = qgraphLayout,
  #                             only.significant = onlySig, threshold = threshold, node.label = allNodeLabels, not.show.variables = c("BM_HC_predictor", "nonBM_HC_predictor", "ln_MAT", "ln_AP"))
  # dev.off()
  # print(paste0(results.dir, "/figures/VCV_networks/total_phylogenetic_conservatism.pdf"))
  
  ### Residual (non-phylogenetic) correlation
  
  # pdf(paste0(results.dir, "/figures/VCV_networks/total_coordinated_radiation.pdf"), width = w, height = h)
  # RC <- plotNetwork(covariance.results = covariancePartition.list$traitsCovariancePartitionResults$covarianceResults,
  #                             variance.results = variancePartition.list$traitsVariancePartitionResults$varianceResults,
  #                             variance.type = "Total_phylogenetic_conservatism", covariance.type = "Total_coordinated_radiation",
  #                             group.variables = grouping_vars, edge.label = T, layout = qgraphLayout,
  #                             only.significant = onlySig, node.label = allNodeLabels, not.show.variables = c("BM_HC_predictor", "nonBM_HC_predictor"))
  # dev.off()
  # print(paste0(results.dir, "/figures/VCV_networks/total_coordinated_radiation.pdf"))
  
  
  #### PARTIAL VCV NETWORK PLOTS ----------------------------------------------- ####
  
  for(predictor in predictors){
    
    ### Total correlation
    
    pdf(paste0(results.dir, "/figures/VCV_networks/total_coordination_", predictor, ".pdf"), width = w, height = h)
    TC <- plotNetwork(covariance.results = covariancePartition.list$traitsCovariancePartitionResults,
                      covariance.type = "Total_coordination", group.variables = grouping_vars, edge.label = edgeLab, layout = qgraphLayout,
                      only.significant = onlySig, threshold = threshold, node.label = allNodeLabels, not.show.variables = predictors[-which(predictors == predictor)], 
                      networkMetrixTextSize = txtSize)
    dev.off()
    print(paste0(results.dir, "/figures/VCV_networks/total_coordination_", predictor, ".pdf"))
    
    
    ### Pure phylogenetic conservatism ####
    
    pdf(paste0(results.dir, "/figures/VCV_networks/pure_phylogenetic_conservatism_", predictor, ".pdf"), width = w, height = h)
    PPC <- plotNetwork(covariance.results = covariancePartition.list[[paste0("traitsCovariancePartitionResults_", predictor)]],
                                variance.results = variancePartition.list[[paste0("traitsVariancePartitionResults_", predictor)]],
                                variance.type = "Pure_phylogenetic_conservatism", covariance.type = "Pure_coordinated_phylogenetic_conservatism", 
                                group.variables = grouping_vars, edge.label = edgeLab, layout = qgraphLayout,
                                only.significant = onlySig, threshold = threshold, node.label = traitsNodeLabels, networkMetrixTextSize = txtSize)
    dev.off()
    print(paste0(results.dir, "/figures/VCV_networks/pure_phylogenetic_conservatism_", predictor, ".pdf"))
    
    
    ### phylogenetic niche conservatism ####
    
    pdf(paste0(results.dir, "/figures/VCV_networks/phylogenetic_niche_conservatism_", predictor, ".pdf"), width = w, height = h)
    PNC <- plotNetwork(covariance.results = covariancePartition.list[[paste0("traitsCovariancePartitionResults_", predictor)]],
                       variance.results = variancePartition.list[[paste0("traitsVariancePartitionResults_", predictor)]],
                       variance.type = "Phylogenetic_niche_conservatism", covariance.type = "Coordinated_phylogenetic_niche_conservatism", 
                       group.variables = grouping_vars, edge.label = edgeLab, layout = qgraphLayout,
                       only.significant = onlySig, threshold = threshold, node.label = traitsNodeLabels, networkMetrixTextSize = txtSize)
    dev.off()
    print(paste0(results.dir, "/figures/VCV_networks/phylogenetic_niche_conservatism_", predictor, ".pdf"))

    
    ### total environmental coordiantion ####
    # 
    # pdf(paste0(results.dir, "/figures/VCV_networks/total_environmental_coordination_", predictor, ".pdf"), width = w, height = h)
    # TE <- plotNetwork(covariance.results = covariancePartition.list[[paste0("traitsCovariancePartitionResults", predictor)]]$covarianceResults,
    #                   variance.results = variancePartition.list[[paste0("traitsVariancePartitionResults", predictor)]]$varianceResults,
    #                   variance.type = "Total_environmental", covariance.type = "Total_environmental_coordination", 
    #                   group.variables = grouping_vars, edge.label = F, layout = qgraphLayout,
    #                   only.significant = onlySig, threshold = threshold, node.label = traitsNodeLabels)
    # dev.off()
    # print(paste0(results.dir, "/figures/VCV_networks/total_environmental_coordination_", predictor, ".pdf"))
    
    
    ### pure environmental coordiantion ####
    
    pdf(paste0(results.dir, "/figures/VCV_networks/pure_environmental_coordination_", predictor, ".pdf"), width = w, height = h)
    PE <- plotNetwork(covariance.results = covariancePartition.list[[paste0("traitsCovariancePartitionResults_", predictor)]],
                       variance.results = variancePartition.list[[paste0("traitsVariancePartitionResults_", predictor)]],
                       variance.type = "Pure_environmental", covariance.type = "Pure_environmental_coordination", 
                       group.variables = grouping_vars, edge.label = edgeLab, layout = qgraphLayout,
                       only.significant = onlySig, threshold = threshold,  node.label = traitsNodeLabels, networkMetrixTextSize = txtSize)
    dev.off()
    print(paste0(results.dir, "/figures/VCV_networks/pure_environmental_coordination_", predictor, ".pdf"))
    
    ### residual coordiantion ####
    
    pdf(paste0(results.dir, "/figures/VCV_networks/residual_coordination_", predictor, ".pdf"), width = w, height = h)
    RES <- plotNetwork(covariance.results = covariancePartition.list[[paste0("traitsCovariancePartitionResults_", predictor)]],
                      variance.results = variancePartition.list[[paste0("traitsVariancePartitionResults_", predictor)]],
                      variance.type = "Residual", covariance.type = "Residual_coordination", 
                      group.variables = grouping_vars, edge.label = edgeLab, layout = qgraphLayout,
                      only.significant = onlySig, threshold = threshold,  node.label = traitsNodeLabels, networkMetrixTextSize = txtSize)
    dev.off()
    print(paste0(results.dir, "/figures/VCV_networks/residual_coordination_", predictor, ".pdf"))
  }
  
  
  #### CORRELATION IMPORTANCES ------------------------------------------------- ####

  # pdf(paste0(results.dir, "/figures/VCV_networks/proportion_correlations.pdf"), width = w, height = h)
  # plotVcv(correlations = correlation.list$correlationsResults$correlation.results, phylogenetic.signal = phyloSignal.list$phylogeneticSignalResults$phylogenetic.signal.results, 
  #         order_vars = grouping_vars$variable)
  # dev.off()
  # print(paste0(results.dir, "/figures/VCV_networks/proportion_correlations.pdf"))
  
}
