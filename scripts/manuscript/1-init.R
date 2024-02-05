####################### ANALYSES INITIALIZAITON #########################################################################################

# Sanchez-Martinez Pablo February 2022


print("loading workspace ...")

#### PACKAGES ------------------------------------------------------------------ ####

# packagesFun <- function(package){
#   if(!require(package, character.only = T)){
#     install.packages(package, dependencies = T)
#     library(package, character.only = T)
#   } else{
#     library(package, character.only = T)
#   }
#
# }
#
# packages <- c("plyr", "dplyr", "tidyverse", "caret", "stringr", "MCMCglmm", "bayestestR", "phytools", "parallel", "foreach", "RColorBrewer", "ggplot2",
# "ggpubr", "ggtree", "ggstance", "ape", "GGally", "network", "sna", "corrr", "intergraph", "qgraph", "igraph", "CINNA", "factoextra", "corrplot",
# "faux", "PCMBaseCpp", "mvSLOUCH", "OUwie", "bayou", "missForest", "doParallel", "gtools")
#
# for(package in packages){
#   packagesFun(package)
# }
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("ggtree")

# remove.packages("TrEvol")
# devtools::install_github("https://github.com/pablosanchezmart/TrEvol")

library(TrEvol)

library(readr)
library(ggtree)
library(dplyr)
library(tidyverse)
library(qgraph)
library(ape)
library(phytools)
library(gtools)
library(ggpubr)

library(parallel)

library(mice)
library(Rphylopars)
library(phytools)
# remotes::install_github("saskiaotto/INDperform")
library(INDperform)

initializeTrEvo(folder.name = dataGroup)

forceRunModels <- T

### GENERAL FUNCTIONS ---------------------------------------------------------- ####

# Function to plot scatterplot with phylogeny

plotTraitsPhylo <- function (tree, coords, rotate = TRUE, colourVar, ...) 
{
  if (hasArg(xlim)) {
    xlim <- list(...)$xlim
  } else range(coords[, 1])
  if (hasArg(ylim)){
    ylim <- list(...)$ylim
  } else ylim <- range(coords[, 2])* 2
  if (hasArg(fsize)){
    fsize <- list(...)$fsize
  } else fsize <- 1
  if (hasArg(split)){
    split <- list(...)$split
  } else split <- c(0.4, 0.6)
  if (hasArg(psize)){
    psize <- list(...)$psize
  }else psize <- 1
  if (hasArg(cex.points)) {
    cex.points <- list(...)$cex.points
    if (length(cex.points) == 1) 
      cex.points <- c(0.6 * cex.points, cex.points)
  } else cex.points <- c(0.6 * psize, psize)
  if (hasArg(mar)){
    mar <- list(...)$mar
  } else mar <- rep(0, 4)
  if (hasArg(asp)){
    asp <- list(...)$asp
  } else asp <- 1
  if (hasArg(ftype)){
    ftype <- list(...)$ftype
  } else ftype <- "reg"
  ftype <- which(c("off", "reg", "b", "i", "bi") == ftype) - 
    1
  if (!ftype){
    fsize = 0
  }
  
  if (hasArg(from.tip)){
    from.tip <- list(...)$from.tip
  } else from.tip <- FALSE
  if (hasArg(colors)){
    colors <- list(...)$colors
  } else colors <- "red"
  if (length(colors) == 1){
    colors <- rep(colors[1], 2)
  }
  
  type <- "phylogram"
  
  if (length(colors) == 2 && type == "phylogram") {
    colors <- matrix(rep(colors, nrow(coords)), nrow(coords), 
                     2, byrow = TRUE)
    rownames(colors) <- rownames(coords)
  } else if (is.vector(colors) && (length(colors) == Ntip(tree))) {
    COLS <- matrix("red", nrow(coords), 2, dimnames = list(rownames(coords)))
    for (i in 1:length(colors)) COLS[which(rownames(COLS) == 
                                             names(colors)[i]), 1:2] <- colors[i]
    colors <- COLS
  }
  
  if (hasArg(pch)){
    pch <- list(...)$pch
  } else pch <- 21
  if (length(pch) == 1){
    pch <- rep(pch, 2)
  } 
  
  if (hasArg(lwd)){
    lwd <- list(...)$lwd
  } else lwd <- c(2, 1)
  if (length(lwd) == 1){
    lwd <- rep(lwd, 2)
  } 
  
  if (hasArg(lty)){
    lty <- list(...)$lty
  } else lty <- "dashed"
  if (hasArg(pts)){
    pts <- list(...)$pts
  } else pts <- TRUE
  if (hasArg(col.edge)){
    col.edge <- list(...)$col.edge
  } else col.edge <- rep(par()$fg, nrow(tree$edge))
  if (hasArg(map.bg)){
    map.bg <- list(...)$map.bg
  } else map.bg <- "gray95"
  
  ylim <- c(ylim[1], ylim[2] + 0.03 * diff(ylim))
  ylim <- c(ylim[1], ylim[2] + split[1]/split[2] * (ylim[2] - ylim[1]))
  
  ylim[1] <- min(coords[, 2]) - 0.5
  
  if (all(mar == 0)){
    mar <- mar + 3
  }
  
  plot.new()
  par(mar = mar)
  plot.window(xlim = xlim, ylim = ylim, asp = asp)
  
  cw <- reorder(tree, "cladewise")
  n <- Ntip(cw)
  
  dx <- abs(diff(xlim))
  rect(xlim[1] - 1.04 * dx, ylim[2] - split[1] * (ylim[2] - 
                                                    ylim[1]), xlim[2] + 1.04 * dx, ylim[2], col = par()$bg, border = par()$bg)
  pdin <- par()$din[2]
  sh <- (fsize * strwidth(paste(" ", cw$tip.label, sep = "")) + 0.3 * fsize * strwidth("W")) * (par()$din[1]/par()$din[2]) * (diff(par()$usr[3:4])/diff(par()$usr[1:2]))
  cw$edge.length <- cw$edge.length/max(nodeHeights(cw)) *   (split[1] * (ylim[2] - ylim[1]) - max(sh))
  pw <- reorder(cw, "postorder")
  x <- vector(length = n + cw$Nnode)
  x[cw$edge[cw$edge[, 2] <= n, 2]] <- 0:(n - 1)/(n - 1) * (xlim[2] - xlim[1]) + xlim[1]
  nn <- unique(pw$edge[, 1])
  for (i in 1:length(nn)) {
    xx <- x[pw$edge[which(pw$edge[, 1] == nn[i]), 
                    2]]
    x[nn[i]] <- mean(range(xx))
  }
  Y <- ylim[2] - nodeHeights(cw)
  coords <- coords[, 2:1]
  for (i in 1:nrow(coords)) {
    tip.i <- which(cw$tip.label == rownames(coords)[i])
    lines(c(x[tip.i], coords[i, 1]), c(Y[which(cw$edge[, 
                                                       2] == tip.i), 2] - if (from.tip) 0 else sh[tip.i], 
                                       coords[i, 2]), col = colors[i, 1], lty = lty, 
          lwd = lwd[2])
  }
  
  color <- as.character(cut(colourVar, breaks = quantile(colourVar), include.lowest = TRUE, labels = colorspace::sequential_hcl(4)))
  
  points(coords, pch = 19, cex = cex.points[2], col = color)
  
  axis1 <- c(round(quantile(coords)))
  axis(1, at = axis1)
  mtext(stringr::str_replace_all(colnames(coords)[1], "_", " "), side=1, line=2, cex.lab=1, las=1, at = axis1[3])
  
  
  if(dim(coords)[2] > 1){
    axis(2, at= axis1) 
    mtext(stringr::str_replace_all(colnames(coords)[2], "_", " "), side=2, line=2, cex.lab=1, las=3, at = axis1[3])
  }
  
  for (i in 1:nrow(Y)) lines(rep(x[cw$edge[i, 2]], 
                                 2), Y[i, ], lwd = lwd[1], lend = 2)
  for (i in 1:cw$Nnode + n) {
    ii <- which(cw$edge[, 1] == i)
    if (length(ii) > 1) 
      lines(range(x[cw$edge[ii, 2]]), Y[ii, 1], lwd = lwd[1], 
            lend = 2)
  }
  for (i in 1:n) {
    if (ftype) 
      text(x[i], Y[which(cw$edge[, 2] == i), 2], 
           paste(" ", sub("_", " ", cw$tip.label[i]), 
                 sep = ""), pos = 4, offset = c(0, 1), srt = -90, 
           cex = fsize, font = ftype)
  }
}

# Function to calculate model performance

modelPerformanceCalculation <- function(miss, imp, true, type, trait){
  
  # r2 calculation
  matObs <- as.matrix(true[, trait])
  matNA <- as.matrix(miss[, trait])
  observed <- matObs[which(!is.na(matObs) & is.na(matNA))]
  predicted <- as.matrix(imp[, trait])
  predicted <- predicted[which(!is.na(matObs) & is.na(matNA))]
  
  # rmse calculation
  
  nrmse <- INDperform::nrmse(pred = predicted, obs = observed)
  
  # R2
  
  r2 <- summary(lm(observed ~ predicted))$adj.r.squared
  
  imp_performance <- data.frame("Type" = type,
                                "Variable" = trait,
                                "NRMSE" = as.numeric(round(nrmse, 3)),
                                "R2" = as.numeric(round(r2, 3)))
  
  return(imp_performance)
}

### GENERAL SPECIFICATIONS ----------------------------------------------------- ####

options(scipen = 999, digits = 3)

# Rerun models even if previous results exist. If TRUE, will overwrite previous results.
forceRunModels <- T

# Figures dimension

h <- 3
w <- 3

qgraphLayout <- "circle" # which kind of network layout do we want.
# qgraphLayout <- "spring"


# set the different paths for results using a subset of the data depending on the grouping variable

if(!exists("dataGroup")){
  stop("Define data set (traits or simulations)")
}

print(paste0("Running models and extracting results for ", dataGroup))

nClusters <- 2


#### DATA ---------------------------------------------------------------------- ####

if(dataGroup == "simulations"){
  
  set.seed(1)
  # tr <- rcoal(100)
  # plot(tr)
  tr <- pbtree(n = 100)
  
  # pdf("exampleTree.pdf", width = w, height = h)
  # plot(tr, show.tip.label = F)
  # dev.off()
  
  
  vcvMat <- matrix(c(1, 0.9, 0.8, 0, 0.1, 0.2,
                     0.9, 1, 0.8, 0, 0.1, 0.2,
                     0.8, 0.8, 1, 0, 0.1, 0.2,
                     0, 0, 0, 1, -0.9, -0.8,
                     0.1, 0.1, 0.1, -0.9, 1, 0.8,
                     0.2, 0.2, 0.2, -0.8, 0.8, 1), ncol = 6)
  
  simulatedDataOutput <- simulateDataSet(tr,vcvMatrix = vcvMat)
  
  if(!file.exists("results/results_simulations/tables/simulatedTraitsVcvMatrix.csv")){
    write.csv(simulatedDataOutput$vcvMatrix, "results/results_simulations/tables/simulatedTraitsVcvMatrix.csv", row.names = T)
    print("results/results_simulations/tables/simulatedTraitsVcvMatrix.csv")
  }
  
  df <- simulatedDataOutput$data
  names(df)
  
  
  
  traits <- c("phylo_G1_trait1", "phylo_G1_trait2", "phylo_G2_trait1", "phylo_G2_trait2", 
              "nonPhylo_G1_trait1", "nonPhylo_G1_trait2", "nonPhylo_G2_trait1", "nonPhylo_G2_trait2")
  predictors <- c("phylo_G1_env", "nonPhylo_G1_env")
  
  # traits <- traits[which(!is.na(str_extract(traits, "G1")))]
  # predictors <- predictors[which(!is.na(str_extract(predictors, "G1")))]
  
  vars <- c(traits, predictors)
  
  allNodeLabels <- data.frame(vars, str_replace_all(vars, "_", " "))
  traitsNodeLabels <- data.frame(traits, str_replace_all(traits, "_", " "))
}

if(dataGroup == "traits"){
  # traits data
  df <- read_csv("data_processed/traits_database_Pablo_2022.csv") 
  
  df <- df %>% rename(ln_abs_Pmin = ln_negMinWP_md, ln_abs_Ptlp = "ln_negPItlp")
  # names(df)
  # Species-level phylogeny
  tr <- read.tree("data_processed/species_phylogeny_Pablo_2022.nex")
  
  # traits <- c("ln_SLA", "ln_Aarea", "ln_N", "ln_LL", "ln_Kl", "ln_negPItlp") # include Kl and Ptlp?
  # predictors <- c("ln_AP",  "ln_MAT", "ln_AI", "ln_Prec_Seasonality", "ln_TSeasonality")
  
  traits <- c("ln_SLA", "ln_Aarea", "ln_N", "ln_LL")
  htraits <- c("ln_Kl", "ln_abs_Ptlp", "ln_abs_Pmin")
  
  traits <- c(traits, htraits)
  
  predictors <- c("ln_AI", "ln_Prec_Seasonality")
  vars <- c(traits, predictors)
  
  allNodeLabels <- data.frame(vars, str_replace_all(vars, "_", " "))
  traitsNodeLabels <- data.frame(traits, str_replace_all(traits, "_", " "))
  
}

# df %>% dplyr::filter(!is.na(LL)) %>% pull(Species) %>% length()

#### MODEL SPECIFICATIONS ------------------------------------------------------ ####

# uninformative priors
modelSpecifications <- defineModelsSpecifications(number.iterations = 100000, burning = 1000, thinning = 30,
                                                     # One response variable, phylogeny and residuals
                                                    uniresponse.prior = list(
                                                       R = list(V = 1, nu = 0.002), 
                                                       G = list(G1 = list(V = 1, nu = 0.002))
                                                     ),
                                                     # Two response variables, phylogeny and residuals
                                                  multiresponse.prior = list(
                                                       R=list(V=diag(2),nu=0.0005), 
                                                       G=list(G1=list(V=diag(2),nu=0.0005))
                                                     ))
                                                  

modelSpecifications <- defineModelsSpecifications(number.iterations = 1000000, burning = 10000, thinning = 20)

# modelSpecifications <- defineModelsSpecifications(number.iterations = 100, burning = 10, thinning = 5)


### IMPUTATIONS SPECIFICATIONS ------------------------------------------------- ####

clNumber <- ceiling(detectCores()/2)

# propNA <- 0.7
numberIterations <- 10
forceRunImputation <- T
nPhyloCoord <- 10
numberOfImpRounds <- 3
txtSize <- 0.8
