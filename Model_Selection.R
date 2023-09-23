rm(list = ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readxl, gtsummary, dplyr, 
               tidyr, ggplot2, rotl, DescTools, stringr, ape, 
               emmeans, patchwork, latex2exp, metafor, brms, 
               flextable, phytools, MCMCglmm, metaAidR, orchaRd, robumeta)

##### Model Selection (Log Response Ratio) #####
# Importing Data Set
data <- read.csv("./3.Data_Analysis/2.Outputs/Data/Final_Data.csv")
data$obs <- 1:nrow(data)
data$Scientific_Name <- sub(" ", "_", data$Scientific_Name)
data$phylo <- data$Scientific_Name

# Phylogenetic covariance matrix
tree <- ape::read.tree("./3.Data_Analysis/2.Outputs/Phylogeny/tree")
phy <- ape::compute.brlen(tree, method = "Grafen", power = 1)
A <- ape::vcv.phylo(phy)
row.names(A) <- colnames(A) <- row.names(A)
A_cor <- ape::vcv.phylo(phy, corr = TRUE)

# All Possible Random Effects
run <- FALSE
system.time( #  80ish minutes
if(run){
  All <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                       ~1|Shared_Animal_Number, ~1|Shared_Control_Number, 
                                       ~1|Measurement), 
                         R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                         control=list(rel.tol=1e-9))
  saveRDS(All, "./3.Data_Analysis/2.Outputs/Models/All_Randoms.rds")
} else {
  All <- readRDS("./3.Data_Analysis/2.Outputs/Models/All_Randoms.rds")})

All_i2 <- data.frame(round(orchaRd::i2_ml(All), 2))
All_aic <- fitstats(All)

# No Measurement Random Effect
run <- FALSE
system.time( #  52ish minutes
  if(run){
    No_Measurement <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                    ~1|Shared_Animal_Number, ~1|Shared_Control_Number), 
                                      R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                      control=list(rel.tol=1e-9))
    saveRDS(No_Measurement, "./3.Data_Analysis/2.Outputs/Models/No_Measurement.rds")
  } else {
    No_Measurement <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement.rds")})

No_Measurement_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement), 2))
No_Measurement_aic <- fitstats(No_Measurement)

# No Shared Control Number Random Effect
run <- FALSE
system.time( #  47ish minutes
  if(run){
    No_Control <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                    ~1|Shared_Animal_Number, ~1|Measurement), 
                                      R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                      control=list(rel.tol=1e-9))
    saveRDS(No_Control, "./3.Data_Analysis/2.Outputs/Models/No_Control.rds")
  } else {
    No_Control <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Control.rds")})

No_Control_i2 <- data.frame(round(orchaRd::i2_ml(No_Control), 2))
No_Control_aic <- fitstats(No_Control)

# No Shared Animal Number Random Effect
run <- FALSE
system.time( #  45ish minutes
  if(run){
    No_Animal <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                               ~1|Shared_Control_Number, ~1|Measurement), 
                                 R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                 control=list(rel.tol=1e-9))
    saveRDS(No_Animal, "./3.Data_Analysis/2.Outputs/Models/No_Animal.rds")
  } else {
    No_Animal <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Animal.rds")})

No_Animal_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal), 2))
No_Animal_aic <- fitstats(No_Animal)

# No Species Random Effect
run <- FALSE
system.time( #  50ish minutes
  if(run){
    No_Species <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Shared_Animal_Number, 
                                               ~1|Shared_Control_Number, ~1|Measurement), 
                                 R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                 control=list(rel.tol=1e-9))
    saveRDS(No_Species, "./3.Data_Analysis/2.Outputs/Models/No_Species.rds")
  } else {
    No_Species <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Species.rds")})

No_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Species), 2))
No_Species_aic <- fitstats(No_Species)

# No Measurement or Shared Control Number Random Effects
run <- FALSE
system.time( #  35ish minutes
  if(run){
    No_Measurement_Control <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                         ~1|Shared_Animal_Number), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Control, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control.rds")
  } else {
    No_Measurement_Control <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control.rds")})

No_Measurement_Control_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Control), 2))
No_Measurement_Control_aic <- fitstats(No_Measurement_Control)

# No Measurement or Shared Animal Number Random Effects
run <- FALSE
system.time( #  38ish minutes
  if(run){
    No_Measurement_Animal <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                         ~1|Shared_Control_Number), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Animal, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Animal.rds")
  } else {
    No_Measurement_Animal <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Animal.rds")})

No_Measurement_Animal_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Animal), 2))
No_Measurement_Animal_aic <- fitstats(No_Measurement_Animal)

# No Measurement of Species Random Effects
run <- FALSE
system.time( #  34ish minutes
  if(run){
    No_Measurement_Species <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs,
                                         ~1|Shared_Animal_Number, ~1|Shared_Control_Number), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Species, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Species.rds")
  } else {
    No_Measurement_Species <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Species.rds")})

No_Measurement_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Species), 2))
No_Measurement_Species_aic <- fitstats(No_Measurement_Species)

# No Shared Control or Shared Animal Random Effects
run <- FALSE
system.time( #  34ish minutes
  if(run){
    No_Control_Animal <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                         ~1|Measurement), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(No_Control_Animal, "./3.Data_Analysis/2.Outputs/Models/No_Control_Animal.rds")
  } else {
    No_Control_Animal <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Control_Animal.rds")})

No_Control_Animal_i2 <- data.frame(round(orchaRd::i2_ml(No_Control_Animal), 2))
No_Control_Animal_aic <- fitstats(No_Control_Animal)

# No Shared Control or Species Random Effects
run <- FALSE
system.time( #  31ish minutes
  if(run){
    No_Control_Species <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                         ~1|Shared_Animal_Number, ~1|Measurement), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(No_Control_Species, "./3.Data_Analysis/2.Outputs/Models/No_Control_Species.rds")
  } else {
    No_Control_Species <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Control_Species.rds")})

No_Control_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Control_Species), 2))
No_Control_Species_aic <- fitstats(No_Control_Species)

# No Shared Animal or Species Random Effects
run <- FALSE
system.time( #  36ish minutes
  if(run){
    No_Animal_Species <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Shared_Control_Number, 
                                         ~1|Measurement), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(No_Animal_Species, "./3.Data_Analysis/2.Outputs/Models/No_Animal_Species.rds")
  } else {
    No_Animal_Species <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Animal_Species.rds")})

No_Animal_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal_Species), 2))
No_Animal_Species_aic <- fitstats(No_Animal_Species)

# No Measurement, Shared Control or Shared Animal Random Effects
run <- FALSE
system.time( #  27ish minutes
  if(run){
    No_Measurement_Control_Animal <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Control_Animal, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control_Animal.rds")
  } else {
    No_Measurement_Control_Animal <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control_Animal.rds")})

No_Measurement_Control_Animal_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Control_Animal), 2))
No_Measurement_Control_Animal_aic <- fitstats(No_Measurement_Control_Animal)

# No Measurement, Shared Control or Species Random Effects
run <- FALSE
system.time( #  28ish minutes
  if(run){
    No_Measurement_Control_Species <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Shared_Animal_Number), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Control_Species, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control_Species.rds")
  } else {
    No_Measurement_Control_Species <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control_Species.rds")})

No_Measurement_Control_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Control_Species), 2))
No_Measurement_Control_Species_aic <- fitstats(No_Measurement_Control_Species)

# No Measurement, Shared Animal or Species Random Effects
run <- FALSE
system.time( #  26ish minutes
  if(run){
    No_Measurement_Animal_Species <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Shared_Control_Number), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Animal_Species, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Animal_Species.rds")
  } else {
    No_Measurement_Animal_Species <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Animal_Species.rds")})

No_Measurement_Animal_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Animal_Species), 2))
No_Measurement_Animal_Species_aic <- fitstats(No_Measurement_Animal_Species)

# No Shared Control, Shared Animal or Species Random Effects
run <- FALSE
system.time( #  27ish minutes
  if(run){
    No_Control_Animal_Species <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Measurement), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(No_Control_Animal_Species, "./3.Data_Analysis/2.Outputs/Models/No_Control_Animal_Species.rds")
  } else {
    No_Control_Animal_Species <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Control_Animal_Species.rds")})

No_Control_Animal_Species_i2 <- data.frame(round(orchaRd::i2_ml(No_Control_Animal_Species), 2))
No_Control_Animal_Species_aic <- fitstats(No_Control_Animal_Species)

# The Base Model
run <- FALSE
system.time( #  20ish minutes
  if(run){
    Base <- metafor::rma.mv(InRR_Transformed ~ 1, V = v_InRR, test = "t", dfs = "contain",
                           random = list(~1|phylo, ~1|Study_ID, ~1|obs), 
                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                           control=list(rel.tol=1e-9))
    saveRDS(Base, "./3.Data_Analysis/2.Outputs/Models/Base.rds")
  } else {
    Base <- readRDS("./3.Data_Analysis/2.Outputs/Models/Base.rds")})

Base_i2 <- data.frame(round(orchaRd::i2_ml(Base), 2))
Base_aic <- fitstats(Base)

# AICc Summary
AICc <- data.frame("Models" = c("All", "No_Measurement", "No_Control", "No_Animal", "No_Species", 
                                "No_Measurement_Control", "No_Measurement_Animal", "No_Measurement_Species", 
                                "No_Control_Animal", "No_Control_Species", "No_Animal_Species", 
                                "No_Measurement_Control_Animal", "No_Measurement_Control_Species", 
                                "No_Measurement_Animal_Species", "No_Control_Animal_Species", "Base"), 
                   "AICc" = c(All_aic[5], No_Measurement_aic[5], No_Control_aic[5], No_Animal_aic[5], No_Species_aic[5], 
                              No_Measurement_Control_aic[5], No_Measurement_Animal_aic[5], No_Measurement_Species_aic[5], 
                              No_Control_Animal_aic[5], No_Control_Species_aic[5], No_Animal_Species_aic[5], 
                              No_Measurement_Control_Animal_aic[5], No_Measurement_Control_Species_aic[5], 
                              No_Measurement_Animal_Species_aic[5], No_Control_Animal_Species_aic[5], Base_aic[5]))

# Heterogeneity Summary
Heterogeneity <- data.frame("Random Effect" = c("Phylogeny", "Scientific_Name", "Shared_Animal_Number", "Shared_Control_Number", "Study_ID", "Measurement", "Obs", "Total"), 
                            "All" = c(All_i2["I2_phylo", 1], All_i2["I2_Scientific_Name", 1], All_i2["I2_Shared_Animal_Number", 1], All_i2["I2_Shared_Control_Number", 1], 
                                      All_i2["I2_Study_ID", 1], All_i2["I2_Measurement", 1], All_i2["I2_obs", 1], All_i2["I2_Total", 1]), 
                            "No_Measurement" = c(No_Measurement_i2["I2_phylo", 1], No_Measurement_i2["I2_Scientific_Name", 1], No_Measurement_i2["I2_Shared_Animal_Number", 1], 
                                                 No_Measurement_i2["I2_Shared_Control_Number", 1], No_Measurement_i2["I2_Study_ID", 1], NA, 
                                                 No_Measurement_i2["I2_obs", 1], No_Measurement_i2["I2_Total", 1]),
                            "No_Control" = c(No_Control_i2["I2_phylo", 1], No_Control_i2["I2_Scientific_Name", 1], No_Control_i2["I2_Shared_Animal_Number", 1], 
                                             NA, No_Control_i2["I2_Study_ID", 1], No_Control_i2["I2_Measurement", 1], 
                                             No_Control_i2["I2_obs", 1], No_Control_i2["I2_Total", 1]),
                            "No_Animal" = c(No_Animal_i2["I2_phylo", 1], No_Animal_i2["I2_Scientific_Name", 1], NA, 
                                            No_Animal_i2["I2_Shared_Control_Number", 1], No_Animal_i2["I2_Study_ID", 1], No_Animal_i2["I2_Measurement", 1], 
                                            No_Animal_i2["I2_obs", 1], No_Animal_i2["I2_Total", 1]),
                            "No_Species" = c(No_Species_i2["I2_phylo", 1], NA, No_Species_i2["I2_Shared_Animal_Number", 1], 
                                             No_Species_i2["I2_Shared_Control_Number", 1], No_Species_i2["I2_Study_ID", 1], No_Species_i2["I2_Measurement", 1], 
                                             No_Species_i2["I2_obs", 1], No_Species_i2["I2_Total", 1]),
                            "No_Measurement_Control" = c(No_Measurement_Control_i2["I2_phylo", 1], No_Measurement_Control_i2["I2_Scientific_Name", 1], 
                                                         No_Measurement_Control_i2["I2_Shared_Animal_Number", 1], NA, 
                                                         No_Measurement_Control_i2["I2_Study_ID", 1], NA, 
                                                         No_Measurement_Control_i2["I2_obs", 1], No_Measurement_Control_i2["I2_Total", 1]), 
                            "No_Measurement_Animal" = c(No_Measurement_Animal_i2["I2_phylo", 1], No_Measurement_Animal_i2["I2_Scientific_Name", 1], 
                                                        NA, No_Measurement_Animal_i2["I2_Shared_Control_Number", 1], 
                                                        No_Measurement_Animal_i2["I2_Study_ID", 1], NA, 
                                                        No_Measurement_Animal_i2["I2_obs", 1], No_Measurement_Animal_i2["I2_Total", 1]), 
                            "No_Measurement_Species" = c(No_Measurement_Species_i2["I2_phylo", 1], NA, 
                                                         No_Measurement_Species_i2["I2_Shared_Animal_Number", 1], No_Measurement_Species_i2["I2_Shared_Control_Number", 1], 
                                                         No_Measurement_Species_i2["I2_Study_ID", 1], NA, 
                                                         No_Measurement_Species_i2["I2_obs", 1], No_Measurement_Species_i2["I2_Total", 1]), 
                            "No_Control_Animal" = c(No_Control_Animal_i2["I2_phylo", 1], No_Control_Animal_i2["I2_Scientific_Name", 1], 
                                                    NA, NA, 
                                                    No_Control_Animal_i2["I2_Study_ID", 1], No_Control_Animal_i2["I2_Measurement", 1], 
                                                    No_Control_Animal_i2["I2_obs", 1], No_Control_Animal_i2["I2_Total", 1]),
                            "No_Control_Species" = c(No_Control_Species_i2["I2_phylo", 1], NA, 
                                                     No_Control_Species_i2["I2_Shared_Animal_Number", 1], NA, 
                                                     No_Control_Species_i2["I2_Study_ID", 1], No_Control_Species_i2["I2_Measurement", 1], 
                                                     No_Control_Species_i2["I2_obs", 1], No_Control_Species_i2["I2_Total", 1]),
                            "No_Animal_Species" = c(No_Animal_Species_i2["I2_phylo", 1], NA, 
                                                    NA, No_Animal_Species_i2["I2_Shared_Control_Number", 1], 
                                                    No_Animal_Species_i2["I2_Study_ID", 1], No_Animal_Species_i2["I2_Measurement", 1], 
                                                    No_Animal_Species_i2["I2_obs", 1], No_Animal_Species_i2["I2_Total", 1]),
                            "No_Measurement_Control_Animal" = c(No_Measurement_Control_Animal_i2["I2_phylo", 1], No_Measurement_Control_Animal_i2["I2_Scientific_Name", 1], 
                                                                NA, NA, 
                                                                No_Measurement_Control_Animal_i2["I2_Study_ID", 1], NA, 
                                                                No_Measurement_Control_Animal_i2["I2_obs", 1], No_Measurement_Control_Animal_i2["I2_Total", 1]),
                            "No_Measurement_Control_Species" = c(No_Measurement_Control_Species_i2["I2_phylo", 1], NA, 
                                                                 No_Measurement_Control_Species_i2["I2_Shared_Animal_Number", 1], NA, 
                                                                 No_Measurement_Control_Species_i2["I2_Study_ID", 1], NA, 
                                                                 No_Measurement_Control_Species_i2["I2_obs", 1], No_Measurement_Control_Species_i2["I2_Total", 1]),
                            "No_Measurement_Animal_Species" = c(No_Measurement_Animal_Species_i2["I2_phylo", 1], NA, 
                                                                NA, No_Measurement_Animal_Species_i2["I2_Shared_Control_Number", 1], 
                                                                No_Measurement_Animal_Species_i2["I2_Study_ID", 1], NA, 
                                                                No_Measurement_Animal_Species_i2["I2_obs", 1], No_Measurement_Animal_Species_i2["I2_Total", 1]),
                            "No_Control_Animal_Species" = c(No_Control_Animal_Species_i2["I2_phylo", 1], NA,
                                                            NA, NA, 
                                                            No_Control_Animal_Species_i2["I2_Study_ID", 1], No_Control_Animal_Species_i2["I2_Measurement", 1], 
                                                            No_Control_Animal_Species_i2["I2_obs", 1], No_Control_Animal_Species_i2["I2_Total", 1]),
                            "Base" = c(Base_i2["I2_phylo", 1], NA, NA, NA, 
                                       Base_i2["I2_Study_ID", 1], NA, Base_i2["I2_obs", 1], Base_i2["I2_Total", 1]))


##### Model Selection (Log Coefficient of Variation Ratio) #####
# All Possible Random Effects
run <- FALSE
system.time( #  44ish minutes
  if(run){
    All_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                               random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                             ~1|Shared_Animal_Number, ~1|Shared_Control_Number, 
                                             ~1|Measurement), 
                               R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                               control=list(rel.tol=1e-9))
    saveRDS(All_CVR, "./3.Data_Analysis/2.Outputs/Models/All_CVR.rds")
  } else {
    All_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/All_CVR.rds")})

All_CVR_i2 <- data.frame(round(orchaRd::i2_ml(All_CVR), 2))
All_CVR_aic <- fitstats(All_CVR)

# No Measurement Random Effect
run <- FALSE
system.time( #  36ish minutes
  if(run){
    No_Measurement_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                    ~1|Shared_Animal_Number, ~1|Shared_Control_Number), 
                                      R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                      control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_CVR.rds")
  } else {
    No_Measurement_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_CVR.rds")})

No_Measurement_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_CVR), 2))
No_Measurement_CVR_aic <- fitstats(No_Measurement_CVR)

# No Shared Control Number Random Effect
run <- FALSE
system.time( #  36ish minutes
  if(run){
    No_Control_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                  random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                ~1|Shared_Animal_Number, ~1|Measurement), 
                                  R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                  control=list(rel.tol=1e-9))
    saveRDS(No_Control_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Control_CVR.rds")
  } else {
    No_Control_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Control_CVR.rds")})

No_Control_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Control_CVR), 2))
No_Control_CVR_aic <- fitstats(No_Control_CVR)

# No Shared Animal Number Random Effect
run <- FALSE
system.time( #  38ish minutes
  if(run){
    No_Animal_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                               ~1|Shared_Control_Number, ~1|Measurement), 
                                 R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                 control=list(rel.tol=1e-9))
    saveRDS(No_Animal_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Animal_CVR.rds")
  } else {
    No_Animal_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Animal_CVR.rds")})

No_Animal_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal_CVR), 2))
No_Animal_CVR_aic <- fitstats(No_Animal_CVR)

# No Species Random Effect
run <- FALSE
system.time( #  38ish minutes
  if(run){
    No_Species_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                  random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Shared_Animal_Number, 
                                                ~1|Shared_Control_Number, ~1|Measurement), 
                                  R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                  control=list(rel.tol=1e-9))
    saveRDS(No_Species_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Species_CVR.rds")
  } else {
    No_Species_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Species_CVR.rds")})

No_Species_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Species_CVR), 2))
No_Species_CVR_aic <- fitstats(No_Species_CVR)

# No Measurement or Shared Control Number Random Effects
run <- FALSE
system.time( #  32ish minutes
  if(run){
    No_Measurement_Control_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                              random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                            ~1|Shared_Animal_Number), 
                                              R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                              control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Control_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control_CVR.rds")
  } else {
    No_Measurement_Control_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control_CVR.rds")})

No_Measurement_Control_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Control_CVR), 2))
No_Measurement_Control_CVR_aic <- fitstats(No_Measurement_Control_CVR)

# No Measurement or Shared Animal Number Random Effects
run <- FALSE
system.time( #  33ish minutes
  if(run){
    No_Measurement_Animal_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                             random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                           ~1|Shared_Control_Number), 
                                             R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                             control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Animal_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Animal_CVR.rds")
  } else {
    No_Measurement_Animal_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Animal_CVR.rds")})

No_Measurement_Animal_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Animal_CVR), 2))
No_Measurement_Animal_CVR_aic <- fitstats(No_Measurement_Animal_CVR)

# No Measurement of Species Random Effects
run <- FALSE
system.time( #  32ish minutes
  if(run){
    No_Measurement_Species_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                              random = list(~1|phylo, ~1|Study_ID, ~1|obs,
                                                            ~1|Shared_Animal_Number, ~1|Shared_Control_Number), 
                                              R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                              control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Species_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Species_CVR.rds")
  } else {
    No_Measurement_Species_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Species_CVR.rds")})

No_Measurement_Species_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Species_CVR), 2))
No_Measurement_Species_CVR_aic <- fitstats(No_Measurement_Species_CVR)

# No Shared Control or Shared Animal Random Effects
run <- FALSE
system.time( #  34ish minutes
  if(run){
    No_Control_Animal_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                       ~1|Measurement), 
                                         R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                         control=list(rel.tol=1e-9))
    saveRDS(No_Control_Animal_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Control_Animal_CVR.rds")
  } else {
    No_Control_Animal_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Control_Animal_CVR.rds")})

No_Control_Animal_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Control_Animal_CVR), 2))
No_Control_Animal_CVR_aic <- fitstats(No_Control_Animal_CVR)

# No Shared Control or Species Random Effects
run <- FALSE
system.time( #  34ish minutes
  if(run){
    No_Control_Species_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                          random = list(~1|phylo, ~1|Study_ID, ~1|obs, 
                                                        ~1|Shared_Animal_Number, ~1|Measurement), 
                                          R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                          control=list(rel.tol=1e-9))
    saveRDS(No_Control_Species_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Control_Species_CVR.rds")
  } else {
    No_Control_Species_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Control_Species_CVR.rds")})

No_Control_Species_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Control_Species_CVR), 2))
No_Control_Species_CVR_aic <- fitstats(No_Control_Species_CVR)

# No Shared Animal or Species Random Effects
run <- FALSE
system.time( #  34ish minutes
  if(run){
    No_Animal_Species_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Shared_Control_Number, 
                                                       ~1|Measurement), 
                                         R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                         control=list(rel.tol=1e-9))
    saveRDS(No_Animal_Species_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Animal_Species_CVR.rds")
  } else {
    No_Animal_Species_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Animal_Species_CVR.rds")})

No_Animal_Species_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Animal_Species_CVR), 2))
No_Animal_Species_CVR_aic <- fitstats(No_Animal_Species_CVR)

# No Measurement, Shared Control or Shared Animal Random Effects
run <- FALSE
system.time( #  26ish minutes
  if(run){
    No_Measurement_Control_Animal_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name), 
                                                     R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                                     control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Control_Animal_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control_Animal_CVR.rds")
  } else {
    No_Measurement_Control_Animal_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control_Animal_CVR.rds")})

No_Measurement_Control_Animal_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Control_Animal_CVR), 2))
No_Measurement_Control_Animal_CVR_aic <- fitstats(No_Measurement_Control_Animal_CVR)

# No Measurement, Shared Control or Species Random Effects
run <- FALSE
system.time( #  26ish minutes
  if(run){
    No_Measurement_Control_Species_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Shared_Animal_Number), 
                                                      R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                                      control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Control_Species_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control_Species_CVR.rds")
  } else {
    No_Measurement_Control_Species_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Control_Species_CVR.rds")})

No_Measurement_Control_Species_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Control_Species_CVR), 2))
No_Measurement_Control_Species_CVR_aic <- fitstats(No_Measurement_Control_Species_CVR)

# No Measurement, Shared Animal or Species Random Effects
run <- FALSE
system.time( #  28ish minutes
  if(run){
    No_Measurement_Animal_Species_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Shared_Control_Number), 
                                                     R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                                     control=list(rel.tol=1e-9))
    saveRDS(No_Measurement_Animal_Species_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Measurement_Animal_Species_CVR.rds")
  } else {
    No_Measurement_Animal_Species_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Measurement_Animal_Species_CVR.rds")})

No_Measurement_Animal_Species_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Measurement_Animal_Species_CVR), 2))
No_Measurement_Animal_Species_CVR_aic <- fitstats(No_Measurement_Animal_Species_CVR)

# No Shared Control, Shared Animal or Species Random Effects
run <- FALSE
system.time( #  26ish minutes
  if(run){
    No_Control_Animal_Species_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Measurement), 
                                                 R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                                 control=list(rel.tol=1e-9))
    saveRDS(No_Control_Animal_Species_CVR, "./3.Data_Analysis/2.Outputs/Models/No_Control_Animal_Species_CVR.rds")
  } else {
    No_Control_Animal_Species_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/No_Control_Animal_Species_CVR.rds")})

No_Control_Animal_Species_CVR_i2 <- data.frame(round(orchaRd::i2_ml(No_Control_Animal_Species_CVR), 2))
No_Control_Animal_Species_CVR_aic <- fitstats(No_Control_Animal_Species_CVR)

# The Base Model
run <- FALSE
system.time( #  19ish minutes
  if(run){
    Base_CVR <- metafor::rma.mv(InCVR ~ 1, V = v_InCVR, test = "t", dfs = "contain",
                            random = list(~1|phylo, ~1|Study_ID, ~1|obs), 
                            R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                            control=list(rel.tol=1e-9))
    saveRDS(Base_CVR, "./3.Data_Analysis/2.Outputs/Models/Base_CVR.rds")
  } else {
    Base_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Base_CVR.rds")})

Base_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Base_CVR), 2))
Base_CVR_aic <- fitstats(Base_CVR)

# AICc Summary
AICc_CVR <- data.frame("Models" = c("All_CVR", "No_Measurement_CVR", "No_Control_CVR", "No_Animal_CVR", "No_Species_CVR", 
                                "No_Measurement_Control_CVR", "No_Measurement_Animal_CVR", "No_Measurement_Species_CVR", 
                                "No_Control_Animal_CVR", "No_Control_Species_CVR", "No_Animal_Species_CVR", 
                                "No_Measurement_Control_Animal_CVR", "No_Measurement_Control_Species_CVR", 
                                "No_Measurement_Animal_Species_CVR", "No_Control_Animal_Species_CVR", "Base_CVR"), 
                   "AICc" = c(All_CVR_aic[5], No_Measurement_CVR_aic[5], No_Control_CVR_aic[5], No_Animal_CVR_aic[5], No_Species_CVR_aic[5], 
                              No_Measurement_Control_CVR_aic[5], No_Measurement_Animal_CVR_aic[5], No_Measurement_Species_CVR_aic[5], 
                              No_Control_Animal_CVR_aic[5], No_Control_Species_CVR_aic[5], No_Animal_Species_CVR_aic[5], 
                              No_Measurement_Control_Animal_CVR_aic[5], No_Measurement_Control_Species_CVR_aic[5], 
                              No_Measurement_Animal_Species_CVR_aic[5], No_Control_Animal_Species_CVR_aic[5], Base_CVR_aic[5]))

# Heterogeneity Summary
Heterogeneity_CVR <- data.frame("Random Effect" = c("Phylogeny", "Scientific_Name", "Shared_Animal_Number", "Shared_Control_Number", "Study_ID", "Measurement", "Obs", "Total"), 
                            "All" = c(All_CVR_i2["I2_phylo", 1], All_CVR_i2["I2_Scientific_Name", 1], All_CVR_i2["I2_Shared_Animal_Number", 1], All_CVR_i2["I2_Shared_Control_Number", 1], 
                                      All_CVR_i2["I2_Study_ID", 1], All_CVR_i2["I2_Measurement", 1], All_CVR_i2["I2_obs", 1], All_CVR_i2["I2_Total", 1]), 
                            "No_Measurement" = c(No_Measurement_CVR_i2["I2_phylo", 1], No_Measurement_CVR_i2["I2_Scientific_Name", 1], No_Measurement_CVR_i2["I2_Shared_Animal_Number", 1], 
                                                 No_Measurement_CVR_i2["I2_Shared_Control_Number", 1], No_Measurement_CVR_i2["I2_Study_ID", 1], NA, 
                                                 No_Measurement_CVR_i2["I2_obs", 1], No_Measurement_CVR_i2["I2_Total", 1]),
                            "No_Control" = c(No_Control_CVR_i2["I2_phylo", 1], No_Control_CVR_i2["I2_Scientific_Name", 1], No_Control_CVR_i2["I2_Shared_Animal_Number", 1], 
                                             NA, No_Control_CVR_i2["I2_Study_ID", 1], No_Control_CVR_i2["I2_Measurement", 1], 
                                             No_Control_CVR_i2["I2_obs", 1], No_Control_CVR_i2["I2_Total", 1]),
                            "No_Animal" = c(No_Animal_CVR_i2["I2_phylo", 1], No_Animal_CVR_i2["I2_Scientific_Name", 1], NA, 
                                            No_Animal_CVR_i2["I2_Shared_Control_Number", 1], No_Animal_CVR_i2["I2_Study_ID", 1], No_Animal_CVR_i2["I2_Measurement", 1], 
                                            No_Animal_CVR_i2["I2_obs", 1], No_Animal_CVR_i2["I2_Total", 1]),
                            "No_Species" = c(No_Species_CVR_i2["I2_phylo", 1], NA, No_Species_CVR_i2["I2_Shared_Animal_Number", 1], 
                                             No_Species_CVR_i2["I2_Shared_Control_Number", 1], No_Species_CVR_i2["I2_Study_ID", 1], No_Species_CVR_i2["I2_Measurement", 1], 
                                             No_Species_CVR_i2["I2_obs", 1], No_Species_CVR_i2["I2_Total", 1]),
                            "No_Measurement_Control" = c(No_Measurement_Control_CVR_i2["I2_phylo", 1], No_Measurement_Control_CVR_i2["I2_Scientific_Name", 1], 
                                                         No_Measurement_Control_CVR_i2["I2_Shared_Animal_Number", 1], NA, 
                                                         No_Measurement_Control_CVR_i2["I2_Study_ID", 1], NA, 
                                                         No_Measurement_Control_CVR_i2["I2_obs", 1], No_Measurement_Control_CVR_i2["I2_Total", 1]), 
                            "No_Measurement_Animal" = c(No_Measurement_Animal_CVR_i2["I2_phylo", 1], No_Measurement_Animal_CVR_i2["I2_Scientific_Name", 1], 
                                                        NA, No_Measurement_Animal_CVR_i2["I2_Shared_Control_Number", 1], 
                                                        No_Measurement_Animal_CVR_i2["I2_Study_ID", 1], NA, 
                                                        No_Measurement_Animal_CVR_i2["I2_obs", 1], No_Measurement_Animal_CVR_i2["I2_Total", 1]), 
                            "No_Measurement_Species" = c(No_Measurement_Species_CVR_i2["I2_phylo", 1], NA, 
                                                         No_Measurement_Species_CVR_i2["I2_Shared_Animal_Number", 1], No_Measurement_Species_CVR_i2["I2_Shared_Control_Number", 1], 
                                                         No_Measurement_Species_CVR_i2["I2_Study_ID", 1], NA, 
                                                         No_Measurement_Species_CVR_i2["I2_obs", 1], No_Measurement_Species_CVR_i2["I2_Total", 1]), 
                            "No_Control_Animal" = c(No_Control_Animal_CVR_i2["I2_phylo", 1], No_Control_Animal_CVR_i2["I2_Scientific_Name", 1], 
                                                    NA, NA, 
                                                    No_Control_Animal_CVR_i2["I2_Study_ID", 1], No_Control_Animal_CVR_i2["I2_Measurement", 1], 
                                                    No_Control_Animal_CVR_i2["I2_obs", 1], No_Control_Animal_CVR_i2["I2_Total", 1]),
                            "No_Control_Species" = c(No_Control_Species_CVR_i2["I2_phylo", 1], NA, 
                                                     No_Control_Species_CVR_i2["I2_Shared_Animal_Number", 1], NA, 
                                                     No_Control_Species_CVR_i2["I2_Study_ID", 1], No_Control_Species_CVR_i2["I2_Measurement", 1], 
                                                     No_Control_Species_CVR_i2["I2_obs", 1], No_Control_Species_CVR_i2["I2_Total", 1]),
                            "No_Animal_Species" = c(No_Animal_Species_CVR_i2["I2_phylo", 1], NA, 
                                                    NA, No_Animal_Species_CVR_i2["I2_Shared_Control_Number", 1], 
                                                    No_Animal_Species_CVR_i2["I2_Study_ID", 1], No_Animal_Species_CVR_i2["I2_Measurement", 1], 
                                                    No_Animal_Species_CVR_i2["I2_obs", 1], No_Animal_Species_CVR_i2["I2_Total", 1]),
                            "No_Measurement_Control_Animal" = c(No_Measurement_Control_Animal_CVR_i2["I2_phylo", 1], No_Measurement_Control_Animal_CVR_i2["I2_Scientific_Name", 1], 
                                                                NA, NA, 
                                                                No_Measurement_Control_Animal_CVR_i2["I2_Study_ID", 1], NA, 
                                                                No_Measurement_Control_Animal_CVR_i2["I2_obs", 1], No_Measurement_Control_Animal_CVR_i2["I2_Total", 1]),
                            "No_Measurement_Control_Species" = c(No_Measurement_Control_Species_CVR_i2["I2_phylo", 1], NA, 
                                                                 No_Measurement_Control_Species_CVR_i2["I2_Shared_Animal_Number", 1], NA, 
                                                                 No_Measurement_Control_Species_CVR_i2["I2_Study_ID", 1], NA, 
                                                                 No_Measurement_Control_Species_CVR_i2["I2_obs", 1], No_Measurement_Control_Species_CVR_i2["I2_Total", 1]),
                            "No_Measurement_Animal_Species" = c(No_Measurement_Animal_Species_CVR_i2["I2_phylo", 1], NA, 
                                                                NA, No_Measurement_Animal_Species_CVR_i2["I2_Shared_Control_Number", 1], 
                                                                No_Measurement_Animal_Species_CVR_i2["I2_Study_ID", 1], NA, 
                                                                No_Measurement_Animal_Species_CVR_i2["I2_obs", 1], No_Measurement_Animal_Species_CVR_i2["I2_Total", 1]),
                            "No_Control_Animal_Species" = c(No_Control_Animal_Species_CVR_i2["I2_phylo", 1], NA,
                                                            NA, NA, 
                                                            No_Control_Animal_Species_CVR_i2["I2_Study_ID", 1], No_Control_Animal_Species_CVR_i2["I2_Measurement", 1], 
                                                            No_Control_Animal_Species_CVR_i2["I2_obs", 1], No_Control_Animal_Species_CVR_i2["I2_Total", 1]),
                            "Base" = c(Base_CVR_i2["I2_phylo", 1], NA, NA, NA, 
                                       Base_CVR_i2["I2_Study_ID", 1], NA, Base_CVR_i2["I2_obs", 1], Base_CVR_i2["I2_Total", 1]))