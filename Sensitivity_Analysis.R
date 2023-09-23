rm(list = ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readxl, gtsummary, dplyr, 
               tidyr, ggplot2, rotl, DescTools, stringr, ape, 
               emmeans, patchwork, latex2exp, metafor, brms, 
               flextable, phytools, MCMCglmm, metaAidR, orchaRd, 
               robumeta, ggpmisc, ggpubr)

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

##### Publication Bias #####

# Model of the Residuals from the Overall Model. 
Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Overall_Model.rds")
Residuals <- rstandard(Model)
Residuals[c("slab", "digits")] = NULL
Residuals_df <- do.call("rbind", Residuals)
Residuals_df2 <- t(Residuals_df)
colnames(Residuals_df2) <- c("yi", "vi", "z")

Residuals_Model <- rma(yi = yi, vi = vi, data = Residuals_df2)

# Publication Year Graph
Graph_Data <- data
Graph_Data <- Graph_Data %>% mutate(n_category = ifelse(n_R1.1 <= 50, "50", 
                                                 ifelse(n_R1.1 > 50 & n_R1.1 <= 100, "100", 
                                                 ifelse(n_R1.1 > 100 & n_R1.1 <= 150, "150", "> 150"))))


Publication_Graph <- ggplot(Graph_Data, aes(x = Year, y = InRR_Transformed)) + 
                     geom_point(aes(x = Year, y = InRR_Transformed, 
                                size = fct_relevel(n_category, c("50", "100", "150", "> 150"))), 
                                shape = 21, fill = "#4292c6", alpha = 0.5) + 
                     labs(x = "Publication Year", y = "Effect Size (InRR)", 
                          size = "Sample Size") +
                     theme_bw() +
                     theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                     theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                     theme(legend.position = "bottom", legend.direction = "horizontal") + 
                     geom_hline(yintercept = Model$b, lty = 2) + 
                     geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                     stat_poly_eq(formula = y ~ x, 
                                  aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                  parse = TRUE) #+
                     #coord_cartesian(xlim = c(0, 25), 
                     #                ylim = c(-5, 5))
Publication_Graph #(750x500)

# Funnel Plot
par(mfrow = c(1,2))
Funnel_Plot <- funnel(Residuals_Model, yaxis = "seinv",
                      ylab = "Inverse Standard Error (1/SE)", xlab = "Observed Outcome Residuals", 
                      pch = 21, back = "#D3DDEB", bg = "#183357")
box(lwd = 2)

# Trim and Fill
Trim_Fill <- trimfill(Residuals_Model, estimator = "R0")
Trim_Fill_Plot <- funnel(Trim_Fill, yaxis = "seinv", ylab = "",
                         xlab = "Observed Outcome Residuals", pch = 21, back = "#D3DDEB", bg = "#183357")
box(lwd = 2)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("top", legend = c("0.05 < p ≤ 1.00", "0 < p ≤ 0.05", "Studies", "Filled Studies"), 
       pch = c(22, 22, 21, 21), pt.bg = c("#FFFFFF","#D3DDEB", "#183357", "#FFFFFF"), box.lwd = 2)

# Egger's Regression Test (< 0.05 publication might be present)
Eggers_Test <- regtest(Residuals_Model, model = "lm")

# Time-lag Bias
run <- FALSE
system.time( #  36ish minutes
  if(run){
TL_Model <- metafor::rma.mv(InRR_Transformed, V = v_InRR, test = "t", dfs = "contain",
                            mods = ~ Year_Z + Precision - 1,
                            random = list(~1|phylo, ~1|Study_ID, ~1|obs,
                                          ~1|Shared_Animal_Number, ~1|Shared_Control_Number), 
                            R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                            control=list(rel.tol=1e-9))
saveRDS(TL_Model, "./3.Data_Analysis/2.Outputs/Models/TL_Model.rds")
  } else {
        TL_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/TL_Model.rds")})

##### Sensitivity Analysis #####

# Cooks Distance
run <- FALSE
system.time( #  9.5 Days
  if(run){
    Overall_Cooks <- cooks.distance(Model)
    saveRDS(Overall_Cooks, "./3.Data_Analysis/2.Outputs/Models/Overall_Cooks.rds")
  } else {
    Overall_Cooks <- readRDS("./3.Data_Analysis/2.Outputs/Models/Overall_Cooks.rds")})

dev.off()
Cooks_Plot <- plot(Overall_Cooks, type = "o", pch = 21, xlab = "Observed Outcome", 
                   ylab = "Cook's Distance", bg = "#183357")
box(lwd = 2)

# Untransformed and Hedges g Models
Model_Estimates <- data.frame(estimate = Model$b, ci.lb = Model$ci.lb, ci.ub = Model$ci.ub)
Model_i2 <- data.frame(round(orchaRd::i2_ml(Model), 2))

run <- FALSE
system.time( #  28ish minutes
  if(run){
    Untransformed_Model <- metafor::rma.mv(InRR_Untransformed ~ 1, V = v_InRR_Untransformed, test = "t", dfs = "contain",
                                           random = list(~1|phylo, ~1|Study_ID, ~1|obs,
                                                         ~1|Shared_Animal_Number, ~1|Shared_Control_Number), 
                                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                           control=list(rel.tol=1e-9))
    saveRDS(Untransformed_Model, "./3.Data_Analysis/2.Outputs/Models/Untransformed_Model.rds")
  } else {
    Untransformed_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/Untransformed_Model.rds")})

Untransformed_Model_Estimates <- data.frame(estimate = Untransformed_Model$b, ci.lb = Untransformed_Model$ci.lb, ci.ub = Untransformed_Model$ci.ub)
Untransformed_Model_i2 <- data.frame(round(orchaRd::i2_ml(Untransformed_Model), 2))

run <- FALSE
system.time( #  49ish minutes
  if(run){
    SMD_Model <- metafor::rma.mv(SMD ~ 1, V = v_SMD, test = "t", dfs = "contain",
                                 random = list(~1|phylo, ~1|Study_ID, ~1|obs,
                                               ~1|Shared_Animal_Number, ~1|Shared_Control_Number), 
                                 R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                 control=list(rel.tol=1e-9))
    saveRDS(SMD_Model, "./3.Data_Analysis/2.Outputs/Models/SMD_Model.rds")
  } else {
    SMD_Model <- readRDS("./3.Data_Analysis/2.Outputs/Models/SMD_Model.rds")})

SMD_Model_Estimates <- data.frame(estimate = SMD_Model$b, ci.lb = SMD_Model$ci.lb, ci.ub = SMD_Model$ci.ub)
SMD_Model_i2 <- data.frame(round(orchaRd::i2_ml(SMD_Model), 2))

##### Publication Bias (Log Coefficient of Variation Ratio) #####

# Model of the Residuals from the Overall Model. 
Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Overall_Model_CVR.rds")
Residuals_CVR <- rstandard(Model_CVR)
Residuals_CVR[c("slab", "digits")] = NULL
Residuals_df_CVR <- do.call("rbind", Residuals_CVR)
Residuals_df2_CVR <- t(Residuals_df_CVR)
colnames(Residuals_df2_CVR) <- c("yi", "vi", "z")

Residuals_Model_CVR <- rma(yi = yi, vi = vi, data = Residuals_df2_CVR)
par(mfrow = c(1,2))

# Publication Year Graph
Publication_CVR_Graph <- ggplot(Graph_Data, aes(x = Year, y = InCVR)) + 
                         geom_point(aes(x = Year, y = InCVR, 
                                    size = fct_relevel(n_category, c("50", "100", "150", "> 150"))), 
                                    shape = 21, fill = "#4292c6", alpha = 0.5) + 
                         labs(x = "Publication Year", y = "Effect Size (InCVR)", 
                              size = "Sample Size") +
                              theme_bw() +
                         theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                         theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                         theme(legend.position = "bottom", legend.direction = "horizontal") + 
                         geom_hline(yintercept = Model_CVR$b, lty = 2) + 
                         geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                         stat_poly_eq(formula = y ~ x, 
                                      aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                      parse = TRUE) #+
                         #coord_cartesian(xlim = c(0, 25), 
                         #                ylim = c(-5, 5))

Publication_CVR_Graph #(750x500)

# Funnel Plot
par(mfrow = c(1,2))
Funnel_Plot_CVR <- funnel(Residuals_Model_CVR, yaxis = "seinv",
                          ylab = "Inverse Standard Error (1/SE)", xlab = "Observed Outcome Residuals", 
                          pch = 21, back = "#D3DDEB", bg = "#183357")
box(lwd = 2)

# Trim and Fill
Trim_Fill_CVR <- trimfill(Residuals_Model_CVR, estimator = "R0")
Trim_Fill_Plot_CVR <- funnel(Trim_Fill_CVR, yaxis = "seinv", ylab = "",
                             xlab = "Observed Outcome Residuals", pch = 21, back = "#D3DDEB", bg = "#183357")
box(lwd = 2)

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("top", legend = c("0.05 < p ≤ 1.00", "0 < p ≤ 0.05", "Studies", "Filled Studies"), 
       pch = c(22, 22, 21, 21), pt.bg = c("#FFFFFF","#D3DDEB", "#183357", "#FFFFFF"), box.lwd = 2)

# Egger's Regression Test (< 0.05 publication might be present)
Eggers_Test_CVR <- regtest(Residuals_Model_CVR, model = "lm")

# Time-lag Bias
run <- TRUE
system.time( #  37ish minutes
  if(run){
    TL_Model_CVR <- metafor::rma.mv(InCVR, V = v_InCVR, test = "t", dfs = "contain",
                                          mods = ~ Year_Z + Precision_CVR - 1,
                                          random = list(~1|phylo, ~1|Study_ID, ~1|obs,
                                                        ~1|Shared_Animal_Number, ~1|Shared_Control_Number), 
                                          R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                          control=list(rel.tol=1e-9))
    saveRDS(TL_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/TL_Model_CVR.rds")
  } else {
            TL_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/TL_Model_CVR.rds")})

##### Sensitivity Analysis (Log Coefficient of Variation Ratio) #####

# Cooks Distance
run <- FALSE
system.time( #  9ish Days
  if(run){
    Overall_Cooks_CVR <- cooks.distance(Model_CVR)
    saveRDS(Overall_Cooks_CVR, "./3.Data_Analysis/2.Outputs/Models/Overall_Cooks_CVR.rds")
  } else {
    Overall_Cooks_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Overall_Cooks_CVR.rds")})

dev.off()
Cooks_CVR_Plot <- plot(Overall_Cooks_CVR, type = "o", pch = 21, xlab = "Observed Outcome", 
                   ylab = "Cook's Distance", bg = "#183357")
box(lwd = 2)

# Untransformed Model
Model_CVR_Estimates <- data.frame(estimate = Model_CVR$b, ci.lb = Model_CVR$ci.lb, ci.ub = Model_CVR$ci.ub)
Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Model_CVR), 2))

run <- FALSE
system.time( #  34ish minutes
  if(run){
    Untransformed_Model_CVR <- metafor::rma.mv(InCVR_Untransformed ~ 1, V = v_InCVR_Untransformed, test = "t", dfs = "contain",
                                           random = list(~1|phylo, ~1|Study_ID, ~1|obs,
                                                         ~1|Shared_Animal_Number, ~1|Shared_Control_Number), 
                                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                           control=list(rel.tol=1e-9))
    saveRDS(Untransformed_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Untransformed_Model_CVR.rds")
  } else {
            Untransformed_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Untransformed_Model_CVR.rds")})

Untransformed_Model_CVR_Estimates <- data.frame(estimate = Untransformed_Model_CVR$b, ci.lb = Untransformed_Model_CVR$ci.lb, ci.ub = Untransformed_Model_CVR$ci.ub)
Untransformed_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Untransformed_Model_CVR), 2))