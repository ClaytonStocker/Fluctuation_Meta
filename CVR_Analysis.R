rm(list = ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readxl, gtsummary, dplyr, 
               tidyr, ggplot2, rotl, DescTools, stringr, ape, 
               emmeans, patchwork, latex2exp, metafor, brms, 
               flextable, phytools, MCMCglmm, metaAidR, orchaRd, 
               robumeta, ggpmisc, ggridges, ggbeeswarm, gridExtra, phylobase)

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

# Variance matrix (Shared Control)

VCV_InCVR <- make_VCV_matrix(data, V = "v_InCVR", cluster = "Shared_Control_Number")

##### Overall Model - CVR #####
run <- FALSE
system.time( #  34ish minutes
  if(run){
    Overall_Model_CVR <- metafor::rma.mv(InCVR ~ 1, V = VCV_InCVR, test = "t", dfs = "contain",
                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                       ~1|Shared_Animal_Number, ~1|Measurement), 
                                         R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                         control=list(rel.tol=1e-9))
    saveRDS(Overall_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Overall_Model_CVR.rds")
  } else {
            Overall_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Overall_Model_CVR.rds")})

Overall_Model_CVR_rob <- robust(Overall_Model_CVR, cluster = data$Study_ID, adjust = TRUE)

Overall_Model_CVR_Estimates <- data.frame(estimate = Overall_Model_CVR$b, ci.lb = Overall_Model_CVR$ci.lb, 
                                          ci.ub = Overall_Model_CVR$ci.ub)
Overall_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Overall_Model_CVR), 2))

#### Overall Model - Fluctuation Amplitude Meta-Regression - CVR ####
run <- FALSE
system.time( #  14ish minutes
  if(run){
    Amplitude_Model_CVR <- metafor::rma.mv(InCVR, V = VCV_InCVR, test = "t", dfs = "contain",
                                           mods = ~ T2_Magnitude - 1,
                                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                         ~1|Shared_Animal_Number, ~1|Measurement), 
                                           R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                           control=list(rel.tol=1e-9))
    saveRDS(Amplitude_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Amplitude_Model_CVR.rds")
  } else {
            Amplitude_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Amplitude_Model_CVR.rds")})

Amplitude_Model_CVR_rob <- robust(Amplitude_Model_CVR, cluster = data$Study_ID, adjust = TRUE)

Amplitude_Model_CVR_Estimates <- data.frame(estimate = Amplitude_Model_CVR$b, ci.lb = Amplitude_Model_CVR$ci.lb, 
                                            ci.ub = Amplitude_Model_CVR$ci.ub)
Amplitude_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Amplitude_Model_CVR), 2))

Plot_Data <- data
Plot_Data <- Plot_Data %>% mutate(n_category = ifelse(n_R1.1 <= 25, "25", 
                                               ifelse(n_R1.1 > 25 & n_R1.1 <= 50, "50", 
                                               ifelse(n_R1.1 > 50 & n_R1.1 <= 75, "75", "> 75"))))

Amplitude_Plot_CVR <- ggplot(Plot_Data, aes(x = T2_Magnitude, y = InCVR)) + 
                      geom_point(aes(x = T2_Magnitude, y = InCVR, 
                                 size = fct_relevel(n_category, c("25", "50", "75", "> 75"))), 
                                 shape = 21, fill = "#4292c6", alpha = 0.5) + 
                      labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (lnCVR)", 
                           size = "Sample Size", title = "Overall Analysis") +
                      theme_bw() +
                      theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                      theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                      theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                      theme(legend.position = "bottom", legend.direction = "horizontal") + 
                      geom_hline(yintercept = Overall_Model_CVR_Estimates$estimate, lty = 2) + 
                      geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                      stat_poly_eq(formula = y ~ x, 
                                   aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                   parse = TRUE) +
                      coord_cartesian(xlim = c(0, 30), 
                                      ylim = c(-2.5, 2.5))

Amplitude_Plot_CVR #(400x400)

#### Overall Model - Type of Fluctuation Meta-Regression - CVR ####
Fluctuation_Data <- data %>% filter(!is.na(Fluctuation_Category))

Fluctuation_Exploration <- Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Fluctuation_Exploration) <- Fluctuation_Exploration$Fluctuation_Category

Fluctuation_Species_Count <- Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% 
                             table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                             select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Fluctuation_Species_Count) <- Fluctuation_Species_Count$Fluctuation_Category

Fluctuation_Study_Count <- Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% 
                           table() %>% data.frame() %>% filter(`Freq` != 0) %>% 
                           select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Fluctuation_Study_Count) <- Fluctuation_Study_Count$Fluctuation_Category

Fluctuation_Species <- Fluctuation_Data %>% select("phylo") %>% unique()

Fluctuation_A_cor <- as.data.frame(A_cor)
Fluctuation_A_cor <- Fluctuation_A_cor[c(Fluctuation_Species$phylo), c(Fluctuation_Species$phylo)]
Fluctuation_A_cor <- as.matrix(Fluctuation_A_cor)

Fluctuation_VCV_InCVR <- make_VCV_matrix(Fluctuation_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  29ish minutes
  if(run){
    Fluctuation_Model_CVR <- metafor::rma.mv(InCVR, V = Fluctuation_VCV_InCVR, test = "t", dfs = "contain",
                                             mods = ~ Fluctuation_Category - 1,
                                             random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                           ~1|Shared_Animal_Number, ~1|Measurement), 
                                             R = list(phylo=Fluctuation_A_cor), data = Fluctuation_Data, method = "REML", sparse = TRUE, 
                                             control=list(rel.tol=1e-9))
    saveRDS(Fluctuation_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Fluctuation_Model_CVR.rds")
  } else {
            Fluctuation_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Fluctuation_Model_CVR.rds")})

Fluctuation_Model_CVR_rob <- robust(Fluctuation_Model_CVR, cluster = Fluctuation_Data$Study_ID, adjust = TRUE)

Fluctuation_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Fluctuation_Model_CVR$b), 21, 100),
                                              estimate = Fluctuation_Model_CVR$b, ci.lb = Fluctuation_Model_CVR$ci.lb, 
                                              ci.ub = Fluctuation_Model_CVR$ci.ub)
rownames(Fluctuation_Model_CVR_Estimates) <- Fluctuation_Model_CVR_Estimates$Category
Fluctuation_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Fluctuation_Model_CVR), 2))


# Preparing Graph - Combined

fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise", "Stochastic")

fluctuation_k <- data.frame("k" = c(Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                    Fluctuation_Exploration["Alternating", "Freq"], 
                                    Fluctuation_Exploration["Stepwise", "Freq"], 
                                    Fluctuation_Exploration["Stochastic", "Freq"]), 
                            row.names = fluctuation_rnames)

fluctuation_group_no <- data.frame("Spp No." = c(Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                 Fluctuation_Species_Count["Alternating", "Freq"], 
                                                 Fluctuation_Species_Count["Stepwise", "Freq"], 
                                                 Fluctuation_Species_Count["Stochastic", "Freq"]), 
                                   row.names = fluctuation_rnames)

fluctuation_study <- data.frame("Study" = c(Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                            Fluctuation_Study_Count["Alternating", "Freq"], 
                                            Fluctuation_Study_Count["Stepwise", "Freq"], 
                                            Fluctuation_Study_Count["Stochastic", "Freq"]), 
                                   row.names = fluctuation_rnames)

Fluctuation_Model_CVR_Estimates_Reorder <- Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise", "Stochastic"), ]

fluctuation_table <- data.frame(estimate = Fluctuation_Model_CVR_Estimates_Reorder[,"estimate"], 
                                lowerCL = Fluctuation_Model_CVR_Estimates_Reorder[,"ci.lb"], 
                                upperCL = Fluctuation_Model_CVR_Estimates_Reorder[,"ci.ub"], 
                                K = fluctuation_k[,1], 
                                group_no = fluctuation_group_no[,1], 
                                row.names = fluctuation_rnames)
fluctuation_table$name <- row.names(fluctuation_table)

fluctuation_raw_mean <- c(unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                          select("InCVR"))), 
                          unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                          select("InCVR"))), 
                          unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                          select("InCVR"))), 
                          unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stochastic") %>% 
                                          select("InCVR"))))

fluctuation_raw_name <- c(replicate(541, "Sinusoidal (Sine Curve)"), 
                          replicate(607, "Alternating"), 
                          replicate(218, "Stepwise"), 
                          replicate(21, "Stochastic"))

fluctuation_raw_df <- data.frame("Model" = fluctuation_raw_name, 
                                 "Effect" = fluctuation_raw_mean)

# Graph code - Combined

Fluctuation_Order <- c("Stochastic", "Stepwise", 
                       "Alternating", "Sinusoidal (Sine Curve)")

density_fluctuation_CVR <- fluctuation_table %>% mutate(name = fct_relevel(name, Fluctuation_Order)) %>%
                           ggplot() +
                           geom_density_ridges(data = fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Fluctuation_Order)), 
                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                           geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                          size = 1) +
                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                           size = 1, fatten = 2) +
                           theme_bw() +
                           guides(fill = "none", colour = "none") +
                           labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                           theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                            vjust = c(-2.7, -2.7, -2.7, -0.8))) +
                           theme(axis.text.x = element_text(margin = margin(b = 5))) +
                           theme(axis.ticks = element_blank()) +
                           theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                           theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                           scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                           scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                           scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                           coord_cartesian(xlim = c(-1, 1)) +
                           annotate('text',  x = 1, y = (seq(1, dim(fluctuation_table)[1], 1)+0.4),
                           label= paste("italic(k)==", c(fluctuation_table["Stochastic", "K"], 
                                                         fluctuation_table["Stepwise", "K"], 
                                                         fluctuation_table["Alternating", "K"], 
                                                         fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                       c(fluctuation_table["Stochastic", "group_no"], 
                                                         fluctuation_table["Stepwise", "group_no"], 
                                                         fluctuation_table["Alternating", "group_no"], 
                                                         fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                           geom_label(aes(label=c(paste(format(round(mean(exp(Fluctuation_Model_CVR_Estimates["Stochastic", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                  paste(format(round(mean(exp(Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                  paste(format(round(mean(exp(Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                  paste(format(round(mean(exp(Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                      x = -0.75, y = (seq(1, dim(fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_fluctuation_CVR #(400x400)

# Preparing Graph - Part 1

fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)", "Alternating")

fluctuation_k_1 <- data.frame("k" = c(Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                      Fluctuation_Exploration["Alternating", "Freq"]), 
                              row.names = fluctuation_rnames_1)

fluctuation_group_no_1 <- data.frame("Spp No." = c(Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                   Fluctuation_Species_Count["Alternating", "Freq"]), 
                                     row.names = fluctuation_rnames_1)

fluctuation_study_1 <- data.frame("Study" = c(Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                              Fluctuation_Study_Count["Alternating", "Freq"]), 
                                  row.names = fluctuation_rnames_1)

Fluctuation_Model_CVR_Estimates_Reorder_1 <- Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

fluctuation_table_1 <- data.frame(estimate = Fluctuation_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                  lowerCL = Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                  upperCL = Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                  K = fluctuation_k_1[,1], 
                                  group_no = fluctuation_group_no_1[,1], 
                                  row.names = fluctuation_rnames_1)
fluctuation_table_1$name <- row.names(fluctuation_table_1)

fluctuation_raw_mean_1 <- c(unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                          select("InCVR"))), 
                            unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                          select("InCVR"))))

fluctuation_raw_name_1 <- c(replicate(541, "Sinusoidal (Sine Curve)"), 
                            replicate(607, "Alternating"))

fluctuation_raw_df_1 <- data.frame("Model" = fluctuation_raw_name_1, 
                                   "Effect" = fluctuation_raw_mean_1)

# Graph code - Part 1

Fluctuation_Order_1 <- c("Alternating", "Sinusoidal (Sine Curve)")

density_fluctuation_CVR_1 <- fluctuation_table_1 %>% mutate(name = fct_relevel(name, Fluctuation_Order_1)) %>%
                             ggplot() +
                             geom_density_ridges(data = fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Fluctuation_Order_1)), 
                                                 aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                             geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                            size = 1) +
                             geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                             size = 1, fatten = 2) +
                             theme_bw() +
                             guides(fill = "none", colour = "none") +
                             labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                             theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                   vjust = c(-2.7, -0.8))) +
                             theme(axis.text.x = element_text(margin = margin(b = 5))) +
                             theme(axis.ticks = element_blank()) +
                             theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                             scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                             scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                             coord_cartesian(xlim = c(-1, 1)) +
                             annotate('text',  x = 1, y = (seq(1, dim(fluctuation_table_1)[1], 1)+0.4),
                             label= paste("italic(k)==", c(fluctuation_table_1["Alternating", "K"], 
                                                           fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                         c(fluctuation_table_1["Alternating", "group_no"], 
                                                           fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                          ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                        x = -0.75, y = (seq(1, dim(fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_fluctuation_CVR_1 #(400x240)

# Preparing Graph - Part 2

fluctuation_rnames_2 <- c("Stepwise", "Stochastic")

fluctuation_k_2 <- data.frame("k" = c(Fluctuation_Exploration["Stepwise", "Freq"], 
                                      Fluctuation_Exploration["Stochastic", "Freq"]), 
                              row.names = fluctuation_rnames_2)

fluctuation_group_no_2 <- data.frame("Spp No." = c(Fluctuation_Species_Count["Stepwise", "Freq"], 
                                                   Fluctuation_Species_Count["Stochastic", "Freq"]), 
                                     row.names = fluctuation_rnames_2)

fluctuation_study_2 <- data.frame("Study" = c(Fluctuation_Study_Count["Stepwise", "Freq"], 
                                              Fluctuation_Study_Count["Stochastic", "Freq"]), 
                                  row.names = fluctuation_rnames_2)

Fluctuation_Model_CVR_Estimates_Reorder_2 <- Fluctuation_Model_CVR_Estimates[c("Stepwise", "Stochastic"), ]

fluctuation_table_2 <- data.frame(estimate = Fluctuation_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                  lowerCL = Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                  upperCL = Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                  K = fluctuation_k_2[,1], 
                                  group_no = fluctuation_group_no_2[,1], 
                                  row.names = fluctuation_rnames_2)
fluctuation_table_2$name <- row.names(fluctuation_table_2)

fluctuation_raw_mean_2 <- c(unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                          select("InCVR"))), 
                            unlist(unname(Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stochastic") %>% 
                                          select("InCVR"))))

fluctuation_raw_name_2 <- c(replicate(218, "Stepwise"), 
                            replicate(21, "Stochastic"))

fluctuation_raw_df_2 <- data.frame("Model" = fluctuation_raw_name_2, 
                                   "Effect" = fluctuation_raw_mean_2)

# Graph code - Part 2

Fluctuation_Order_2 <- c("Stochastic", "Stepwise")

density_fluctuation_CVR_2 <- fluctuation_table_2 %>% mutate(name = fct_relevel(name, Fluctuation_Order_2)) %>%
                             ggplot() +
                             geom_density_ridges(data = fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Fluctuation_Order_2)), 
                                                 aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                             geom_linerange(aes(y = rev(seq(1, dim(fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                            size = 1) +
                             geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                             size = 1, fatten = 2) +
                             theme_bw() +
                             guides(fill = "none", colour = "none") +
                             labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                             theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                   vjust = c(-2.7, -2.7))) +
                             theme(axis.text.x = element_text(margin = margin(b = 5))) +
                             theme(axis.ticks = element_blank()) +
                             theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                             scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                             scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                             coord_cartesian(xlim = c(-1, 1)) +
                             annotate('text',  x = 1, y = (seq(1, dim(fluctuation_table_2)[1], 1)+0.4),
                             label= paste("italic(k)==", c(fluctuation_table_2["Stochastic", "K"], 
                                                           fluctuation_table_2["Stepwise", "K"]), "~","(", 
                                                         c(fluctuation_table_2["Stochastic", "group_no"], 
                                                           fluctuation_table_2["Stepwise", "group_no"]), 
                                          ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(Fluctuation_Model_CVR_Estimates["Stochastic", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                        x = -0.75, y = (seq(1, dim(fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_fluctuation_CVR_2 #(400x240)

##### Overall Model - Trait Meta-Regression - CVR #####
Trait_Exploration <- data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Trait_Exploration) <- Trait_Exploration$Trait_Category

Trait_Species_Count <- data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>% 
                       filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Trait_Species_Count) <- Trait_Species_Count$Trait_Category

Trait_Study_Count <- data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>% 
                     filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Trait_Study_Count) <- Trait_Study_Count$Trait_Category

run <- FALSE
system.time( #  43ish minutes
  if(run){
    Trait_Model_CVR <- metafor::rma.mv(InCVR, V = VCV_InCVR, test = "t", dfs = "contain",
                                       mods = ~ Trait_Category - 1,
                                       random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                     ~1|Shared_Animal_Number, ~1|Measurement), 
                                       R = list(phylo=A_cor), data = data, method = "REML", sparse = TRUE, 
                                       control=list(rel.tol=1e-9))
    saveRDS(Trait_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Trait_Model_CVR.rds")
  } else {
            Trait_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Trait_Model_CVR.rds")})

Trait_Model_CVR_rob <- robust(Trait_Model_CVR, cluster = data$Study_ID, adjust = TRUE)

Trait_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Trait_Model_CVR$b), 15, 100),
                                        estimate = Trait_Model_CVR$b, ci.lb = Trait_Model_CVR$ci.lb, 
                                        ci.ub = Trait_Model_CVR$ci.ub)
rownames(Trait_Model_CVR_Estimates) <- Trait_Model_CVR_Estimates$Category
Trait_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Trait_Model_CVR), 2))

# Preparing Graph - Combined

trait_rnames <- c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits", 
                  "Morphology", "Physiological", "Population")

trait_k <- data.frame("k" = c(Trait_Exploration["Behavioural", "Freq"], 
                              Trait_Exploration["Biochemical Assay", "Freq"], 
                              Trait_Exploration["Gene Expression", "Freq"], 
                              Trait_Exploration["Life-History Traits", "Freq"], 
                              Trait_Exploration["Morphology", "Freq"], 
                              Trait_Exploration["Physiological", "Freq"], 
                              Trait_Exploration["Population", "Freq"]), 
                      row.names = trait_rnames)

trait_group_no <- data.frame("Spp No." = c(Trait_Species_Count["Behavioural", "Freq"], 
                                           Trait_Species_Count["Biochemical Assay", "Freq"], 
                                           Trait_Species_Count["Gene Expression", "Freq"], 
                                           Trait_Species_Count["Life-History Traits", "Freq"],
                                           Trait_Species_Count["Morphology", "Freq"],
                                           Trait_Species_Count["Physiological", "Freq"],
                                           Trait_Species_Count["Population", "Freq"]), 
                             row.names = trait_rnames)

trait_study <- data.frame("Study" = c(Trait_Study_Count["Behavioural", "Freq"], 
                                      Trait_Study_Count["Biochemical Assay", "Freq"], 
                                      Trait_Study_Count["Gene Expression", "Freq"], 
                                      Trait_Study_Count["Life-History Traits", "Freq"],
                                      Trait_Study_Count["Morphology", "Freq"],
                                      Trait_Study_Count["Physiological", "Freq"],
                                      Trait_Study_Count["Population", "Freq"]), 
                             row.names = trait_rnames)

trait_table <- data.frame(estimate = Trait_Model_CVR_Estimates[,"estimate"], 
                          lowerCL = Trait_Model_CVR_Estimates[,"ci.lb"], 
                          upperCL = Trait_Model_CVR_Estimates[,"ci.ub"], 
                          K = trait_k[,1], 
                          group_no = trait_group_no[,1], 
                          row.names = trait_rnames)
trait_table$name <- row.names(trait_table)

trait_raw_mean <- c(unlist(unname(data %>% filter(`Trait_Category` == "Behavioural") %>% 
                                    select("InCVR"))), 
                    unlist(unname(data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                    select("InCVR"))), 
                    unlist(unname(data %>% filter(`Trait_Category` == "Gene Expression") %>% 
                                    select("InCVR"))), 
                    unlist(unname(data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                    select("InCVR"))), 
                    unlist(unname(data %>% filter(`Trait_Category` == "Morphology") %>% 
                                    select("InCVR"))),
                    unlist(unname(data %>% filter(`Trait_Category` == "Physiological") %>% 
                                    select("InCVR"))),
                    unlist(unname(data %>% filter(`Trait_Category` == "Population") %>% 
                                    select("InCVR"))))

trait_raw_name <- c(replicate(38, "Behavioural"), 
                    replicate(154, "Biochemical Assay"), 
                    replicate(50, "Gene Expression"), 
                    replicate(508, "Life-history Traits"), 
                    replicate(376, "Morphology"),
                    replicate(198, "Physiological"),
                    replicate(168, "Population"))

trait_raw_df <- data.frame("Model" = trait_raw_name, 
                           "Effect" = trait_raw_mean)

# Graph code - Combined

Trait_Order <- c("Population", "Physiological", "Morphology", "Life-history Traits",  
                 "Gene Expression", "Biochemical Assay", "Behavioural")

density_trait_CVR <- trait_table %>% mutate(name = fct_relevel(name, Trait_Order)) %>%
                     ggplot() +
                     geom_density_ridges(data = trait_raw_df %>% mutate(Model = fct_relevel(Model, Trait_Order)), 
                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                         scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                     geom_linerange(aes(y = rev(seq(1, dim(trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                    size = 1) +
                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                     size = 1, fatten = 2) +
                     theme_bw() +
                     guides(fill = "none", colour = "none") +
                     labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                      vjust = c(-2.7, -2.7, -2.7, -0.8, -0.8, -0.8, -2.7))) +
                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                     theme(axis.ticks = element_blank()) +
                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                     scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                     scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                     coord_cartesian(xlim = c(-1, 1)) +
                     annotate('text',  x = 1, y = (seq(1, dim(trait_table)[1], 1)+0.4),
                     label= paste("italic(k)==", c(trait_table["Population", "K"], 
                                                   trait_table["Physiological", "K"], 
                                                   trait_table["Morphology", "K"], 
                                                   trait_table["Life-history Traits", "K"],
                                                   trait_table["Gene Expression", "K"],
                                                   trait_table["Biochemical Assay", "K"],
                                                   trait_table["Behavioural", "K"]), "~","(", 
                                                 c(trait_table["Population", "group_no"], 
                                                   trait_table["Physiological", "group_no"], 
                                                   trait_table["Morphology", "group_no"], 
                                                   trait_table["Life-history Traits", "group_no"],
                                                   trait_table["Gene Expression", "group_no"],
                                                   trait_table["Biochemical Assay", "group_no"],
                                                   trait_table["Behavioural", "group_no"]), 
                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                     geom_label(aes(label=c(paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Population", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Gene Expression", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Behavioural", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                x = -0.75, y = (seq(1, dim(trait_table)[1], 1)+0.4)), size = 3.5)

density_trait_CVR #(400x640)

# Preparing Graph - Part 1

trait_rnames_1 <- c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits")

trait_k_1 <- data.frame("k" = c(Trait_Exploration["Behavioural", "Freq"], 
                                Trait_Exploration["Biochemical Assay", "Freq"], 
                                Trait_Exploration["Gene Expression", "Freq"], 
                                Trait_Exploration["Life-History Traits", "Freq"]), 
                        row.names = trait_rnames_1)

trait_group_no_1 <- data.frame("Spp No." = c(Trait_Species_Count["Behavioural", "Freq"], 
                                             Trait_Species_Count["Biochemical Assay", "Freq"], 
                                             Trait_Species_Count["Gene Expression", "Freq"], 
                                             Trait_Species_Count["Life-History Traits", "Freq"]), 
                               row.names = trait_rnames_1)

trait_study_1 <- data.frame("Study" = c(Trait_Study_Count["Behavioural", "Freq"], 
                                        Trait_Study_Count["Biochemical Assay", "Freq"], 
                                        Trait_Study_Count["Gene Expression", "Freq"], 
                                        Trait_Study_Count["Life-History Traits", "Freq"]), 
                            row.names = trait_rnames_1)

Trait_Model_CVR_Estimates_Reorder_1 <- Trait_Model_CVR_Estimates[c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-History Traits"), ]

trait_table_1 <- data.frame(estimate = Trait_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                            lowerCL = Trait_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                            upperCL = Trait_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                            K = trait_k_1[,1], 
                            group_no = trait_group_no_1[,1], 
                            row.names = trait_rnames_1)
trait_table_1$name <- row.names(trait_table_1)

trait_raw_mean_1 <- c(unlist(unname(data %>% filter(`Trait_Category` == "Behavioural") %>% 
                                    select("InCVR"))), 
                      unlist(unname(data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                    select("InCVR"))), 
                      unlist(unname(data %>% filter(`Trait_Category` == "Gene Expression") %>% 
                                    select("InCVR"))), 
                      unlist(unname(data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                    select("InCVR"))))

trait_raw_name_1 <- c(replicate(38, "Behavioural"), 
                      replicate(154, "Biochemical Assay"), 
                      replicate(50, "Gene Expression"), 
                      replicate(508, "Life-history Traits"))

trait_raw_df_1 <- data.frame("Model" = trait_raw_name_1, 
                             "Effect" = trait_raw_mean_1)

# Graph code - Part 1

Trait_Order_1 <- c("Life-history Traits", "Gene Expression", "Biochemical Assay", "Behavioural")

density_trait_CVR_1 <- trait_table_1 %>% mutate(name = fct_relevel(name, Trait_Order_1)) %>%
                       ggplot() +
                       geom_density_ridges(data = trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Trait_Order_1)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                       size = 1, fatten = 2) +
                       theme_bw() +
                       guides(fill = "none", colour = "none") +
                       labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                             vjust = c(-0.8, -0.8, -0.8, -2.7))) +
                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                       theme(axis.ticks = element_blank()) +
                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                       scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                       scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                       coord_cartesian(xlim = c(-1, 1)) +
                       annotate('text',  x = 1, y = (seq(1, dim(trait_table_1)[1], 1)+0.4),
                       label= paste("italic(k)==", c(trait_table["Life-history Traits", "K"],
                                                     trait_table["Gene Expression", "K"],
                                                     trait_table["Biochemical Assay", "K"],
                                                     trait_table["Behavioural", "K"]), "~","(", 
                                                   c(trait_table["Life-history Traits", "group_no"],
                                                     trait_table["Gene Expression", "group_no"],
                                                     trait_table["Biochemical Assay", "group_no"],
                                                     trait_table["Behavioural", "group_no"]), 
                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                         geom_label(aes(label=c(paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Gene Expression", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Behavioural", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                        x = -0.75, y = (seq(1, dim(trait_table_1)[1], 1)+0.4)), size = 3.5)

density_trait_CVR_1 #(400x400)

# Preparing Graph - Part 2

trait_rnames_2 <- c("Morphology", "Physiological", "Population")

trait_k_2 <- data.frame("k" = c(Trait_Exploration["Morphology", "Freq"], 
                                Trait_Exploration["Physiological", "Freq"], 
                                Trait_Exploration["Population", "Freq"]), 
                        row.names = trait_rnames_2)

trait_group_no_2 <- data.frame("Spp No." = c(Trait_Species_Count["Morphology", "Freq"],
                                             Trait_Species_Count["Physiological", "Freq"],
                                             Trait_Species_Count["Population", "Freq"]), 
                               row.names = trait_rnames_2)

trait_study_2 <- data.frame("Study" = c(Trait_Study_Count["Morphology", "Freq"],
                                        Trait_Study_Count["Physiological", "Freq"],
                                        Trait_Study_Count["Population", "Freq"]), 
                            row.names = trait_rnames_2)

Trait_Model_CVR_Estimates_Reorder_2 <- Trait_Model_CVR_Estimates[c("Morphology", "Physiological", "Population"), ]

trait_table_2 <- data.frame(estimate = Trait_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                            lowerCL = Trait_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                            upperCL = Trait_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                            K = trait_k_2[,1], 
                            group_no = trait_group_no_2[,1], 
                            row.names = trait_rnames_2)
trait_table_2$name <- row.names(trait_table_2)

trait_raw_mean_2 <- c(unlist(unname(data %>% filter(`Trait_Category` == "Morphology") %>% 
                                    select("InCVR"))),
                      unlist(unname(data %>% filter(`Trait_Category` == "Physiological") %>% 
                                    select("InCVR"))),
                      unlist(unname(data %>% filter(`Trait_Category` == "Population") %>% 
                                    select("InCVR"))))

trait_raw_name_2 <- c(replicate(376, "Morphology"),
                      replicate(198, "Physiological"),
                      replicate(168, "Population"))

trait_raw_df_2 <- data.frame("Model" = trait_raw_name_2, 
                             "Effect" = trait_raw_mean_2)

# Graph code - Part 2

Trait_Order_2 <- c("Population", "Physiological", "Morphology")

density_trait_CVR_2 <- trait_table_2 %>% mutate(name = fct_relevel(name, Trait_Order_2)) %>%
                       ggplot() +
                       geom_density_ridges(data = trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Trait_Order_2)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                       size = 1, fatten = 2) +
                       theme_bw() +
                       guides(fill = "none", colour = "none") +
                       labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                             vjust = c(-2.7, -2.7, -2.7))) +
                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                       theme(axis.ticks = element_blank()) +
                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                       scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                       scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                       coord_cartesian(xlim = c(-1, 1)) +
                       annotate('text',  x = 1, y = (seq(1, dim(trait_table_2)[1], 1)+0.4),
                       label= paste("italic(k)==", c(trait_table_2["Population", "K"], 
                                                     trait_table_2["Physiological", "K"], 
                                                     trait_table_2["Morphology", "K"]), "~","(", 
                                                   c(trait_table_2["Population", "group_no"], 
                                                     trait_table_2["Physiological", "group_no"], 
                                                     trait_table_2["Morphology", "group_no"]), 
                                   ")"), parse = TRUE, hjust = "right", size = 3.5) +
                       geom_label(aes(label=c(paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Population", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                              paste(format(round(mean(exp(Trait_Model_CVR_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                  x = -0.75, y = (seq(1, dim(trait_table_2)[1], 1)+0.4)), size = 3.5)

density_trait_CVR_2 #(400x320)

##### Overall Model - Class Meta-Regression - CVR #####
Class_Exploration <- data %>% select("Class") %>% table() %>% data.frame()
rownames(Class_Exploration) <- Class_Exploration$Class

Class_Data <- data %>% filter(Class != "Clitellata" &
                              Class != "Collembola")

Class_Species_Count <- Class_Data %>% select("Scientific_Name", "Class") %>% table() %>% data.frame() %>% 
                       filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Class_Species_Count) <- Class_Species_Count$Class

Class_Study_Count <- Class_Data %>% select("Study_ID", "Class") %>% table() %>% data.frame() %>% 
                     filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Class_Study_Count) <- Class_Study_Count$Class

Class_Species <- Class_Data %>% select("phylo") %>% unique()

Class_A_cor <- as.data.frame(A_cor)
Class_A_cor <- Class_A_cor[c(Class_Species$phylo), c(Class_Species$phylo)]
Class_A_cor <- as.matrix(Class_A_cor)

Class_VCV_InCVR <- make_VCV_matrix(Class_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  23ish minutes
  if(run){
    Class_Model_CVR <- metafor::rma.mv(InCVR, V = Class_VCV_InCVR, test = "t", dfs = "contain",
                                       mods = ~ Class - 1,
                                       random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                     ~1|Shared_Animal_Number, ~1|Measurement), 
                                       R = list(phylo=Class_A_cor), data = Class_Data, method = "REML", sparse = TRUE, 
                                       control=list(rel.tol=1e-9))
    saveRDS(Class_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Class_Model_CVR.rds")
  } else {
            Class_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Class_Model_CVR.rds")})

Class_Model_CVR_rob <- robust(Class_Model_CVR, cluster = Class_Data$Study_ID, adjust = TRUE)

Class_Model_CVR_Estimates <- data.frame(Class = substr(row.names(Class_Model_CVR$b), 6, 100),
                                        estimate = Class_Model_CVR$b, ci.lb = Class_Model_CVR$ci.lb, 
                                        ci.ub = Class_Model_CVR$ci.ub)
rownames(Class_Model_CVR_Estimates) <- Class_Model_CVR_Estimates$Class
Class_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Class_Model_CVR), 2))

# Preparing Graph - Combined

class_rnames <- c("Actinopteri", "Amphibia", "Anthozoa", "Arachnida", 
                  "Bivalvia", "Branchiopoda", "Gastropoda", "Holothuroidea", 
                  "Insecta", "Malacostraca")

class_k <- data.frame("k" = c(Class_Exploration["Actinopteri", "Freq"], 
                              Class_Exploration["Amphibia", "Freq"], 
                              Class_Exploration["Anthozoa", "Freq"], 
                              Class_Exploration["Arachnida", "Freq"], 
                              Class_Exploration["Bivalvia", "Freq"], 
                              Class_Exploration["Branchiopoda", "Freq"], 
                              Class_Exploration["Gastropoda", "Freq"], 
                              Class_Exploration["Holothuroidea", "Freq"],
                              Class_Exploration["Insecta", "Freq"],
                              Class_Exploration["Malacostraca", "Freq"]), 
                      row.names = class_rnames)

class_group_no <- data.frame("Spp No." = c(Class_Species_Count["Actinopteri", "Freq"], 
                                           Class_Species_Count["Amphibia", "Freq"], 
                                           Class_Species_Count["Anthozoa", "Freq"], 
                                           Class_Species_Count["Arachnida", "Freq"],
                                           Class_Species_Count["Bivalvia", "Freq"],
                                           Class_Species_Count["Branchiopoda", "Freq"],
                                           Class_Species_Count["Gastropoda", "Freq"], 
                                           Class_Species_Count["Holothuroidea", "Freq"],
                                           Class_Species_Count["Insecta", "Freq"],
                                           Class_Species_Count["Malacostraca", "Freq"]), 
                             row.names = class_rnames)

class_study <- data.frame("Study" = c(Class_Study_Count["Actinopteri", "Freq"], 
                                      Class_Study_Count["Amphibia", "Freq"], 
                                      Class_Study_Count["Anthozoa", "Freq"], 
                                      Class_Study_Count["Arachnida", "Freq"],
                                      Class_Study_Count["Bivalvia", "Freq"],
                                      Class_Study_Count["Branchiopoda", "Freq"],
                                      Class_Study_Count["Gastropoda", "Freq"], 
                                      Class_Study_Count["Holothuroidea", "Freq"],
                                      Class_Study_Count["Insecta", "Freq"],
                                      Class_Study_Count["Malacostraca", "Freq"]), 
                             row.names = class_rnames)

class_table <- data.frame(estimate = Class_Model_CVR_Estimates[,"estimate"], 
                          lowerCL = Class_Model_CVR_Estimates[,"ci.lb"], 
                          upperCL = Class_Model_CVR_Estimates[,"ci.ub"], 
                          K = class_k[,1], 
                          group_no = class_group_no[,1], 
                          row.names = class_rnames)
class_table$name <- row.names(class_table)

class_raw_mean <- c(unlist(unname(Class_Data %>% filter(`Class` == "Actinopteri") %>% 
                                    select("InCVR"))), 
                    unlist(unname(Class_Data %>% filter(`Class` == "Amphibia") %>% 
                                    select("InCVR"))), 
                    unlist(unname(Class_Data %>% filter(`Class` == "Anthozoa") %>% 
                                    select("InCVR"))), 
                    unlist(unname(Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                    select("InCVR"))), 
                    unlist(unname(Class_Data %>% filter(`Class` == "Bivalvia") %>% 
                                    select("InCVR"))),
                    unlist(unname(Class_Data %>% filter(`Class` == "Branchiopoda") %>% 
                                    select("InCVR"))),
                    unlist(unname(Class_Data %>% filter(`Class` == "Gastropoda") %>% 
                                    select("InCVR"))), 
                    unlist(unname(Class_Data %>% filter(`Class` == "Holothuroidea") %>% 
                                    select("InCVR"))),
                    unlist(unname(Class_Data %>% filter(`Class` == "Insecta") %>% 
                                    select("InCVR"))),
                    unlist(unname(Class_Data %>% filter(`Class` == "Malacostraca") %>% 
                                    select("InCVR"))))

class_raw_name <- c(replicate(150, "Actinopteri"), 
                    replicate(64, "Amphibia"), 
                    replicate(12, "Anthozoa"), 
                    replicate(126, "Arachnida"), 
                    replicate(14, "Bivalvia"),
                    replicate(21, "Branchiopoda"),
                    replicate(21, "Gastropoda"), 
                    replicate(42, "Holothuroidea"),
                    replicate(722, "Insecta"),
                    replicate(60, "Malacostraca"))

class_raw_df <- data.frame("Model" = class_raw_name, 
                           "Effect" = class_raw_mean)

# Graph code - Combined

Class_Order <- c("Malacostraca", "Insecta", "Holothuroidea", "Gastropoda",  
                 "Branchiopoda", "Bivalvia", "Arachnida", "Anthozoa", "Amphibia", "Actinopteri")

density_class_CVR <- class_table %>% mutate(name = fct_relevel(name, Class_Order)) %>%
                     ggplot() +
                     geom_density_ridges(data = class_raw_df %>% mutate(Model = fct_relevel(Model, Class_Order)), 
                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                         scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                     geom_linerange(aes(y = rev(seq(1, dim(class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                    size = 1) +
                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                     size = 1, fatten = 2) +
                     theme_bw() +
                     guides(fill = "none", colour = "none") +
                     labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                      vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7, 
                                                                -2.7, -2.7, -2.7, -2.7, -2.7))) +
                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                     theme(axis.ticks = element_blank()) +
                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                     scale_colour_manual(values = c("#6C87AB", "#6582A9", "#607C9F", "#4A6E9C", "#446692",
                                                    "#3D5E89", "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                     scale_fill_manual(values = c("#6C87AB", "#6582A9", "#607C9F", "#4A6E9C", "#446692",
                                                  "#3D5E89", "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                     coord_cartesian(xlim = c(-1, 1)) +
                     annotate('text',  x = 1, y = (seq(1, dim(class_table)[1], 1)+0.4),
                     label= paste("italic(k)==", c(class_table["Malacostraca", "K"], 
                                                   class_table["Insecta", "K"], 
                                                   class_table["Holothuroidea", "K"], 
                                                   class_table["Gastropoda", "K"],
                                                   class_table["Branchiopoda", "K"],
                                                   class_table["Bivalvia", "K"],
                                                   class_table["Arachnida", "K"], 
                                                   class_table["Anthozoa", "K"],
                                                   class_table["Amphibia", "K"],
                                                   class_table["Actinopteri", "K"]), "~","(", 
                                                 c(class_table["Malacostraca", "group_no"], 
                                                   class_table["Insecta", "group_no"], 
                                                   class_table["Holothuroidea", "group_no"], 
                                                   class_table["Gastropoda", "group_no"],
                                                   class_table["Branchiopoda", "group_no"],
                                                   class_table["Bivalvia", "group_no"],
                                                   class_table["Arachnida", "group_no"], 
                                                   class_table["Anthozoa", "group_no"],
                                                   class_table["Amphibia", "group_no"],
                                                   class_table["Actinopteri", "group_no"]), 
                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                     geom_label(aes(label=c(paste(format(round(mean(exp(Class_Model_CVR_Estimates["Malacostraca", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(Class_Model_CVR_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(Class_Model_CVR_Estimates["Holothuroidea", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(Class_Model_CVR_Estimates["Gastropoda", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(Class_Model_CVR_Estimates["Branchiopoda", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(Class_Model_CVR_Estimates["Bivalvia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(Class_Model_CVR_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                            paste(format(round(mean(exp(Class_Model_CVR_Estimates["Anthozoa", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(Class_Model_CVR_Estimates["Amphibia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                            paste(format(round(mean(exp(Class_Model_CVR_Estimates["Actinopteri", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                x = -0.75, y = (seq(1, dim(class_table)[1], 1)+0.4)), size = 3.5)

density_class_CVR #(400x880)

# Preparing Graph - Part 1

class_rnames_1 <- c("Actinopteri", "Amphibia", "Anthozoa", "Arachnida", "Bivalvia")

class_k_1 <- data.frame("k" = c(Class_Exploration["Actinopteri", "Freq"], 
                                Class_Exploration["Amphibia", "Freq"], 
                                Class_Exploration["Anthozoa", "Freq"], 
                                Class_Exploration["Arachnida", "Freq"], 
                                Class_Exploration["Bivalvia", "Freq"]), 
                        row.names = class_rnames_1)

class_group_no_1 <- data.frame("Spp No." = c(Class_Species_Count["Actinopteri", "Freq"], 
                                             Class_Species_Count["Amphibia", "Freq"], 
                                             Class_Species_Count["Anthozoa", "Freq"], 
                                             Class_Species_Count["Arachnida", "Freq"],
                                             Class_Species_Count["Bivalvia", "Freq"]), 
                               row.names = class_rnames_1)

class_study_1 <- data.frame("Study" = c(Class_Study_Count["Actinopteri", "Freq"], 
                                        Class_Study_Count["Amphibia", "Freq"], 
                                        Class_Study_Count["Anthozoa", "Freq"], 
                                        Class_Study_Count["Arachnida", "Freq"],
                                        Class_Study_Count["Bivalvia", "Freq"]), 
                            row.names = class_rnames_1)

Class_Model_CVR_Estimates_Reorder_1 <- Class_Model_CVR_Estimates[c("Actinopteri", "Amphibia", "Anthozoa", "Arachnida", "Bivalvia"), ]

class_table_1 <- data.frame(estimate = Class_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                            lowerCL = Class_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                            upperCL = Class_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                            K = class_k_1[,1], 
                            group_no = class_group_no_1[,1], 
                            row.names = class_rnames_1)
class_table_1$name <- row.names(class_table_1)

class_raw_mean_1 <- c(unlist(unname(Class_Data %>% filter(`Class` == "Actinopteri") %>% 
                                    select("InCVR"))), 
                      unlist(unname(Class_Data %>% filter(`Class` == "Amphibia") %>% 
                                    select("InCVR"))), 
                      unlist(unname(Class_Data %>% filter(`Class` == "Anthozoa") %>% 
                                    select("InCVR"))), 
                      unlist(unname(Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                    select("InCVR"))), 
                      unlist(unname(Class_Data %>% filter(`Class` == "Bivalvia") %>% 
                                    select("InCVR"))))

class_raw_name_1 <- c(replicate(150, "Actinopteri"), 
                      replicate(64, "Amphibia"), 
                      replicate(12, "Anthozoa"), 
                      replicate(126, "Arachnida"), 
                      replicate(14, "Bivalvia"))

class_raw_df_1 <- data.frame("Model" = class_raw_name_1, 
                             "Effect" = class_raw_mean_1)

# Graph code - Part 1

Class_Order_1 <- c("Bivalvia", "Arachnida", "Anthozoa", "Amphibia", "Actinopteri")

density_class_CVR_1 <- class_table_1 %>% mutate(name = fct_relevel(name, Class_Order_1)) %>%
                       ggplot() +
                       geom_density_ridges(data = class_raw_df_1 %>% mutate(Model = fct_relevel(Model, Class_Order_1)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                       size = 1, fatten = 2) +
                       theme_bw() +
                       guides(fill = "none", colour = "none") +
                       labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                             vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7))) +
                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                       theme(axis.ticks = element_blank()) +
                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                       scale_colour_manual(values = c("#3D5E89", "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                       scale_fill_manual(values = c("#3D5E89", "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                       coord_cartesian(xlim = c(-1, 1)) +
                       annotate('text',  x = 1, y = (seq(1, dim(class_table_1)[1], 1)+0.4),
                       label= paste("italic(k)==", c(class_table_1["Bivalvia", "K"],
                                                     class_table_1["Arachnida", "K"], 
                                                     class_table_1["Anthozoa", "K"],
                                                     class_table_1["Amphibia", "K"],
                                                     class_table_1["Actinopteri", "K"]), "~","(", 
                                                   c(class_table_1["Bivalvia", "group_no"],
                                                     class_table_1["Arachnida", "group_no"], 
                                                     class_table_1["Anthozoa", "group_no"],
                                                     class_table_1["Amphibia", "group_no"],
                                                     class_table_1["Actinopteri", "group_no"]), 
                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                       geom_label(aes(label=c(paste(format(round(mean(exp(Class_Model_CVR_Estimates["Bivalvia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                              paste(format(round(mean(exp(Class_Model_CVR_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(Class_Model_CVR_Estimates["Anthozoa", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                              paste(format(round(mean(exp(Class_Model_CVR_Estimates["Amphibia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                              paste(format(round(mean(exp(Class_Model_CVR_Estimates["Actinopteri", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                  x = -0.75, y = (seq(1, dim(class_table_1)[1], 1)+0.4)), size = 3.5)

density_class_CVR_1 #(400x480)

# Preparing Graph - Part 2

class_rnames_2 <- c("Branchiopoda", "Gastropoda", "Holothuroidea", "Insecta", "Malacostraca")

class_k_2 <- data.frame("k" = c(Class_Exploration["Branchiopoda", "Freq"], 
                                Class_Exploration["Gastropoda", "Freq"], 
                                Class_Exploration["Holothuroidea", "Freq"],
                                Class_Exploration["Insecta", "Freq"],
                                Class_Exploration["Malacostraca", "Freq"]), 
                        row.names = class_rnames_2)

class_group_no_2 <- data.frame("Spp No." = c(Class_Species_Count["Branchiopoda", "Freq"],
                                             Class_Species_Count["Gastropoda", "Freq"], 
                                             Class_Species_Count["Holothuroidea", "Freq"],
                                             Class_Species_Count["Insecta", "Freq"],
                                             Class_Species_Count["Malacostraca", "Freq"]), 
                               row.names = class_rnames_2)

class_study_2 <- data.frame("Study" = c(Class_Study_Count["Branchiopoda", "Freq"],
                                        Class_Study_Count["Gastropoda", "Freq"], 
                                        Class_Study_Count["Holothuroidea", "Freq"],
                                        Class_Study_Count["Insecta", "Freq"],
                                        Class_Study_Count["Malacostraca", "Freq"]), 
                            row.names = class_rnames_2)

Class_Model_CVR_Estimates_Reorder_2 <- Class_Model_CVR_Estimates[c("Branchiopoda", "Gastropoda", "Holothuroidea", "Insecta", "Malacostraca"), ]

class_table_2 <- data.frame(estimate = Class_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                            lowerCL = Class_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                            upperCL = Class_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                            K = class_k_2[,1], 
                            group_no = class_group_no_2[,1], 
                            row.names = class_rnames_2)
class_table_2$name <- row.names(class_table_2)

class_raw_mean_2 <- c(unlist(unname(Class_Data %>% filter(`Class` == "Branchiopoda") %>% 
                                    select("InCVR"))),
                      unlist(unname(Class_Data %>% filter(`Class` == "Gastropoda") %>% 
                                    select("InCVR"))), 
                      unlist(unname(Class_Data %>% filter(`Class` == "Holothuroidea") %>% 
                                    select("InCVR"))),
                      unlist(unname(Class_Data %>% filter(`Class` == "Insecta") %>% 
                                    select("InCVR"))),
                      unlist(unname(Class_Data %>% filter(`Class` == "Malacostraca") %>% 
                                    select("InCVR"))))

class_raw_name_2 <- c(replicate(21, "Branchiopoda"),
                      replicate(21, "Gastropoda"), 
                      replicate(42, "Holothuroidea"),
                      replicate(722, "Insecta"),
                      replicate(60, "Malacostraca"))

class_raw_df_2 <- data.frame("Model" = class_raw_name_2, 
                             "Effect" = class_raw_mean_2)

# Graph code - Part 2

Class_Order_2 <- c("Malacostraca", "Insecta", "Holothuroidea", "Gastropoda", "Branchiopoda")

density_class_CVR_2 <- class_table_2 %>% mutate(name = fct_relevel(name, Class_Order_2)) %>%
                       ggplot() +
                       geom_density_ridges(data = class_raw_df_2 %>% mutate(Model = fct_relevel(Model, Class_Order_2)), 
                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                       geom_linerange(aes(y = rev(seq(1, dim(class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                      size = 1) +
                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                       size = 1, fatten = 2) +
                       theme_bw() +
                       guides(fill = "none", colour = "none") +
                       labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                             vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7))) +
                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                       theme(axis.ticks = element_blank()) +
                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                       scale_colour_manual(values = c("#6C87AB", "#6582A9", "#607C9F", "#4A6E9C", "#446692")) +
                       scale_fill_manual(values = c("#6C87AB", "#6582A9", "#607C9F", "#4A6E9C", "#446692")) +
                       coord_cartesian(xlim = c(-1, 1)) +
                       annotate('text',  x = 1, y = (seq(1, dim(class_table_2)[1], 1)+0.4),
                       label= paste("italic(k)==", c(class_table_2["Malacostraca", "K"], 
                                                     class_table_2["Insecta", "K"], 
                                                     class_table_2["Holothuroidea", "K"], 
                                                     class_table_2["Gastropoda", "K"],
                                                     class_table_2["Branchiopoda", "K"]), "~","(", 
                                                   c(class_table_2["Malacostraca", "group_no"], 
                                                     class_table_2["Insecta", "group_no"], 
                                                     class_table_2["Holothuroidea", "group_no"], 
                                                     class_table_2["Gastropoda", "group_no"],
                                                     class_table_2["Branchiopoda", "group_no"]), 
                                   ")"), parse = TRUE, hjust = "right", size = 3.5) +
                       geom_label(aes(label=c(paste(format(round(mean(exp(Class_Model_CVR_Estimates["Malacostraca", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(Class_Model_CVR_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                              paste(format(round(mean(exp(Class_Model_CVR_Estimates["Holothuroidea", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(Class_Model_CVR_Estimates["Gastropoda", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                              paste(format(round(mean(exp(Class_Model_CVR_Estimates["Branchiopoda", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                  x = -0.75, y = (seq(1, dim(class_table_2)[1], 1)+0.4)), size = 3.5)

density_class_CVR_2 #(400x480)

##### Overall Model - Specific Trait Meta-Regression - CVR #####
Specific_Trait_Exploration <- data %>% select("Measurement") %>% table() %>% data.frame()
Specific_Trait_Exploration <- Specific_Trait_Exploration %>% filter(Freq > 10)
rownames(Specific_Trait_Exploration) <- Specific_Trait_Exploration$Measurement

Specific_Trait_Data <- data %>% filter(Measurement == "Apparent Digestability Coefficient"| 
                                       Measurement == "Catalase Activity"|
                                       Measurement == "Cortisol"|
                                       Measurement == "Development Time"|
                                       Measurement == "Fecundity"|
                                       Measurement == "Food Consumption"|
                                       Measurement == "Gender"|
                                       Measurement == "Head Width"|
                                       Measurement == "hsp70"|
                                       Measurement == "Immune Defense"|
                                       Measurement == "Length"|
                                       Measurement == "Locomotor Performance"|
                                       Measurement == "Longevity"|
                                       Measurement == "Mass"|
                                       Measurement == "Metabolic Rate"|
                                       Measurement == "Mortality"|
                                       Measurement == "PO Activity"|
                                       Measurement == "Reproductive Rate"|
                                       Measurement == "SOD Activity"|
                                       Measurement == "Survival"|
                                       Measurement == "Tail Length"|
                                       Measurement == "Triglyceride")

Specific_Trait_Species_Count <- Specific_Trait_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>% 
                                filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Specific_Trait_Species_Count) <- Specific_Trait_Species_Count$Measurement

Specific_Trait_Study_Count <- Specific_Trait_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>% 
                              filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Specific_Trait_Study_Count) <- Specific_Trait_Study_Count$Measurement

Specific_Trait_Species <- Specific_Trait_Data %>% select("phylo") %>% unique()

Specific_Trait_A_cor <- as.data.frame(A_cor)
Specific_Trait_A_cor <- Specific_Trait_A_cor[c(Specific_Trait_Species$phylo), c(Specific_Trait_Species$phylo)]
Specific_Trait_A_cor <- as.matrix(Specific_Trait_A_cor)

Specific_Trait_VCV_InCVR <- make_VCV_matrix(Specific_Trait_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  17ish minutes
  if(run){
    Specific_Trait_Model_CVR <- metafor::rma.mv(InCVR, V = Specific_Trait_VCV_InCVR, test = "t", dfs = "contain",
                                                mods = ~ Measurement - 1,
                                                random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                              ~1|Shared_Animal_Number), 
                                                R = list(phylo=Specific_Trait_A_cor), data = Specific_Trait_Data, method = "REML", sparse = TRUE, 
                                                control=list(rel.tol=1e-9))
    saveRDS(Specific_Trait_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Specific_Trait_Model_CVR.rds")
  } else {
            Specific_Trait_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Specific_Trait_Model_CVR.rds")})

Specific_Trait_Model_CVR_rob <- robust(Specific_Trait_Model_CVR, cluster = Specific_Trait_Data$Study_ID, adjust = TRUE)

Specific_Trait_Model_CVR_Estimates <- data.frame(Trait = substr(row.names(Specific_Trait_Model_CVR$b), 12, 100),
                                                 estimate = Specific_Trait_Model_CVR$b, 
                                                 ci.lb = Specific_Trait_Model_CVR$ci.lb, 
                                                 ci.ub = Specific_Trait_Model_CVR$ci.ub)
rownames(Specific_Trait_Model_CVR_Estimates) <- Specific_Trait_Model_CVR_Estimates$Trait
Specific_Trait_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Specific_Trait_Model_CVR), 2))

# Preparing Graph - Part 1

specific_trait_rnames <- c("Apparent Digestibility Coefficient", "Catalase Activity", "Cortisol Levels", 
                           "Development Time", "Fecundity", "Food Consumption", "Head Width", 
                           "hsp70", "Immune Defense", "Length", "Locomotor Performance")

specific_trait_k <- data.frame("k" = c(Specific_Trait_Exploration["Apparent Digestability Coefficient", "Freq"], 
                                       Specific_Trait_Exploration["Catalase Activity", "Freq"], 
                                       Specific_Trait_Exploration["Cortisol", "Freq"], 
                                       Specific_Trait_Exploration["Development Time", "Freq"], 
                                       Specific_Trait_Exploration["Fecundity", "Freq"], 
                                       Specific_Trait_Exploration["Food Consumption", "Freq"], 
                                       Specific_Trait_Exploration["Head Width", "Freq"],
                                       Specific_Trait_Exploration["hsp70", "Freq"],
                                       Specific_Trait_Exploration["Immune Defense", "Freq"], 
                                       Specific_Trait_Exploration["Length", "Freq"], 
                                       Specific_Trait_Exploration["Locomotor Performance", "Freq"]), 
                               row.names = specific_trait_rnames)

specific_trait_group_no <- data.frame("Spp No." = c(Specific_Trait_Species_Count["Apparent Digestability Coefficient", "Freq"], 
                                                    Specific_Trait_Species_Count["Catalase Activity", "Freq"], 
                                                    Specific_Trait_Species_Count["Cortisol", "Freq"], 
                                                    Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                    Specific_Trait_Species_Count["Fecundity", "Freq"], 
                                                    Specific_Trait_Species_Count["Food Consumption", "Freq"], 
                                                    Specific_Trait_Species_Count["Head Width", "Freq"],
                                                    Specific_Trait_Species_Count["hsp70", "Freq"],
                                                    Specific_Trait_Species_Count["Immune Defense", "Freq"], 
                                                    Specific_Trait_Species_Count["Length", "Freq"], 
                                                    Specific_Trait_Species_Count["Locomotor Performance", "Freq"]), 
                                      row.names = specific_trait_rnames)

specific_trait_study <- data.frame("Study" = c(Specific_Trait_Study_Count["Apparent Digestability Coefficient", "Freq"], 
                                               Specific_Trait_Study_Count["Catalase Activity", "Freq"], 
                                               Specific_Trait_Study_Count["Cortisol", "Freq"], 
                                               Specific_Trait_Study_Count["Development Time", "Freq"], 
                                               Specific_Trait_Study_Count["Fecundity", "Freq"], 
                                               Specific_Trait_Study_Count["Food Consumption", "Freq"], 
                                               Specific_Trait_Study_Count["Head Width", "Freq"],
                                               Specific_Trait_Study_Count["hsp70", "Freq"],
                                               Specific_Trait_Study_Count["Immune Defense", "Freq"], 
                                               Specific_Trait_Study_Count["Length", "Freq"], 
                                               Specific_Trait_Study_Count["Locomotor Performance", "Freq"]), 
                                      row.names = specific_trait_rnames)

Specific_Trait_Model_CVR_EStimates_Reorder <- Specific_Trait_Model_CVR_Estimates[c("Apparent Digestability Coefficient", "Catalase Activity", "Cortisol", 
                                                                                   "Development Time", "Fecundity", "Food Consumption", "Head Width", 
                                                                                   "hsp70", "Immune Defense", "Length", "Locomotor Performance"), ]

specific_trait_table <- data.frame(estimate = Specific_Trait_Model_CVR_EStimates_Reorder[,"estimate"], 
                                   lowerCL = Specific_Trait_Model_CVR_EStimates_Reorder[,"ci.lb"], 
                                   upperCL = Specific_Trait_Model_CVR_EStimates_Reorder[,"ci.ub"], 
                                   K = specific_trait_k[,1], 
                                   group_no = specific_trait_group_no[,1], 
                                   row.names = specific_trait_rnames)
specific_trait_table$name <- row.names(specific_trait_table)

specific_trait_raw_mean <- c(unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Apparent Digestability Coefficient") %>% 
                                             select("InCVR"))), 
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Catalase Activity") %>% 
                                             select("InCVR"))), 
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Cortisol") %>% 
                                             select("InCVR"))), 
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                             select("InCVR"))), 
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Fecundity") %>% 
                                             select("InCVR"))),
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Food Consumption") %>% 
                                             select("InCVR"))),
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Head Width") %>% 
                                             select("InCVR"))),
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "hsp70") %>% 
                                             select("InCVR"))),
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Immune Defense") %>% 
                                             select("InCVR"))), 
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                             select("InCVR"))), 
                             unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Locomotor Performance") %>% 
                                             select("InCVR"))))

specific_trait_raw_name <- c(replicate(16, "Apparent Digestibility Coefficient"), 
                             replicate(11, "Catalase Activity"), 
                             replicate(13, "Cortisol Levels"), 
                             replicate(294, "Development Time"), 
                             replicate(72, "Fecundity"),
                             replicate(24, "Food Consumption"),
                             replicate(14, "Head Width"),
                             replicate(12, "hsp70"),
                             replicate(14, "Immune Defense"), 
                             replicate(89, "Length"), 
                             replicate(35, "Locomotor Performance"))

specific_trait_raw_df <- data.frame("Model" = specific_trait_raw_name, 
                                    "Effect" = specific_trait_raw_mean)

# Graph code - Part 1

Specific_Trait_Order <- c("Locomotor Performance", "Length", "Immune Defense", 
                          "hsp70", "Head Width", "Food Consumption", "Fecundity", 
                          "Development Time", "Cortisol Levels", "Catalase Activity", 
                          "Apparent Digestibility Coefficient")

density_specific_trait_CVR <- specific_trait_table %>% mutate(name = fct_relevel(name, Specific_Trait_Order)) %>%
                              ggplot() +
                              geom_density_ridges(data = specific_trait_raw_df %>% mutate(Model = fct_relevel(Model, Specific_Trait_Order)), 
                                                  aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                  scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                              geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                              geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                              size = 1, fatten = 2) +
                              theme_bw() +
                              guides(fill = "none", colour = "none") +
                              labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                              theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                               vjust = c(-0.8, -2.7, -0.8, -2.7, -2.7, -0.8, 
                                                                         -2.7, -0.8, -0.8, -0.8, -0.4))) +
                              theme(axis.text.x = element_text(margin = margin(b = 5))) +
                              theme(axis.ticks = element_blank()) +
                              theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                              theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                              scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                              scale_colour_manual(values = c("#3D5E89", "#355784", "#2B4E7A", "#234572", "#234373", "#1B3B6B", 
                                                             "#1C375F", "#163058", "#0D2A51", "#102C50", "#0F2643")) +
                              scale_fill_manual(values = c("#3D5E89", "#355784", "#2B4E7A", "#234572", "#234373", "#1B3B6B", 
                                                           "#1C375F", "#163058", "#0D2A51", "#102C50", "#0F2643")) +
                              coord_cartesian(xlim = c(-1, 1)) +
                              annotate('text',  x = 1, y = (seq(1, dim(specific_trait_table)[1], 1)+0.4),
                              label= paste("italic(k)==", c(specific_trait_table["Locomotor Performance", "K"], 
                                                            specific_trait_table["Length", "K"], 
                                                            specific_trait_table["Immune Defense", "K"], 
                                                            specific_trait_table["hsp70", "K"],
                                                            specific_trait_table["Head Width", "K"],
                                                            specific_trait_table["Food Consumption", "K"], 
                                                            specific_trait_table["Fecundity", "K"],
                                                            specific_trait_table["Development Time", "K"],
                                                            specific_trait_table["Cortisol Levels", "K"], 
                                                            specific_trait_table["Catalase Activity", "K"],
                                                            specific_trait_table["Apparent Digestibility Coefficient", "K"]), "~","(", 
                                                          c(specific_trait_table["Locomotor Performance", "group_no"], 
                                                            specific_trait_table["Length", "group_no"], 
                                                            specific_trait_table["Immune Defense", "group_no"], 
                                                            specific_trait_table["hsp70", "group_no"],
                                                            specific_trait_table["Head Width", "group_no"],
                                                            specific_trait_table["Food Consumption", "group_no"], 
                                                            specific_trait_table["Fecundity", "group_no"],
                                                            specific_trait_table["Development Time", "group_no"],
                                                            specific_trait_table["Cortisol Levels", "group_no"], 
                                                            specific_trait_table["Catalase Activity", "group_no"],
                                                            specific_trait_table["Apparent Digestibility Coefficient", "group_no"]), 
                                           ")"), parse = TRUE, hjust = "right", size = 3.5) +
                              geom_label(aes(label=c(paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Locomotor Performance", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                     paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                     paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Immune Defense", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                     paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["hsp70", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                     paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Head Width", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                     paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Food Consumption", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                     paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Fecundity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                     paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                     paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Cortisol", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                     paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Catalase Activity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                     paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Apparent Digestability Coefficient", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                         x = -0.75, y = (seq(1, dim(specific_trait_table)[1], 1)+0.4)), size = 3.5)

density_specific_trait_CVR #(400x960)

# Preparing Graph - Part 2

specific_trait_rnames_2 <- c("Longevity", "Mass", "Metabolic Rate", "Mortality", "PO Activity", "Reproductive Rate", 
                             "Sex", "SOD Activity", "Survival", "Tail Length", "Triglyceride")

specific_trait_k_2 <- data.frame("k" = c(Specific_Trait_Exploration["Longevity", "Freq"], 
                                         Specific_Trait_Exploration["Mass", "Freq"], 
                                         Specific_Trait_Exploration["Metabolic Rate", "Freq"], 
                                         Specific_Trait_Exploration["Mortality", "Freq"], 
                                         Specific_Trait_Exploration["PO Activity", "Freq"], 
                                         Specific_Trait_Exploration["Reproductive Rate", "Freq"],
                                         Specific_Trait_Exploration["Gender", "Freq"], 
                                         Specific_Trait_Exploration["SOD Activity", "Freq"],
                                         Specific_Trait_Exploration["Survival", "Freq"], 
                                         Specific_Trait_Exploration["Tail Length", "Freq"],
                                         Specific_Trait_Exploration["Triglyceride", "Freq"]), 
                                 row.names = specific_trait_rnames_2)

specific_trait_group_no_2 <- data.frame("Spp No." = c(Specific_Trait_Species_Count["Longevity", "Freq"], 
                                                      Specific_Trait_Species_Count["Mass", "Freq"], 
                                                      Specific_Trait_Species_Count["Metabolic Rate", "Freq"], 
                                                      Specific_Trait_Species_Count["Mortality", "Freq"], 
                                                      Specific_Trait_Species_Count["PO Activity", "Freq"], 
                                                      Specific_Trait_Species_Count["Reproductive Rate", "Freq"],
                                                      Specific_Trait_Species_Count["Gender", "Freq"], 
                                                      Specific_Trait_Species_Count["SOD Activity", "Freq"],
                                                      Specific_Trait_Species_Count["Survival", "Freq"], 
                                                      Specific_Trait_Species_Count["Tail Length", "Freq"],
                                                      Specific_Trait_Species_Count["Triglyceride", "Freq"]), 
                                        row.names = specific_trait_rnames_2)

specific_trait_study_2 <- data.frame("Study" = c(Specific_Trait_Study_Count["Longevity", "Freq"], 
                                                 Specific_Trait_Study_Count["Mass", "Freq"], 
                                                 Specific_Trait_Study_Count["Metabolic Rate", "Freq"], 
                                                 Specific_Trait_Study_Count["Mortality", "Freq"], 
                                                 Specific_Trait_Study_Count["PO Activity", "Freq"], 
                                                 Specific_Trait_Study_Count["Reproductive Rate", "Freq"],
                                                 Specific_Trait_Study_Count["Gender", "Freq"], 
                                                 Specific_Trait_Study_Count["SOD Activity", "Freq"],
                                                 Specific_Trait_Study_Count["Survival", "Freq"], 
                                                 Specific_Trait_Study_Count["Tail Length", "Freq"],
                                                 Specific_Trait_Study_Count["Triglyceride", "Freq"]), 
                                        row.names = specific_trait_rnames_2)

Specific_Trait_Model_CVR_EStimates_Reorder_2 <- Specific_Trait_Model_CVR_Estimates[c("Longevity", "Mass", "Metabolic Rate", "Mortality", "PO Activity", "Reproductive Rate", 
                                                                                     "Gender", "SOD Activity", "Survival", "Tail Length", "Triglyceride"), ]

specific_trait_table_2 <- data.frame(estimate = Specific_Trait_Model_CVR_EStimates_Reorder_2[,"estimate"], 
                                     lowerCL = Specific_Trait_Model_CVR_EStimates_Reorder_2[,"ci.lb"], 
                                     upperCL = Specific_Trait_Model_CVR_EStimates_Reorder_2[,"ci.ub"], 
                                     K = specific_trait_k_2[,1], 
                                     group_no = specific_trait_group_no_2[,1], 
                                     row.names = specific_trait_rnames_2)
specific_trait_table_2$name <- row.names(specific_trait_table_2)

specific_trait_raw_mean_2 <- c(unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Longevity") %>% 
                                             select("InCVR"))), 
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                             select("InCVR"))), 
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Metabolic Rate") %>% 
                                             select("InCVR"))),
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Mortality") %>% 
                                             select("InCVR"))),
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "PO Activity") %>% 
                                             select("InCVR"))), 
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Reproductive Rate") %>% 
                                             select("InCVR"))),
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Gender") %>% 
                                             select("InCVR"))), 
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "SOD Activity") %>% 
                                             select("InCVR"))),
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Survival") %>% 
                                             select("InCVR"))), 
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Tail Length") %>% 
                                             select("InCVR"))),
                               unlist(unname(Specific_Trait_Data %>% filter(`Measurement` == "Triglyceride") %>% 
                                             select("InCVR"))))

specific_trait_raw_name_2 <- c(replicate(93, "Longevity"), 
                               replicate(116, "Mass"), 
                               replicate(46, "Metabolic Rate"),
                               replicate(13, "Mortality"),
                               replicate(14, "PO Activity"), 
                               replicate(17, "Reproductive Rate"),
                               replicate(15, "Sex"), 
                               replicate(11, "SOD Activity"),
                               replicate(127, "Survival"), 
                               replicate(26, "Tail Length"),
                               replicate(14, "Triglyceride"))

specific_trait_raw_df_2 <- data.frame("Model" = specific_trait_raw_name_2, 
                                      "Effect" = specific_trait_raw_mean_2)

# Graph code - Part 2

Specific_Trait_Order_2 <- c("Triglyceride", "Tail Length", "Survival", "SOD Activity", "Sex",  
                            "Reproductive Rate", "PO Activity", "Mortality", "Metabolic Rate", 
                            "Mass", "Longevity", "Locomotor Performance")

density_specific_trait_CVR_2 <- specific_trait_table_2 %>% mutate(name = fct_relevel(name, Specific_Trait_Order_2)) %>%
                                ggplot() +
                                geom_density_ridges(data = specific_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Specific_Trait_Order_2)), 
                                                    aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                    scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                geom_linerange(aes(y = rev(seq(1, dim(specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                               size = 1) +
                                geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                size = 1, fatten = 2) +
                                theme_bw() +
                                guides(fill = "none", colour = "none") +
                                labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                      vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7, -0.8, -2.7, 
                                                -2.7, -0.8, -2.7, -2.7))) +
                                theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                theme(axis.ticks = element_blank()) +
                                theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                scale_colour_manual(values = c("#6C87AB", "#6582A9", "#6381A9", "#607C9F", "#5D7AA1", "#5777A1", 
                                                               "#50719C", "#4A6E9C", "#4D6E9A", "#446692", "#3C5F8D")) +
                                scale_fill_manual(values = c("#6C87AB", "#6582A9", "#6381A9", "#607C9F", "#5D7AA1", "#5777A1", 
                                                             "#50719C", "#4A6E9C", "#4D6E9A", "#446692", "#3C5F8D")) +
                                coord_cartesian(xlim = c(-1, 1)) +
                                annotate('text',  x = 1, y = (seq(1, dim(specific_trait_table_2)[1], 1)+0.4),
                                label= paste("italic(k)==", c(specific_trait_table_2["Triglyceride", "K"], 
                                                              specific_trait_table_2["Tail Length", "K"], 
                                                              specific_trait_table_2["Survival", "K"], 
                                                              specific_trait_table_2["SOD Activity", "K"],
                                                              specific_trait_table_2["Sex", "K"],
                                                              specific_trait_table_2["Reproductive Rate", "K"],
                                                              specific_trait_table_2["PO Activity", "K"],
                                                              specific_trait_table_2["Mortality", "K"], 
                                                              specific_trait_table_2["Metabolic Rate", "K"],
                                                              specific_trait_table_2["Mass", "K"],
                                                              specific_trait_table_2["Longevity", "K"]), "~","(", 
                                                            c(specific_trait_table_2["Triglyceride", "group_no"], 
                                                              specific_trait_table_2["Tail Length", "group_no"], 
                                                              specific_trait_table_2["Survival", "group_no"], 
                                                              specific_trait_table_2["SOD Activity", "group_no"],
                                                              specific_trait_table_2["Sex", "group_no"],
                                                              specific_trait_table_2["Reproductive Rate", "group_no"],
                                                              specific_trait_table_2["PO Activity", "group_no"],
                                                              specific_trait_table_2["Mortality", "group_no"], 
                                                              specific_trait_table_2["Metabolic Rate", "group_no"],
                                                              specific_trait_table_2["Mass", "group_no"],
                                                              specific_trait_table_2["Longevity", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Triglyceride", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Tail Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Survival", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["SOD Activity", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Gender", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Reproductive Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["PO Activity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Mortality", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Specific_Trait_Model_CVR_Estimates["Longevity", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                          x = -0.75, y = (seq(1, dim(specific_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_specific_trait_CVR_2 #(400x960)

##### Summary of Overall Plots #####

Overall_Layout <- rbind(c(1, 2), 
                        c(1, 2),
                        c(1, 2),
                        c(1, 2),
                        c(1, 2),
                        c(1, 2),
                        c(1, 2),
                        c(1, 2),
                        c(1, 2),
                        c(1, 2),
                        c(1, 3),
                        c(1, 3),
                        c(1, 3),
                        c(1, 3),
                        c(1, 3),
                        c(1, 3),
                        c(1, 3), 
                        c(1, 4), 
                        c(1, 4),
                        c(1, 4),
                        c(1, 4))

Overall_Combined_CVR <- grid.arrange(density_specific_trait_CVR, density_class_CVR, density_trait_CVR, density_fluctuation_CVR, 
                                     layout_matrix = Overall_Layout)

Overall_Combined_CVR #(850 x 1500 - does not include amplitude plot)

##### Population-Level Subset Model - CVR #####
Population_Subset_Data <- data %>% filter(Trait_Category == "Population")
Population_Species <- Population_Subset_Data %>% select("phylo") %>% unique()

Population_A_cor <- as.data.frame(A_cor)
Population_A_cor <- Population_A_cor[c(Population_Species$phylo), c(Population_Species$phylo)]
Population_A_cor <- as.matrix(Population_A_cor)

Population_VCV_InCVR <- make_VCV_matrix(Population_Subset_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Population_Model_CVR <- metafor::rma.mv(InCVR ~ 1, V = Population_VCV_InCVR, test = "t", dfs = "contain",
                                            random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                          ~1|Shared_Animal_Number, ~1|Measurement), 
                                            R = list(phylo=Population_A_cor), data = Population_Subset_Data, method = "REML", sparse = TRUE, 
                                            control=list(rel.tol=1e-9))
    saveRDS(Population_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Population_Model_CVR.rds")
  } else {
            Population_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Population_Model_CVR.rds")})

Population_Model_CVR_rob <- robust(Population_Model_CVR, cluster = Population_Subset_Data$Study_ID, adjust = TRUE)

Population_Model_CVR_Estimates <- data.frame(estimate = Population_Model_CVR$b, 
                                             ci.lb = Population_Model_CVR$ci.lb, 
                                             ci.ub = Population_Model_CVR$ci.ub)
Population_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Population_Model_CVR), 2))

#### Population-Level Subset Model - Fluctuation Amplitude Meta-Regression - CVR ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Population_Amplitude_Model_CVR <- metafor::rma.mv(InCVR, V = Population_VCV_InCVR, test = "t", dfs = "contain",
                                                      mods = ~ T2_Magnitude - 1,
                                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                    ~1|Shared_Animal_Number, ~1|Measurement), 
                                                      R = list(phylo=Population_A_cor), data = Population_Subset_Data, method = "REML", sparse = TRUE, 
                                                      control=list(rel.tol=1e-9))
    saveRDS(Population_Amplitude_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Population_Amplitude_Model_CVR.rds")
  } else {
            Population_Amplitude_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Population_Amplitude_Model_CVR.rds")})

Population_Amplitude_Model_CVR_rob <- robust(Population_Amplitude_Model_CVR, cluster = Population_Subset_Data$Study_ID, adjust = TRUE)

Population_Amplitude_Model_CVR_Estimates <- data.frame(estimate = Population_Amplitude_Model_CVR$b, 
                                                       ci.lb = Population_Amplitude_Model_CVR$ci.lb, 
                                                       ci.ub = Population_Amplitude_Model_CVR$ci.ub)
Population_Amplitude_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Population_Amplitude_Model_CVR), 2))

# Preparing Graph

Population_Plot_Data <- Population_Subset_Data
Population_Plot_Data <- Population_Plot_Data %>% mutate(n_category = ifelse(n_R1.1 <= 25, "25", 
                                                                     ifelse(n_R1.1 > 25 & n_R1.1 <= 50, "50", 
                                                                     ifelse(n_R1.1 > 50 & n_R1.1 <= 75, "75", "> 75"))))

# Graph Code

Population_Amplitude_Plot <- ggplot(Population_Plot_Data, aes(x = T2_Magnitude, y = InCVR)) + 
                             geom_point(aes(x = T2_Magnitude, y = InCVR, 
                             size = fct_relevel(n_category, c("25", "50", "75", "> 75"))), 
                             shape = 21, fill = "#4292c6", alpha = 0.5) + 
                             labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (lnCVR)", 
                             size = "Sample Size", title = "Population-level Traits") +
                             theme_bw() +
                             theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                             theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                             theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                             theme(legend.position = "bottom", legend.direction = "horizontal") + 
                             geom_hline(yintercept = Population_Amplitude_Model_CVR_Estimates$estimate, lty = 2) + 
                             geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                             stat_poly_eq(formula = y ~ x, 
                                          aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                          parse = TRUE) +
                             coord_cartesian(xlim = c(0, 30), 
                                             ylim = c(-2.5, 2.5))

Population_Amplitude_Plot #(400x400)

#### Population-Level Subset Model - Type of Fluctuation Meta-Regression - CVR ####
Population_Fluctuation_Data <- Population_Subset_Data %>% filter(!is.na(Fluctuation_Category))
Population_Fluctuation_Exploration <- Population_Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Population_Fluctuation_Exploration) <- Population_Fluctuation_Exploration$Fluctuation_Category

Population_Fluctuation_Data <- Population_Subset_Data %>% filter(Fluctuation_Category != "Stochastic")

Population_Fluctuation_Species_Count <- Population_Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% 
                                        table() %>% data.frame() %>% 
                                        filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Population_Fluctuation_Species_Count) <- Population_Fluctuation_Species_Count$Fluctuation_Category

Population_Fluctuation_Study_Count <- Population_Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% 
                                      table() %>% data.frame() %>% 
                                      filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Population_Fluctuation_Study_Count) <- Population_Fluctuation_Study_Count$Fluctuation_Category

Population_Fluctuation_Species <- Population_Fluctuation_Data %>% select("phylo") %>% unique()

Population_Fluctuation_A_cor <- as.data.frame(A_cor)
Population_Fluctuation_A_cor <- Population_Fluctuation_A_cor[c(Population_Fluctuation_Species$phylo), c(Population_Fluctuation_Species$phylo)]
Population_Fluctuation_A_cor <- as.matrix(Population_Fluctuation_A_cor)

Population_Fluctuation_VCV_InCVR <- make_VCV_matrix(Population_Fluctuation_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Population_Fluctuation_Model_CVR <- metafor::rma.mv(InCVR, V = Population_Fluctuation_VCV_InCVR, test = "t", dfs = "contain",
                                                        mods = ~ Fluctuation_Category - 1,
                                                        random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                      ~1|Shared_Animal_Number, ~1|Measurement), 
                                                        R = list(phylo=Population_Fluctuation_A_cor), data = Population_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                        control=list(rel.tol=1e-9))
    saveRDS(Population_Fluctuation_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Population_Fluctuation_Model_CVR.rds")
  } else {
            Population_Fluctuation_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Population_Fluctuation_Model_CVR.rds")})

Population_Fluctuation_Model_CVR_rob <- robust(Population_Fluctuation_Model_CVR, cluster = Population_Fluctuation_Data$Study_ID, adjust = TRUE)

Population_Fluctuation_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Population_Fluctuation_Model_CVR$b), 21, 100),
                                                         estimate = Population_Fluctuation_Model_CVR$b, 
                                                         ci.lb = Population_Fluctuation_Model_CVR$ci.lb, 
                                                         ci.ub = Population_Fluctuation_Model_CVR$ci.ub)
rownames(Population_Fluctuation_Model_CVR_Estimates) <- Population_Fluctuation_Model_CVR_Estimates$Category
Population_Fluctuation_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Population_Fluctuation_Model_CVR), 2))

# Preparing Graph - Combined

population_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise")

population_fluctuation_k <- data.frame("k" = c(Population_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                               Population_Fluctuation_Exploration["Alternating", "Freq"], 
                                               Population_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                       row.names = population_fluctuation_rnames)

population_fluctuation_group_no <- data.frame("Spp No." = c(Population_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                            Population_Fluctuation_Species_Count["Alternating", "Freq"], 
                                                            Population_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                              row.names = population_fluctuation_rnames)

population_fluctuation_study <- data.frame("Study" = c(Population_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                       Population_Fluctuation_Study_Count["Alternating", "Freq"], 
                                                       Population_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                              row.names = population_fluctuation_rnames)

Population_Fluctuation_Model_CVR_Estimates_Reorder <- Population_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise"), ]

population_fluctuation_table <- data.frame(estimate = Population_Fluctuation_Model_CVR_Estimates_Reorder[,"estimate"], 
                                           lowerCL = Population_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.lb"], 
                                           upperCL = Population_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.ub"], 
                                           K = population_fluctuation_k[,1], 
                                           group_no = population_fluctuation_group_no[,1], 
                                           row.names = population_fluctuation_rnames)
population_fluctuation_table$name <- row.names(population_fluctuation_table)

population_fluctuation_raw_mean <- c(unlist(unname(Population_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Population_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Population_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                     select("InCVR"))))

population_fluctuation_raw_name <- c(replicate(60, "Sinusoidal (Sine Curve)"), 
                                     replicate(74, "Alternating"), 
                                     replicate(22, "Stepwise"))

population_fluctuation_raw_df <- data.frame("Model" = population_fluctuation_raw_name, 
                                            "Effect" = population_fluctuation_raw_mean)

# Graph code - Combined

Population_Fluctuation_Order <- c("Stepwise", "Alternating", "Sinusoidal (Sine Curve)")

density_population_fluctuation_CVR <- population_fluctuation_table %>% mutate(name = fct_relevel(name, Population_Fluctuation_Order)) %>%
                                      ggplot() +
                                      geom_density_ridges(data = population_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Population_Fluctuation_Order)), 
                                                          aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                          scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                      geom_linerange(aes(y = rev(seq(1, dim(population_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                     size = 1) +
                                      geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(population_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                      size = 1, fatten = 2) +
                                      theme_bw() +
                                      guides(fill = "none", colour = "none") +
                                      labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                      theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                       vjust = c(-2.7, -2.7, -0.8))) +
                                      theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                      theme(axis.ticks = element_blank()) +
                                      theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                      scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                      scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                      coord_cartesian(xlim = c(-1, 1)) +
                                      annotate('text',  x = 1, y = (seq(1, dim(population_fluctuation_table)[1], 1)+0.4),
                                      label= paste("italic(k)==", c(population_fluctuation_table["Stepwise", "K"], 
                                                                    population_fluctuation_table["Alternating", "K"], 
                                                                    population_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                  c(population_fluctuation_table["Stepwise", "group_no"], 
                                                                    population_fluctuation_table["Alternating", "group_no"], 
                                                                    population_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                                   ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                      geom_label(aes(label=c(paste(format(round(mean(exp(Population_Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                             paste(format(round(mean(exp(Population_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                             paste(format(round(mean(exp(Population_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                 x = -0.75, y = (seq(1, dim(population_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_population_fluctuation_CVR #(400x320)

# Preparing Graph - Part 1

population_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)", "Alternating")

population_fluctuation_k_1 <- data.frame("k" = c(Population_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                 Population_Fluctuation_Exploration["Alternating", "Freq"]), 
                                         row.names = population_fluctuation_rnames_1)

population_fluctuation_group_no_1 <- data.frame("Spp No." = c(Population_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                              Population_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                                row.names = population_fluctuation_rnames_1)

population_fluctuation_study_1 <- data.frame("Study" = c(Population_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                         Population_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                             row.names = population_fluctuation_rnames_1)

Population_Fluctuation_Model_CVR_Estimates_Reorder_1 <- Population_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

population_fluctuation_table_1 <- data.frame(estimate = Population_Fluctuation_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                             lowerCL = Population_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                             upperCL = Population_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                             K = population_fluctuation_k_1[,1], 
                                             group_no = population_fluctuation_group_no_1[,1], 
                                             row.names = population_fluctuation_rnames_1)
population_fluctuation_table_1$name <- row.names(population_fluctuation_table_1)

population_fluctuation_raw_mean_1 <- c(unlist(unname(Population_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                     select("InCVR"))), 
                                       unlist(unname(Population_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                     select("InCVR"))))

population_fluctuation_raw_name_1 <- c(replicate(60, "Sinusoidal (Sine Curve)"), 
                                       replicate(74, "Alternating"))

population_fluctuation_raw_df_1 <- data.frame("Model" = population_fluctuation_raw_name_1, 
                                              "Effect" = population_fluctuation_raw_mean_1)

# Graph code - Part 1

Population_Fluctuation_Order_1 <- c("Alternating", "Sinusoidal (Sine Curve)")

density_population_fluctuation_CVR_1 <- population_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Population_Fluctuation_Order_1)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = population_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Population_Fluctuation_Order_1)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(population_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(population_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                        size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                              vjust = c(-2.7, -0.8))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                        scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                        coord_cartesian(xlim = c(-1, 1)) +
                                        annotate('text',  x = 1, y = (seq(1, dim(population_fluctuation_table_1)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(population_fluctuation_table_1["Alternating", "K"], 
                                                                      population_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                    c(population_fluctuation_table_1["Alternating", "group_no"], 
                                                                      population_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Population_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                               paste(format(round(mean(exp(Population_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.75, y = (seq(1, dim(population_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_population_fluctuation_CVR_1 #(400x240)

# Preparing Graph - Part 2

population_fluctuation_rnames_2 <- c("Stepwise")

population_fluctuation_k_2 <- data.frame("k" = c(Population_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                         row.names = population_fluctuation_rnames_2)

population_fluctuation_group_no_2 <- data.frame("Spp No." = c(Population_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                                row.names = population_fluctuation_rnames_2)

population_fluctuation_study_2 <- data.frame("Study" = c(Population_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                             row.names = population_fluctuation_rnames_2)

Population_Fluctuation_Model_CVR_Estimates_Reorder_2 <- Population_Fluctuation_Model_CVR_Estimates[c("Stepwise"), ]

population_fluctuation_table_2 <- data.frame(estimate = Population_Fluctuation_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                             lowerCL = Population_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                             upperCL = Population_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                             K = population_fluctuation_k_2[,1], 
                                             group_no = population_fluctuation_group_no_2[,1], 
                                             row.names = population_fluctuation_rnames_2)
population_fluctuation_table_2$name <- row.names(population_fluctuation_table_2)

population_fluctuation_raw_mean_2 <- c(unlist(unname(Population_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                     select("InCVR"))))

population_fluctuation_raw_name_2 <- c(replicate(22, "Stepwise"))

population_fluctuation_raw_df_2 <- data.frame("Model" = population_fluctuation_raw_name_2, 
                                              "Effect" = population_fluctuation_raw_mean_2)

# Graph code - Part 2

Population_Fluctuation_Order_2 <- c("Stepwise")

density_population_fluctuation_CVR_2 <- population_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Population_Fluctuation_Order_2)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = population_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Population_Fluctuation_Order_2)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.2, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(population_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_linerange(aes(y = rev(seq(1, dim(population_fluctuation_table_2)[1], 1)), xmin = min(population_fluctuation_raw_df_2$Effect)-0.1, xmax = -1.5, colour = name),
                                                       size = 1) +
                                        geom_linerange(aes(y = rev(seq(1, dim(population_fluctuation_table_2)[1], 1)), xmin = max(population_fluctuation_raw_df_2$Effect)+0.1, xmax = 1.5, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(population_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                        size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                              vjust = c(-2.7))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#5D7AA1")) +
                                        scale_fill_manual(values = c("#5D7AA1")) +
                                        coord_cartesian(xlim = c(-1, 1)) +
                                        annotate('text',  x = 1, y = (seq(1, dim(population_fluctuation_table_2)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(population_fluctuation_table_2["Stepwise", "K"]), "~","(", 
                                                                    c(population_fluctuation_table_2["Stepwise", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Population_Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                  x = -0.75, y = (seq(1, dim(population_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_population_fluctuation_CVR_2 #(400x160)

##### Population-Level Subset Model - Class Meta-Regression - CVR #####
Population_Class_Exploration <- Population_Subset_Data %>% select("Class") %>% table() %>% data.frame()
rownames(Population_Class_Exploration) <- Population_Class_Exploration$Class

Population_Class_Data <- Population_Subset_Data %>% filter(Class == "Actinopteri"|
                                                           Class == "Arachnida"|
                                                           Class == "Insecta")

Population_Class_Species_Count <- Population_Class_Data %>% select("Scientific_Name", "Class") %>% table() %>% data.frame() %>%
                                  filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Population_Class_Species_Count) <- Population_Class_Species_Count$Class

Population_Class_Study_Count <- Population_Class_Data %>% select("Study_ID", "Class") %>% table() %>% data.frame() %>%
                                filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Population_Class_Study_Count) <- Population_Class_Study_Count$Class

Population_Class_Species <- Population_Class_Data %>% select("phylo") %>% unique()

Population_Class_A_cor <- as.data.frame(A_cor)
Population_Class_A_cor <- Population_Class_A_cor[c(Population_Class_Species$phylo), c(Population_Class_Species$phylo)]
Population_Class_A_cor <- as.matrix(Population_Class_A_cor)

Population_Class_VCV_InCVR <- make_VCV_matrix(Population_Class_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Population_Class_Model_CVR <- metafor::rma.mv(InCVR, V = Population_Class_VCV_InCVR, test = "t", dfs = "contain",
                                                  mods = ~ Class - 1,
                                                  random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                ~1|Shared_Animal_Number, ~1|Measurement), 
                                                  R = list(phylo=Population_Class_A_cor), data = Population_Class_Data, method = "REML", sparse = TRUE, 
                                                  control=list(rel.tol=1e-9))
    saveRDS(Population_Class_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Population_Class_Model_CVR.rds")
  } else {
            Population_Class_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Population_Class_Model_CVR.rds")})

Population_Class_Model_CVR_rob <- robust(Population_Class_Model_CVR, cluster = Population_Class_Data$Study_ID, adjust = TRUE)

Population_Class_Model_CVR_Estimates <- data.frame(Class = substr(row.names(Population_Class_Model_CVR$b), 6, 100),
                                                   estimate = Population_Class_Model_CVR$b, 
                                                   ci.lb = Population_Class_Model_CVR$ci.lb, 
                                                   ci.ub = Population_Class_Model_CVR$ci.ub)
rownames(Population_Class_Model_CVR_Estimates) <- Population_Class_Model_CVR_Estimates$Class
Population_Class_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Population_Class_Model_CVR), 2))

# Preparing Graph - Combined

population_class_rnames <- c("Actinopteri", "Arachnida", "Insecta")

population_class_k <- data.frame("k" = c(Population_Class_Exploration["Actinopteri", "Freq"], 
                                         Population_Class_Exploration["Arachnida", "Freq"], 
                                         Population_Class_Exploration["Insecta", "Freq"]), 
                                 row.names = population_class_rnames)

population_class_group_no <- data.frame("Spp No." = c(Population_Class_Species_Count["Actinopteri", "Freq"], 
                                                      Population_Class_Species_Count["Arachnida", "Freq"],
                                                      Population_Class_Species_Count["Insecta", "Freq"]), 
                                        row.names = population_class_rnames)

population_class_study <- data.frame("Study" = c(Population_Class_Study_Count["Actinopteri", "Freq"], 
                                                 Population_Class_Study_Count["Arachnida", "Freq"],
                                                 Population_Class_Study_Count["Insecta", "Freq"]), 
                                        row.names = population_class_rnames)

population_class_table <- data.frame(estimate = Population_Class_Model_CVR_Estimates[,"estimate"], 
                                     lowerCL = Population_Class_Model_CVR_Estimates[,"ci.lb"], 
                                     upperCL = Population_Class_Model_CVR_Estimates[,"ci.ub"], 
                                     K = population_class_k[,1], 
                                     group_no = population_class_group_no[,1], 
                                     row.names = population_class_rnames)
population_class_table$name <- row.names(population_class_table)

population_class_raw_mean <- c(unlist(unname(Population_Class_Data %>% filter(`Class` == "Actinopteri") %>% 
                                               select("InCVR"))),
                               unlist(unname(Population_Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                               select("InCVR"))),
                               unlist(unname(Population_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                               select("InCVR"))))

population_class_raw_name <- c(replicate(19, "Actinopteri"), 
                               replicate(24, "Arachnida"),
                               replicate(108, "Insecta"))

population_class_raw_df <- data.frame("Model" = population_class_raw_name, 
                                      "Effect" = population_class_raw_mean)

# Graph code - Combined

Population_Class_Order <- c("Insecta", "Arachnida", "Actinopteri")

density_population_class_CVR <- population_class_table %>% mutate(name = fct_relevel(name, Population_Class_Order)) %>%
                                ggplot() +
                                geom_density_ridges(data = population_class_raw_df %>% mutate(Model = fct_relevel(Model, Population_Class_Order)), 
                                                    aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                    scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                geom_linerange(aes(y = rev(seq(1, dim(population_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                               size = 1) +
                                geom_linerange(aes(y = rev(seq(1, dim(population_class_table)[1], 1)), xmin = min(population_class_raw_df$Effect)-0.11, xmax = -1.5, colour = name),
                                               size = 1) +
                                geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(population_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                size = 1, fatten = 2) +
                                theme_bw() +
                                guides(fill = "none", colour = "none") +
                                labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                 vjust = c(-2.7, -2.7, -2.7))) +
                                theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                theme(axis.ticks = element_blank()) +
                                theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                coord_cartesian(xlim = c(-1, 1)) +
                                annotate('text',  x = 1, y = (seq(1, dim(population_class_table)[1], 1)+0.4),
                                label= paste("italic(k)==", c(population_class_table["Insecta", "K"],
                                                              population_class_table["Arachnida", "K"],
                                                              population_class_table["Actinopteri", "K"]), "~","(", 
                                                            c(population_class_table["Insecta", "group_no"],
                                                              population_class_table["Arachnida", "group_no"],
                                                              population_class_table["Actinopteri", "group_no"]), 
                                             ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(Population_Class_Model_CVR_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Population_Class_Model_CVR_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Population_Class_Model_CVR_Estimates["Actinopteri", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                           x = -0.75, y = (seq(1, dim(population_class_table)[1], 1)+0.4)), size = 3.5)

density_population_class_CVR #(400x320)

# Preparing Graph - Part 1

population_class_rnames_1 <- c("Actinopteri", "Arachnida")

population_class_k_1 <- data.frame("k" = c(Population_Class_Exploration["Actinopteri", "Freq"], 
                                           Population_Class_Exploration["Arachnida", "Freq"]), 
                                   row.names = population_class_rnames_1)

population_class_group_no_1 <- data.frame("Spp No." = c(Population_Class_Species_Count["Actinopteri", "Freq"], 
                                                        Population_Class_Species_Count["Arachnida", "Freq"]), 
                                          row.names = population_class_rnames_1)

population_class_study_1 <- data.frame("Study" = c(Population_Class_Study_Count["Actinopteri", "Freq"], 
                                                   Population_Class_Study_Count["Arachnida", "Freq"]), 
                                       row.names = population_class_rnames_1)

Population_Class_Model_CVR_Estimates_Reorder_1 <- Population_Class_Model_CVR_Estimates[c("Actinopteri", "Arachnida"), ]

population_class_table_1 <- data.frame(estimate = Population_Class_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                       lowerCL = Population_Class_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                       upperCL = Population_Class_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                       K = population_class_k_1[,1], 
                                       group_no = population_class_group_no_1[,1], 
                                       row.names = population_class_rnames_1)
population_class_table_1$name <- row.names(population_class_table_1)

population_class_raw_mean_1 <- c(unlist(unname(Population_Class_Data %>% filter(`Class` == "Actinopteri") %>% 
                                               select("InCVR"))),
                                 unlist(unname(Population_Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                               select("InCVR"))))

population_class_raw_name_1 <- c(replicate(19, "Actinopteri"), 
                                 replicate(24, "Arachnida"))

population_class_raw_df_1 <- data.frame("Model" = population_class_raw_name_1, 
                                        "Effect" = population_class_raw_mean_1)

# Graph code - Part 1

Population_Class_Order_1 <- c("Arachnida", "Actinopteri")

density_population_class_CVR_1 <- population_class_table_1 %>% mutate(name = fct_relevel(name, Population_Class_Order_1)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = population_class_raw_df_1 %>% mutate(Model = fct_relevel(Model, Population_Class_Order_1)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(population_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(population_class_table_1)[1], 1)), xmin = max(population_class_raw_df_1$Effect)+0.1, xmax = 1.5, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(population_class_table_1)[1], 1)), xmin = min(population_class_raw_df_1$Effect)-0.07, xmax = -1.5, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(population_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                  size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                        vjust = c(-2.7, -2.7))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                  scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                  coord_cartesian(xlim = c(-1, 1)) +
                                  annotate('text',  x = 1, y = (seq(1, dim(population_class_table_1)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(population_class_table_1["Arachnida", "K"],
                                                                population_class_table_1["Actinopteri", "K"]), "~","(", 
                                                              c(population_class_table_1["Arachnida", "group_no"],
                                                                population_class_table_1["Actinopteri", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Population_Class_Model_CVR_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(Population_Class_Model_CVR_Estimates["Actinopteri", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.75, y = (seq(1, dim(population_class_table_1)[1], 1)+0.4)), size = 3.5)

density_population_class_CVR_1 #(400x240)

# Preparing Graph - Part 2

population_class_rnames_2 <- c("Insecta")

population_class_k_2 <- data.frame("k" = c(Population_Class_Exploration["Insecta", "Freq"]), 
                                   row.names = population_class_rnames_2)

population_class_group_no_2 <- data.frame("Spp No." = c(Population_Class_Species_Count["Insecta", "Freq"]), 
                                          row.names = population_class_rnames_2)

population_class_study_2 <- data.frame("Study" = c(Population_Class_Study_Count["Insecta", "Freq"]), 
                                       row.names = population_class_rnames_2)

Population_Class_Model_CVR_Estimates_Reorder_2 <- Population_Class_Model_CVR_Estimates[c("Insecta"), ]

population_class_table_2 <- data.frame(estimate = Population_Class_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                       lowerCL = Population_Class_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                       upperCL = Population_Class_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                       K = population_class_k_2[,1], 
                                       group_no = population_class_group_no_2[,1], 
                                       row.names = population_class_rnames_2)
population_class_table_2$name <- row.names(population_class_table_2)

population_class_raw_mean_2 <- c(unlist(unname(Population_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                               select("InCVR"))))

population_class_raw_name_2 <- c(replicate(108, "Insecta"))

population_class_raw_df_2 <- data.frame("Model" = population_class_raw_name_2, 
                                        "Effect" = population_class_raw_mean_2)

# Graph code - Part 2

Population_Class_Order_2 <- c("Insecta")

density_population_class_CVR_2 <- population_class_table_2 %>% mutate(name = fct_relevel(name, Population_Class_Order_2)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = population_class_raw_df_2 %>% mutate(Model = fct_relevel(Model, Population_Class_Order_2)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.1, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(population_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                                  geom_linerange(aes(y = rev(seq(1, dim(population_class_table_2)[1], 1)), xmin = min(population_class_raw_df_2$Effect)-0.2, xmax = -1.5, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(population_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                  size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                        vjust = c(-2.7))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#5D7AA1")) +
                                  scale_fill_manual(values = c("#5D7AA1")) +
                                  coord_cartesian(xlim = c(-1, 1)) +
                                  annotate('text',  x = 1, y = (seq(1, dim(population_class_table_2)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(population_class_table_2["Insecta", "K"]), "~","(", 
                                                              c(population_class_table_2["Insecta", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Population_Class_Model_CVR_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.75, y = (seq(1, dim(population_class_table_2)[1], 1)+0.4)), size = 3.5)

density_population_class_CVR_2 #(400x160)

# Combining Population Plots

Population_Layout <- rbind(c(1, 2))

Population_Combined_CVR <- grid.arrange(density_population_class_CVR, density_population_fluctuation_CVR, 
                                        layout_matrix = Population_Layout)

Population_Combined_CVR #(850 x 300 - does not include amplitude plot)

##### Individual-Level Trait Subset Model - CVR #####
Individual_Subset_Data <- data %>% filter(Trait_Category != "Population")
Individual_Species <- Individual_Subset_Data %>% select("phylo") %>% unique()

Individual_A_cor <- as.data.frame(A_cor)
Individual_A_cor <- Individual_A_cor[c(Individual_Species$phylo), c(Individual_Species$phylo)]
Individual_A_cor <- as.matrix(Individual_A_cor)

Individual_VCV_InCVR <- make_VCV_matrix(Individual_Subset_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  31ish minutes
  if(run){
    Individual_Model_CVR <- metafor::rma.mv(InCVR ~ 1, V = Individual_VCV_InCVR, test = "t", dfs = "contain",
                                            random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                          ~1|Shared_Animal_Number, ~1|Measurement), 
                                            R = list(phylo=Individual_A_cor), data = Individual_Subset_Data, method = "REML", sparse = TRUE,
                                            control=list(rel.tol=1e-9))
    saveRDS(Individual_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Individual_Model_CVR.rds")
  } else {
            Individual_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Individual_Model_CVR.rds")})

Individual_Model_CVR_rob <- robust(Individual_Model_CVR, cluster = Individual_Subset_Data$Study_ID, adjust = TRUE)

Individual_Model_CVR_Estimates <- data.frame(estimate = Individual_Model_CVR$b, 
                                             ci.lb = Individual_Model_CVR$ci.lb, 
                                             ci.ub = Individual_Model_CVR$ci.ub)
Individual_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Individual_Model_CVR), 2))

#### Individual-Level Subset Model - Fluctuation Amplitude Meta-Regression - CVR ####
run <- FALSE
system.time( #  19ish minutes
  if(run){
    Individual_Amplitude_Model_CVR <- metafor::rma.mv(InCVR, V = Individual_VCV_InCVR, test = "t", dfs = "contain",
                                                      mods = ~ T2_Magnitude - 1,
                                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                    ~1|Shared_Animal_Number, ~1|Measurement), 
                                                      R = list(phylo=Individual_A_cor), data = Individual_Subset_Data, method = "REML", sparse = TRUE, 
                                                      control=list(rel.tol=1e-9))
    saveRDS(Individual_Amplitude_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Individual_Amplitude_Model_CVR.rds")
  } else {
            Individual_Amplitude_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Individual_Amplitude_Model_CVR.rds")})

Individual_Amplitude_Model_CVR_rob <- robust(Individual_Amplitude_Model_CVR, cluster = Individual_Subset_Data$Study_ID, adjust = TRUE)

Individual_Amplitude_Model_CVR_Estimates <- data.frame(estimate = Individual_Amplitude_Model_CVR$b, 
                                                       ci.lb = Individual_Amplitude_Model_CVR$ci.lb, 
                                                       ci.ub = Individual_Amplitude_Model_CVR$ci.ub)
Individual_Amplitude_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Individual_Amplitude_Model_CVR), 2))

# Preparing Graph

Individual_Plot_Data <- Individual_Subset_Data
Individual_Plot_Data <- Individual_Plot_Data %>% mutate(n_category = ifelse(n_R1.1 <= 25, "25", 
                                                                     ifelse(n_R1.1 > 25 & n_R1.1 <= 50, "50", 
                                                                     ifelse(n_R1.1 > 50 & n_R1.1 <= 75, "75", "> 75"))))

# Graph Code

Individual_Amplitude_Plot_CVR <- ggplot(Individual_Plot_Data, aes(x = T2_Magnitude, y = InCVR)) + 
                                 geom_point(aes(x = T2_Magnitude, y = InCVR, 
                                            size = fct_relevel(n_category, c("25", "50", "75", "> 75"))), 
                                            shape = 21, fill = "#4292c6", alpha = 0.5) + 
                                 labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (lnCVR)", 
                                      size = "Sample Size", title = "Individual-level Traits") +
                                 theme_bw() +
                                 theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                                 theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                                 theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                                 theme(legend.position = "bottom", legend.direction = "horizontal") + 
                                 geom_hline(yintercept = Individual_Model_CVR_Estimates$estimate, lty = 2) + 
                                 geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                                 stat_poly_eq(formula = y ~ x, 
                                              aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                              parse = TRUE) +
                                 coord_cartesian(xlim = c(0, 30), 
                                                 ylim = c(-2.5, 2.5))

Individual_Amplitude_Plot_CVR #(400x400)

#### Individual-Level Subset Model - Type of Fluctuation Meta-Regression - CVR ####
Individual_Fluctuation_Data <- Individual_Subset_Data %>% filter(!is.na(Fluctuation_Category))

Individual_Fluctuation_Exploration <- Individual_Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Individual_Fluctuation_Exploration) <- Individual_Fluctuation_Exploration$Fluctuation_Category

Individual_Fluctuation_Species_Count <- Individual_Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                        filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Individual_Fluctuation_Species_Count) <- Individual_Fluctuation_Species_Count$Fluctuation_Category

Individual_Fluctuation_Study_Count <- Individual_Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                      filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Individual_Fluctuation_Study_Count) <- Individual_Fluctuation_Study_Count$Fluctuation_Category

Individual_Fluctuation_Species <- Individual_Fluctuation_Data %>% select("phylo") %>% unique()

Individual_Fluctuation_A_cor <- as.data.frame(A_cor)
Individual_Fluctuation_A_cor <- Individual_Fluctuation_A_cor[c(Individual_Fluctuation_Species$phylo), c(Individual_Fluctuation_Species$phylo)]
Individual_Fluctuation_A_cor <- as.matrix(Individual_Fluctuation_A_cor)

Individual_Fluctuation_VCV_InCVR <- make_VCV_matrix(Individual_Fluctuation_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  20ish minutes
  if(run){
    Individual_Fluctuation_Model_CVR <- metafor::rma.mv(InCVR, V = Individual_Fluctuation_VCV_InCVR, test = "t", dfs = "contain",
                                                        mods = ~ Fluctuation_Category - 1,
                                                        random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                      ~1|Shared_Animal_Number, ~1|Measurement), 
                                                        R = list(phylo=Individual_Fluctuation_A_cor), data = Individual_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                        control=list(rel.tol=1e-9))
    saveRDS(Individual_Fluctuation_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Individual_Fluctuation_Model_CVR.rds")
  } else {
            Individual_Fluctuation_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Individual_Fluctuation_Model_CVR.rds")})

Individual_Fluctuation_Model_CVR_rob <- robust(Individual_Fluctuation_Model_CVR, cluster = Individual_Fluctuation_Data$Study_ID, adjust = TRUE)

Individual_Fluctuation_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Individual_Fluctuation_Model_CVR$b), 21, 100),
                                                         estimate = Individual_Fluctuation_Model_CVR$b, 
                                                         ci.lb = Individual_Fluctuation_Model_CVR$ci.lb, 
                                                         ci.ub = Individual_Fluctuation_Model_CVR$ci.ub)
rownames(Individual_Fluctuation_Model_CVR_Estimates) <- Individual_Fluctuation_Model_CVR_Estimates$Category
Individual_Fluctuation_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Individual_Fluctuation_Model_CVR), 2))


# Preparing Graph - Combined

individual_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise", "Stochastic")

individual_fluctuation_k <- data.frame("k" = c(Individual_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                               Individual_Fluctuation_Exploration["Alternating", "Freq"], 
                                               Individual_Fluctuation_Exploration["Stepwise", "Freq"], 
                                               Individual_Fluctuation_Exploration["Stochastic", "Freq"]), 
                                       row.names = individual_fluctuation_rnames)

individual_fluctuation_group_no <- data.frame("Spp No." = c(Individual_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                            Individual_Fluctuation_Species_Count["Alternating", "Freq"], 
                                                            Individual_Fluctuation_Species_Count["Stepwise", "Freq"], 
                                                            Individual_Fluctuation_Species_Count["Stochastic", "Freq"]), 
                                              row.names = individual_fluctuation_rnames)

individual_fluctuation_study <- data.frame("Study" = c(Individual_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                       Individual_Fluctuation_Study_Count["Alternating", "Freq"], 
                                                       Individual_Fluctuation_Study_Count["Stepwise", "Freq"], 
                                                       Individual_Fluctuation_Study_Count["Stochastic", "Freq"]), 
                                              row.names = individual_fluctuation_rnames)

Individual_Fluctuation_Model_CVR_Estimates_Reorder <- Individual_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise", "Stochastic"), ]

individual_fluctuation_table <- data.frame(estimate = Individual_Fluctuation_Model_CVR_Estimates_Reorder[,"estimate"], 
                                           lowerCL = Individual_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.lb"], 
                                           upperCL = Individual_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.ub"], 
                                           K = individual_fluctuation_k[,1], 
                                           group_no = individual_fluctuation_group_no[,1], 
                                           row.names = individual_fluctuation_rnames)
individual_fluctuation_table$name <- row.names(individual_fluctuation_table)

individual_fluctuation_raw_mean <- c(unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stochastic") %>% 
                                                     select("InCVR"))))

individual_fluctuation_raw_name <- c(replicate(481, "Sinusoidal (Sine Curve)"), 
                                     replicate(533, "Alternating"), 
                                     replicate(196, "Stepwise"), 
                                     replicate(20, "Stochastic"))

individual_fluctuation_raw_df <- data.frame("Model" = individual_fluctuation_raw_name, 
                                            "Effect" = individual_fluctuation_raw_mean)

# Graph code - Combined

Individual_Fluctuation_Order <- c("Stochastic", "Stepwise", 
                                  "Alternating", "Sinusoidal (Sine Curve)")

density_individual_fluctuation_CVR <- individual_fluctuation_table %>% mutate(name = fct_relevel(name, Individual_Fluctuation_Order)) %>%
                                      ggplot() +
                                      geom_density_ridges(data = individual_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Individual_Fluctuation_Order)), 
                                                          aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                          scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                      geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                     size = 1) +
                                      geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                      size = 1, fatten = 2) +
                                      theme_bw() +
                                      guides(fill = "none", colour = "none") +
                                      labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                      theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                       vjust = c(-2.7, -2.7, -2.7, -0.8))) +
                                      theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                      theme(axis.ticks = element_blank()) +
                                      theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                      scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                      scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                      coord_cartesian(xlim = c(-1, 1)) +
                                      annotate('text',  x = 1, y = (seq(1, dim(individual_fluctuation_table)[1], 1)+0.4),
                                      label= paste("italic(k)==", c(individual_fluctuation_table["Stochastic", "K"], 
                                                                    individual_fluctuation_table["Stepwise", "K"], 
                                                                    individual_fluctuation_table["Alternating", "K"], 
                                                                    individual_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                  c(individual_fluctuation_table["Stochastic", "group_no"], 
                                                                    individual_fluctuation_table["Stepwise", "group_no"], 
                                                                    individual_fluctuation_table["Alternating", "group_no"], 
                                                                    individual_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                                   ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                      geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Fluctuation_Model_CVR_Estimates["Stochastic", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                             paste(format(round(mean(exp(Individual_Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                             paste(format(round(mean(exp(Individual_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                             paste(format(round(mean(exp(Individual_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                 x = -0.75, y = (seq(1, dim(individual_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_individual_fluctuation_CVR #(400x400)

# Preparing Graph - Part 1

individual_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)", "Alternating")

individual_fluctuation_k_1 <- data.frame("k" = c(Individual_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                 Individual_Fluctuation_Exploration["Alternating", "Freq"]), 
                                         row.names = individual_fluctuation_rnames_1)

individual_fluctuation_group_no_1 <- data.frame("Spp No." = c(Individual_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                              Individual_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                                row.names = individual_fluctuation_rnames_1)

individual_fluctuation_study_1 <- data.frame("Study" = c(Individual_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                         Individual_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                             row.names = individual_fluctuation_rnames_1)

Individual_Fluctuation_Model_CVR_Estimates_Reorder_1 <- Individual_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

individual_fluctuation_table_1 <- data.frame(estimate = Individual_Fluctuation_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                             lowerCL = Individual_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                             upperCL = Individual_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                             K = individual_fluctuation_k_1[,1], 
                                             group_no = individual_fluctuation_group_no_1[,1], 
                                             row.names = individual_fluctuation_rnames_1)
individual_fluctuation_table_1$name <- row.names(individual_fluctuation_table_1)

individual_fluctuation_raw_mean_1 <- c(unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                     select("InCVR"))), 
                                       unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                     select("InCVR"))))

individual_fluctuation_raw_name_1 <- c(replicate(481, "Sinusoidal (Sine Curve)"), 
                                       replicate(533, "Alternating"))

individual_fluctuation_raw_df_1 <- data.frame("Model" = individual_fluctuation_raw_name_1, 
                                              "Effect" = individual_fluctuation_raw_mean_1)

# Graph code - Part 1

Individual_Fluctuation_Order_1 <- c("Alternating", "Sinusoidal (Sine Curve)")

density_individual_fluctuation_CVR_1 <- individual_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Individual_Fluctuation_Order_1)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = individual_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Individual_Fluctuation_Order_1)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                        size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                              vjust = c(-2.7, -0.8))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                        scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                        coord_cartesian(xlim = c(-1, 1)) +
                                        annotate('text',  x = 1, y = (seq(1, dim(individual_fluctuation_table_1)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(individual_fluctuation_table_1["Alternating", "K"], 
                                                                      individual_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                    c(individual_fluctuation_table_1["Alternating", "group_no"], 
                                                                      individual_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                               paste(format(round(mean(exp(Individual_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.75, y = (seq(1, dim(individual_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_individual_fluctuation_CVR_1 #(400x240)

# Preparing Graph - Part 2

individual_fluctuation_rnames_2 <- c("Stepwise", "Stochastic")

individual_fluctuation_k_2 <- data.frame("k" = c(Individual_Fluctuation_Exploration["Stepwise", "Freq"], 
                                                 Individual_Fluctuation_Exploration["Stochastic", "Freq"]), 
                                         row.names = individual_fluctuation_rnames_2)

individual_fluctuation_group_no_2 <- data.frame("Spp No." = c(Individual_Fluctuation_Species_Count["Stepwise", "Freq"], 
                                                              Individual_Fluctuation_Species_Count["Stochastic", "Freq"]), 
                                                row.names = individual_fluctuation_rnames_2)

individual_fluctuation_study_2 <- data.frame("Study" = c(Individual_Fluctuation_Study_Count["Stepwise", "Freq"], 
                                                         Individual_Fluctuation_Study_Count["Stochastic", "Freq"]), 
                                             row.names = individual_fluctuation_rnames_2)

Individual_Fluctuation_Model_CVR_Estimates_Reorder_2 <- Individual_Fluctuation_Model_CVR_Estimates[c("Stepwise", "Stochastic"), ]

individual_fluctuation_table_2 <- data.frame(estimate = Individual_Fluctuation_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                             lowerCL = Individual_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                             upperCL = Individual_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                             K = individual_fluctuation_k_2[,1], 
                                             group_no = individual_fluctuation_group_no_2[,1], 
                                             row.names = individual_fluctuation_rnames_2)
individual_fluctuation_table_2$name <- row.names(individual_fluctuation_table_2)

individual_fluctuation_raw_mean_2 <- c(unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                     select("InCVR"))), 
                                       unlist(unname(Individual_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stochastic") %>% 
                                                     select("InCVR"))))

individual_fluctuation_raw_name_2 <- c(replicate(196, "Stepwise"), 
                                       replicate(20, "Stochastic"))

individual_fluctuation_raw_df_2 <- data.frame("Model" = individual_fluctuation_raw_name_2, 
                                              "Effect" = individual_fluctuation_raw_mean_2)

# Graph code - Part 2

Individual_Fluctuation_Order_2 <- c("Stochastic", "Stepwise")

density_individual_fluctuation_CVR_2 <- individual_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Individual_Fluctuation_Order_2)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = individual_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Individual_Fluctuation_Order_2)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(individual_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                        size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                              vjust = c(-2.7, -2.7))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                        scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                        coord_cartesian(xlim = c(-1, 1)) +
                                        annotate('text',  x = 1, y = (seq(1, dim(individual_fluctuation_table_2)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(individual_fluctuation_table_2["Stochastic", "K"], 
                                                                      individual_fluctuation_table_2["Stepwise", "K"]), "~","(", 
                                                                    c(individual_fluctuation_table_2["Stochastic", "group_no"], 
                                                                      individual_fluctuation_table_2["Stepwise", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Fluctuation_Model_CVR_Estimates["Stochastic", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                               paste(format(round(mean(exp(Individual_Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.75, y = (seq(1, dim(individual_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_individual_fluctuation_CVR_2 #(400x240)

##### Individual-Level Subset Model - Class Meta-Regression - CVR #####
Individual_Class_Exploration <- Individual_Subset_Data %>% select("Class") %>% table() %>% data.frame()
rownames(Individual_Class_Exploration) <- Individual_Class_Exploration$Class

Individual_Class_Data <- Individual_Subset_Data %>% filter(Class == "Actinopteri"|
                                                           Class == "Amphibia"|
                                                           Class == "Arachnida"|
                                                           Class == "Bivalvia"|
                                                           Class == "Branchiopoda"|
                                                           Class == "Gastropoda"|
                                                           Class == "Holothuroidea"|
                                                           Class == "Insecta"|
                                                           Class == "Malacostraca")

Individual_Class_Species_Count <- Individual_Class_Data %>% select("Scientific_Name", "Class") %>% table() %>% data.frame() %>%
                                  filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Individual_Class_Species_Count) <- Individual_Class_Species_Count$Class

Individual_Class_Study_Count <- Individual_Class_Data %>% select("Study_ID", "Class") %>% table() %>% data.frame() %>%
                                filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Individual_Class_Study_Count) <- Individual_Class_Study_Count$Class

Individual_Class_Species <- Individual_Class_Data %>% select("phylo") %>% unique()

Individual_Class_A_cor <- as.data.frame(A_cor)
Individual_Class_A_cor <- Individual_Class_A_cor[c(Individual_Class_Species$phylo), c(Individual_Class_Species$phylo)]
Individual_Class_A_cor <- as.matrix(Individual_Class_A_cor)

Individual_Class_VCV_InCVR <- make_VCV_matrix(Individual_Class_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  12ish minutes
  if(run){
    Individual_Class_Model_CVR <- metafor::rma.mv(InCVR, V = Individual_Class_VCV_InCVR, test = "t", dfs = "contain",
                                                  mods = ~ Class - 1,
                                                  random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                ~1|Shared_Animal_Number, ~1|Measurement), 
                                                  R = list(phylo=Individual_Class_A_cor), data = Individual_Class_Data, method = "REML", sparse = TRUE, 
                                                  control=list(rel.tol=1e-9))
    saveRDS(Individual_Class_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Individual_Class_Model_CVR.rds")
  } else {
            Individual_Class_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Individual_Class_Model_CVR.rds")})

Individual_Class_Model_CVR_rob <- robust(Individual_Class_Model_CVR, cluster = Individual_Class_Data$Study_ID, adjust = TRUE)

Individual_Class_Model_CVR_Estimates <- data.frame(Class = substr(row.names(Individual_Class_Model_CVR$b), 6, 100),
                                                   estimate = Individual_Class_Model_CVR$b, 
                                                   ci.lb = Individual_Class_Model_CVR$ci.lb, 
                                                   ci.ub = Individual_Class_Model_CVR$ci.ub)
rownames(Individual_Class_Model_CVR_Estimates) <- Individual_Class_Model_CVR_Estimates$Class
Individual_Class_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Individual_Class_Model_CVR), 2))

# Preparing Graph - Combined

individual_class_rnames <- c("Actinopteri", "Amphibia", "Arachnida", 
                             "Bivalvia", "Branchiopoda", "Gastropoda", "Holothuroidea", 
                             "Insecta", "Malacostraca")

individual_class_k <- data.frame("k" = c(Individual_Class_Exploration["Actinopteri", "Freq"], 
                                         Individual_Class_Exploration["Amphibia", "Freq"], 
                                         Individual_Class_Exploration["Arachnida", "Freq"], 
                                         Individual_Class_Exploration["Bivalvia", "Freq"], 
                                         Individual_Class_Exploration["Branchiopoda", "Freq"], 
                                         Individual_Class_Exploration["Gastropoda", "Freq"], 
                                         Individual_Class_Exploration["Holothuroidea", "Freq"],
                                         Individual_Class_Exploration["Insecta", "Freq"],
                                         Individual_Class_Exploration["Malacostraca", "Freq"]), 
                                 row.names = individual_class_rnames)

individual_class_group_no <- data.frame("Spp No." = c(Individual_Class_Species_Count["Actinopteri", "Freq"], 
                                                      Individual_Class_Species_Count["Amphibia", "Freq"], 
                                                      Individual_Class_Species_Count["Arachnida", "Freq"],
                                                      Individual_Class_Species_Count["Bivalvia", "Freq"],
                                                      Individual_Class_Species_Count["Branchiopoda", "Freq"],
                                                      Individual_Class_Species_Count["Gastropoda", "Freq"], 
                                                      Individual_Class_Species_Count["Holothuroidea", "Freq"],
                                                      Individual_Class_Species_Count["Insecta", "Freq"],
                                                      Individual_Class_Species_Count["Malacostraca", "Freq"]), 
                                        row.names = individual_class_rnames)

individual_class_study <- data.frame("Study" = c(Individual_Class_Study_Count["Actinopteri", "Freq"], 
                                                 Individual_Class_Study_Count["Amphibia", "Freq"], 
                                                 Individual_Class_Study_Count["Arachnida", "Freq"],
                                                 Individual_Class_Study_Count["Bivalvia", "Freq"],
                                                 Individual_Class_Study_Count["Branchiopoda", "Freq"],
                                                 Individual_Class_Study_Count["Gastropoda", "Freq"], 
                                                 Individual_Class_Study_Count["Holothuroidea", "Freq"],
                                                 Individual_Class_Study_Count["Insecta", "Freq"],
                                                 Individual_Class_Study_Count["Malacostraca", "Freq"]), 
                                        row.names = individual_class_rnames)

individual_class_table <- data.frame(estimate = Individual_Class_Model_CVR_Estimates[,"estimate"], 
                                     lowerCL = Individual_Class_Model_CVR_Estimates[,"ci.lb"], 
                                     upperCL = Individual_Class_Model_CVR_Estimates[,"ci.ub"], 
                                     K = individual_class_k[,1], 
                                     group_no = individual_class_group_no[,1], 
                                     row.names = individual_class_rnames)
individual_class_table$name <- row.names(individual_class_table)

individual_class_raw_mean <- c(unlist(unname(Individual_Class_Data %>% filter(`Class` == "Actinopteri") %>% 
                                               select("InCVR"))), 
                               unlist(unname(Individual_Class_Data %>% filter(`Class` == "Amphibia") %>% 
                                               select("InCVR"))), 
                               unlist(unname(Individual_Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                               select("InCVR"))), 
                               unlist(unname(Individual_Class_Data %>% filter(`Class` == "Bivalvia") %>% 
                                               select("InCVR"))),
                               unlist(unname(Individual_Class_Data %>% filter(`Class` == "Branchiopoda") %>% 
                                               select("InCVR"))),
                               unlist(unname(Individual_Class_Data %>% filter(`Class` == "Gastropoda") %>% 
                                               select("InCVR"))), 
                               unlist(unname(Individual_Class_Data %>% filter(`Class` == "Holothuroidea") %>% 
                                               select("InCVR"))),
                               unlist(unname(Individual_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                               select("InCVR"))),
                               unlist(unname(Individual_Class_Data %>% filter(`Class` == "Malacostraca") %>% 
                                               select("InCVR"))))

individual_class_raw_name <- c(replicate(131, "Actinopteri"), 
                               replicate(64, "Amphibia"), 
                               replicate(102, "Arachnida"), 
                               replicate(12, "Bivalvia"),
                               replicate(21, "Branchiopoda"),
                               replicate(21, "Gastropoda"), 
                               replicate(39, "Holothuroidea"),
                               replicate(614, "Insecta"),
                               replicate(56, "Malacostraca"))

individual_class_raw_df <- data.frame("Model" = individual_class_raw_name, 
                                      "Effect" = individual_class_raw_mean)

# Graph code - Combined

Individual_Class_Order <- c("Malacostraca", "Insecta", "Holothuroidea", "Gastropoda",  
                            "Branchiopoda", "Bivalvia", "Arachnida", "Amphibia", "Actinopteri")

density_individual_class_CVR <- individual_class_table %>% mutate(name = fct_relevel(name, Individual_Class_Order)) %>%
                                ggplot() +
                                geom_density_ridges(data = individual_class_raw_df %>% mutate(Model = fct_relevel(Model, Individual_Class_Order)), 
                                                    aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                    scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                geom_linerange(aes(y = rev(seq(1, dim(individual_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                               size = 1) +
                                geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                size = 1, fatten = 2) +
                                theme_bw() +
                                guides(fill = "none", colour = "none") +
                                labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                 vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7, 
                                                                           -2.7, -2.7, -2.7, -2.7))) +
                                theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                theme(axis.ticks = element_blank()) +
                                theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                scale_colour_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692", "#3D5E89", 
                                                               "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                scale_fill_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692", "#3D5E89", 
                                                             "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                coord_cartesian(xlim = c(-1, 1)) +
                                annotate('text',  x = 1, y = (seq(1, dim(individual_class_table)[1], 1)+0.4),
                                label= paste("italic(k)==", c(individual_class_table["Malacostraca", "K"], 
                                                              individual_class_table["Insecta", "K"], 
                                                              individual_class_table["Holothuroidea", "K"], 
                                                              individual_class_table["Gastropoda", "K"],
                                                              individual_class_table["Branchiopoda", "K"],
                                                              individual_class_table["Bivalvia", "K"],
                                                              individual_class_table["Arachnida", "K"],
                                                              individual_class_table["Amphibia", "K"],
                                                              individual_class_table["Actinopteri", "K"]), "~","(", 
                                                            c(individual_class_table["Malacostraca", "group_no"], 
                                                              individual_class_table["Insecta", "group_no"], 
                                                              individual_class_table["Holothuroidea", "group_no"], 
                                                              individual_class_table["Gastropoda", "group_no"],
                                                              individual_class_table["Branchiopoda", "group_no"],
                                                              individual_class_table["Bivalvia", "group_no"],
                                                              individual_class_table["Arachnida", "group_no"],
                                                              individual_class_table["Amphibia", "group_no"],
                                                              individual_class_table["Actinopteri", "group_no"]), 
                                             ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Malacostraca", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Holothuroidea", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Gastropoda", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                       paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Branchiopoda", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Bivalvia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Amphibia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                       paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Actinopteri", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                           x = -0.75, y = (seq(1, dim(individual_class_table)[1], 1)+0.4)), size = 3.5)

density_individual_class_CVR #(400x800)

# Preparing Graph - Part 1

individual_class_rnames_1 <- c("Actinopteri", "Amphibia", "Arachnida", "Bivalvia", "Branchiopoda")

individual_class_k_1 <- data.frame("k" = c(Individual_Class_Exploration["Actinopteri", "Freq"], 
                                           Individual_Class_Exploration["Amphibia", "Freq"], 
                                           Individual_Class_Exploration["Arachnida", "Freq"], 
                                           Individual_Class_Exploration["Bivalvia", "Freq"], 
                                           Individual_Class_Exploration["Branchiopoda", "Freq"]), 
                                   row.names = individual_class_rnames_1)

individual_class_group_no_1 <- data.frame("Spp No." = c(Individual_Class_Species_Count["Actinopteri", "Freq"], 
                                                        Individual_Class_Species_Count["Amphibia", "Freq"], 
                                                        Individual_Class_Species_Count["Arachnida", "Freq"],
                                                        Individual_Class_Species_Count["Bivalvia", "Freq"],
                                                        Individual_Class_Species_Count["Branchiopoda", "Freq"]), 
                                          row.names = individual_class_rnames_1)

individual_class_study_1 <- data.frame("Study" = c(Individual_Class_Study_Count["Actinopteri", "Freq"], 
                                                   Individual_Class_Study_Count["Amphibia", "Freq"], 
                                                   Individual_Class_Study_Count["Arachnida", "Freq"],
                                                   Individual_Class_Study_Count["Bivalvia", "Freq"],
                                                   Individual_Class_Study_Count["Branchiopoda", "Freq"]), 
                                       row.names = individual_class_rnames_1)

Individual_Class_Model_CVR_Estimates_Reorder_1 <- Individual_Class_Model_CVR_Estimates[c("Actinopteri", "Amphibia", "Arachnida", "Bivalvia", "Branchiopoda"), ]

individual_class_table_1 <- data.frame(estimate = Individual_Class_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                       lowerCL = Individual_Class_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                       upperCL = Individual_Class_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                       K = individual_class_k_1[,1], 
                                       group_no = individual_class_group_no_1[,1], 
                                       row.names = individual_class_rnames_1)
individual_class_table_1$name <- row.names(individual_class_table_1)

individual_class_raw_mean_1 <- c(unlist(unname(Individual_Class_Data %>% filter(`Class` == "Actinopteri") %>% 
                                               select("InCVR"))), 
                                 unlist(unname(Individual_Class_Data %>% filter(`Class` == "Amphibia") %>% 
                                               select("InCVR"))), 
                                 unlist(unname(Individual_Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                               select("InCVR"))), 
                                 unlist(unname(Individual_Class_Data %>% filter(`Class` == "Bivalvia") %>% 
                                               select("InCVR"))),
                                 unlist(unname(Individual_Class_Data %>% filter(`Class` == "Branchiopoda") %>% 
                                               select("InCVR"))))

individual_class_raw_name_1 <- c(replicate(131, "Actinopteri"), 
                                 replicate(64, "Amphibia"), 
                                 replicate(102, "Arachnida"), 
                                 replicate(12, "Bivalvia"),
                                 replicate(21, "Branchiopoda"))

individual_class_raw_df_1 <- data.frame("Model" = individual_class_raw_name_1, 
                                        "Effect" = individual_class_raw_mean_1)

# Graph code - Part 1

Individual_Class_Order_1 <- c("Branchiopoda", "Bivalvia", "Arachnida", "Amphibia", "Actinopteri")

density_individual_class_CVR_1 <- individual_class_table_1 %>% mutate(name = fct_relevel(name, Individual_Class_Order_1)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = individual_class_raw_df_1 %>% mutate(Model = fct_relevel(Model, Individual_Class_Order_1)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(individual_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                  size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                        vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#3D5E89", "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                  scale_fill_manual(values = c("#3D5E89", "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                  coord_cartesian(xlim = c(-1, 1)) +
                                  annotate('text',  x = 1, y = (seq(1, dim(individual_class_table_1)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(individual_class_table_1["Branchiopoda", "K"],
                                                                individual_class_table_1["Bivalvia", "K"],
                                                                individual_class_table_1["Arachnida", "K"],
                                                                individual_class_table_1["Amphibia", "K"],
                                                                individual_class_table_1["Actinopteri", "K"]), "~","(", 
                                                              c(individual_class_table_1["Branchiopoda", "group_no"],
                                                                individual_class_table_1["Bivalvia", "group_no"],
                                                                individual_class_table_1["Arachnida", "group_no"],
                                                                individual_class_table_1["Amphibia", "group_no"],
                                                                individual_class_table_1["Actinopteri", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Branchiopoda", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Bivalvia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Amphibia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Actinopteri", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.75, y = (seq(1, dim(individual_class_table_1)[1], 1)+0.4)), size = 3.5)

density_individual_class_CVR_1 #(400x480)

# Preparing Graph - Part 2

individual_class_rnames_2 <- c("Gastropoda", "Holothuroidea", "Insecta", "Malacostraca")

individual_class_k_2 <- data.frame("k" = c(Individual_Class_Exploration["Gastropoda", "Freq"], 
                                           Individual_Class_Exploration["Holothuroidea", "Freq"],
                                           Individual_Class_Exploration["Insecta", "Freq"],
                                           Individual_Class_Exploration["Malacostraca", "Freq"]), 
                                   row.names = individual_class_rnames_2)

individual_class_group_no_2 <- data.frame("Spp No." = c(Individual_Class_Species_Count["Gastropoda", "Freq"], 
                                                        Individual_Class_Species_Count["Holothuroidea", "Freq"],
                                                        Individual_Class_Species_Count["Insecta", "Freq"],
                                                        Individual_Class_Species_Count["Malacostraca", "Freq"]), 
                                          row.names = individual_class_rnames_2)

individual_class_study_2 <- data.frame("Study" = c(Individual_Class_Study_Count["Gastropoda", "Freq"], 
                                                   Individual_Class_Study_Count["Holothuroidea", "Freq"],
                                                   Individual_Class_Study_Count["Insecta", "Freq"],
                                                   Individual_Class_Study_Count["Malacostraca", "Freq"]), 
                                       row.names = individual_class_rnames_2)

Individual_Class_Model_CVR_Estimates_Reorder_2 <- Individual_Class_Model_CVR_Estimates[c("Gastropoda", "Holothuroidea", "Insecta", "Malacostraca"), ]

individual_class_table_2 <- data.frame(estimate = Individual_Class_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                       lowerCL = Individual_Class_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                       upperCL = Individual_Class_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                       K = individual_class_k_2[,1], 
                                       group_no = individual_class_group_no_2[,1], 
                                       row.names = individual_class_rnames_2)
individual_class_table_2$name <- row.names(individual_class_table_2)

individual_class_raw_mean_2 <- c(unlist(unname(Individual_Class_Data %>% filter(`Class` == "Gastropoda") %>% 
                                               select("InCVR"))), 
                                unlist(unname(Individual_Class_Data %>% filter(`Class` == "Holothuroidea") %>% 
                                               select("InCVR"))),
                                unlist(unname(Individual_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                               select("InCVR"))),
                                unlist(unname(Individual_Class_Data %>% filter(`Class` == "Malacostraca") %>% 
                                               select("InCVR"))))

individual_class_raw_name_2 <- c(replicate(21, "Gastropoda"), 
                                 replicate(39, "Holothuroidea"),
                                 replicate(614, "Insecta"),
                                 replicate(56, "Malacostraca"))

individual_class_raw_df_2 <- data.frame("Model" = individual_class_raw_name_2, 
                                        "Effect" = individual_class_raw_mean_2)

# Graph code - Part 2

Individual_Class_Order_2 <- c("Malacostraca", "Insecta", "Holothuroidea", "Gastropoda")

density_individual_class_CVR_2 <- individual_class_table_2 %>% mutate(name = fct_relevel(name, Individual_Class_Order_2)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = individual_class_raw_df_2 %>% mutate(Model = fct_relevel(Model, Individual_Class_Order_2)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(individual_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(individual_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                  size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                        vjust = c(-2.7, -2.7, -2.7, -2.7))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692")) +
                                  scale_fill_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692")) +
                                  coord_cartesian(xlim = c(-1, 1)) +
                                  annotate('text',  x = 1, y = (seq(1, dim(individual_class_table_2)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(individual_class_table_2["Malacostraca", "K"], 
                                                                individual_class_table_2["Insecta", "K"], 
                                                                individual_class_table_2["Holothuroidea", "K"], 
                                                                individual_class_table_2["Gastropoda", "K"]), "~","(", 
                                                              c(individual_class_table_2["Malacostraca", "group_no"], 
                                                                individual_class_table_2["Insecta", "group_no"], 
                                                                individual_class_table_2["Holothuroidea", "group_no"], 
                                                                individual_class_table_2["Gastropoda", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Malacostraca", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                         paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Holothuroidea", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                         paste(format(round(mean(exp(Individual_Class_Model_CVR_Estimates["Gastropoda", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.75, y = (seq(1, dim(individual_class_table_2)[1], 1)+0.4)), size = 3.5)

density_individual_class_CVR_2 #(400x400)

##### Summary of Individual Plots #####

Individual_Layout <- rbind(c(1, 2),
                           c(1, 2),
                           c(1, 2),
                           c(1, 2),
                           c(1, 3),
                           c(1, 3),
                           c(1, 3),
                           c(1, 3),
                           c(1, 3))

Individual_Combined_CVR <- grid.arrange(density_individual_class_CVR, density_individual_fluctuation_CVR, Individual_Amplitude_Plot_CVR, 
                                        layout_matrix = Individual_Layout)

Individual_Combined_CVR #(850 x 850)

##### Aquatic Subset Model - CVR #####
Aquatic_Subset_Data <- Individual_Subset_Data %>% filter(Ecosystem == "Aquatic")
Aquatic_Species <- Aquatic_Subset_Data %>% select("phylo") %>% unique()

Aquatic_A_cor <- as.data.frame(A_cor)
Aquatic_A_cor <- Aquatic_A_cor[c(Aquatic_Species$phylo), c(Aquatic_Species$phylo)]
Aquatic_A_cor <- as.matrix(Aquatic_A_cor)

Aquatic_VCV_InCVR <- make_VCV_matrix(Aquatic_Subset_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Model_CVR <- metafor::rma.mv(InCVR ~ 1, V = Aquatic_VCV_InCVR, test = "t", dfs = "contain",
                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                       ~1|Shared_Animal_Number, ~1|Measurement), 
                                         R = list(phylo=Aquatic_A_cor), data = Aquatic_Subset_Data, method = "REML", sparse = TRUE, 
                                         control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Aquatic_Model_CVR.rds")
  } else {
            Aquatic_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Aquatic_Model_CVR.rds")})

Aquatic_Model_CVR_rob <- robust(Aquatic_Model_CVR, cluster = Aquatic_Subset_Data$Study_ID, adjust = TRUE)

Aquatic_Model_CVR_Estimates <- data.frame(estimate = Aquatic_Model_CVR$b, 
                                          ci.lb = Aquatic_Model_CVR$ci.lb, 
                                          ci.ub = Aquatic_Model_CVR$ci.ub)
Aquatic_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Model_CVR), 2))

#### Aquatic Subset Model - Fluctuation Amplitude Meta-Regression - CVR ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Amplitude_Model_CVR <- metafor::rma.mv(InCVR, V = Aquatic_VCV_InCVR, test = "t", dfs = "contain",
                                                   mods = ~ T2_Magnitude - 1,
                                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                 ~1|Shared_Animal_Number, ~1|Measurement), 
                                                   R = list(phylo=Aquatic_A_cor), data = Aquatic_Subset_Data, method = "REML", sparse = TRUE, 
                                                   control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Amplitude_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Aquatic_Amplitude_Model_CVR.rds")
  } else {
            Aquatic_Amplitude_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Aquatic_Amplitude_Model_CVR.rds")})

Aquatic_Amplitude_Model_CVR_rob <- robust(Aquatic_Amplitude_Model_CVR, cluster = Aquatic_Subset_Data$Study_ID, adjust = TRUE)

Aquatic_Amplitude_Model_CVR_Estimates <- data.frame(estimate = Aquatic_Amplitude_Model_CVR$b, 
                                                    ci.lb = Aquatic_Amplitude_Model_CVR$ci.lb, 
                                                    ci.ub = Aquatic_Amplitude_Model_CVR$ci.ub)
Aquatic_Amplitude_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Amplitude_Model_CVR), 2))

# Graph Preparing

Aquatic_Plot_Data <- Aquatic_Subset_Data
Aquatic_Plot_Data <- Aquatic_Plot_Data %>% mutate(n_category = ifelse(n_R1.1 <= 25, "25", 
                                                               ifelse(n_R1.1 > 25 & n_R1.1 <= 50, "50", 
                                                               ifelse(n_R1.1 > 50 & n_R1.1 <= 75, "75", "> 75"))))

# Graph Code

Aquatic_Amplitude_Plot_CVR <- ggplot(Aquatic_Plot_Data, aes(x = T2_Magnitude, y = InCVR)) + 
                              geom_point(aes(x = T2_Magnitude, y = InCVR, 
                                         size = fct_relevel(n_category, c("25", "50", "75", "> 75"))), 
                                         shape = 21, fill = "#4292c6", alpha = 0.5) + 
                              labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (lnCVR)", 
                                   size = "Sample Size", title = "Aquatic Organisms") +
                              theme_bw() +
                              theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                              theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                              theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                              theme(legend.position = "bottom", legend.direction = "horizontal") + 
                              geom_hline(yintercept = Aquatic_Model_CVR_Estimates$estimate, lty = 2) + 
                              geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                              stat_poly_eq(formula = y ~ x, 
                              aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                              parse = TRUE) +
                              coord_cartesian(xlim = c(0, 30), 
                                              ylim = c(-2.5, 2.5))

Aquatic_Amplitude_Plot_CVR #(400x400)

#### Aquatic Subset Model - Type of Fluctuation Meta-Regression - CVR ####
Aquatic_Fluctuation_Data <- Aquatic_Subset_Data %>% filter(!is.na(Fluctuation_Category))

Aquatic_Fluctuation_Exploration <- Aquatic_Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Aquatic_Fluctuation_Exploration) <- Aquatic_Fluctuation_Exploration$Fluctuation_Category

Aquatic_Fluctuation_Data <- Aquatic_Fluctuation_Data %>% filter(Fluctuation_Category != "Stepwise" &
                                                                Fluctuation_Category != "Stochastic")

Aquatic_Fluctuation_Species_Count <- Aquatic_Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% table() %>% data.frame() %>% 
                                     filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Aquatic_Fluctuation_Species_Count) <- Aquatic_Fluctuation_Species_Count$Fluctuation_Category

Aquatic_Fluctuation_Study_Count <- Aquatic_Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% table() %>% data.frame() %>% 
                                   filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Aquatic_Fluctuation_Study_Count) <- Aquatic_Fluctuation_Study_Count$Fluctuation_Category

Aquatic_Fluctuation_Species <- Aquatic_Fluctuation_Data %>% select("phylo") %>% unique()

Aquatic_Fluctuation_A_cor <- as.data.frame(A_cor)
Aquatic_Fluctuation_A_cor <- Aquatic_Fluctuation_A_cor[c(Aquatic_Fluctuation_Species$phylo), c(Aquatic_Fluctuation_Species$phylo)]
Aquatic_Fluctuation_A_cor <- as.matrix(Aquatic_Fluctuation_A_cor)

Aquatic_Fluctuation_VCV_InCVR <- make_VCV_matrix(Aquatic_Fluctuation_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Fluctuation_Model_CVR <- metafor::rma.mv(InCVR, V = Aquatic_Fluctuation_VCV_InCVR, test = "t", dfs = "contain",
                                                     mods = ~ Fluctuation_Category - 1,
                                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                   ~1|Shared_Animal_Number, ~1|Measurement), 
                                                     R = list(phylo=Aquatic_Fluctuation_A_cor), data = Aquatic_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                     control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Fluctuation_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Aquatic_Fluctuation_Model_CVR.rds")
  } else {
            Aquatic_Fluctuation_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Aquatic_Fluctuation_Model_CVR.rds")})

Aquatic_Fluctuation_Model_CVR_rob <- robust(Aquatic_Fluctuation_Model_CVR, cluster = Aquatic_Fluctuation_Data$Study_ID, adjust = TRUE)

Aquatic_Fluctuation_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Aquatic_Fluctuation_Model_CVR$b), 21, 100),
                                                      estimate = Aquatic_Fluctuation_Model_CVR$b, 
                                                      ci.lb = Aquatic_Fluctuation_Model_CVR$ci.lb, 
                                                      ci.ub = Aquatic_Fluctuation_Model_CVR$ci.ub)
rownames(Aquatic_Fluctuation_Model_CVR_Estimates) <- Aquatic_Fluctuation_Model_CVR_Estimates$Category
Aquatic_Fluctuation_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Fluctuation_Model_CVR), 2))

# Preparing Graph - Combined

aquatic_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating")

aquatic_fluctuation_k <- data.frame("k" = c(Aquatic_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                            Aquatic_Fluctuation_Exploration["Alternating", "Freq"]), 
                                    row.names = aquatic_fluctuation_rnames)

aquatic_fluctuation_group_no <- data.frame("Spp No." = c(Aquatic_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                         Aquatic_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                           row.names = aquatic_fluctuation_rnames)

aquatic_fluctuation_study <- data.frame("Study" = c(Aquatic_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                    Aquatic_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                           row.names = aquatic_fluctuation_rnames)

Aquatic_Fluctuation_Model_CVR_Estimates_Reorder <- Aquatic_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

aquatic_fluctuation_table <- data.frame(estimate = Aquatic_Fluctuation_Model_CVR_Estimates_Reorder[,"estimate"], 
                                        lowerCL = Aquatic_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.lb"], 
                                        upperCL = Aquatic_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.ub"], 
                                        K = aquatic_fluctuation_k[,1], 
                                        group_no = aquatic_fluctuation_group_no[,1], 
                                        row.names = aquatic_fluctuation_rnames)
aquatic_fluctuation_table$name <- row.names(aquatic_fluctuation_table)

aquatic_fluctuation_raw_mean <- c( unlist(unname(Aquatic_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                   select("InCVR"))), 
                                   unlist(unname(Aquatic_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                  select("InCVR"))))

aquatic_fluctuation_raw_name <- c(replicate(267, "Sinusoidal (Sine Curve)"), 
                                  replicate(98, "Alternating"))

aquatic_fluctuation_raw_df <- data.frame("Model" = aquatic_fluctuation_raw_name, 
                                         "Effect" = aquatic_fluctuation_raw_mean)

# Graph code - Combined

Aquatic_Fluctuation_Order <- c("Alternating", "Sinusoidal (Sine Curve)")

density_aquatic_fluctuation_CVR <- aquatic_fluctuation_table %>% mutate(name = fct_relevel(name, Aquatic_Fluctuation_Order)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = aquatic_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Fluctuation_Order)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                    vjust = c(-2.7, -0.8))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                   scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(aquatic_fluctuation_table)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(aquatic_fluctuation_table["Alternating", "K"], 
                                                                 aquatic_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                               c(aquatic_fluctuation_table["Alternating", "group_no"], 
                                                                 aquatic_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                          paste(format(round(mean(exp(Aquatic_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(aquatic_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_fluctuation_CVR #(400x240)

# Preparing Graph - Part 1

aquatic_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)")

aquatic_fluctuation_k_1 <- data.frame("k" = c(Aquatic_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"]), 
                                      row.names = aquatic_fluctuation_rnames_1)

aquatic_fluctuation_group_no_1 <- data.frame("Spp No." = c(Aquatic_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"]), 
                                             row.names = aquatic_fluctuation_rnames_1)

aquatic_fluctuation_study_1 <- data.frame("Study" = c(Aquatic_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"]), 
                                          row.names = aquatic_fluctuation_rnames_1)

Aquatic_Fluctuation_Model_CVR_Estimates_Reorder_1 <- Aquatic_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)"), ]

aquatic_fluctuation_table_1 <- data.frame(estimate = Aquatic_Fluctuation_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                          lowerCL = Aquatic_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                          upperCL = Aquatic_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                          K = aquatic_fluctuation_k_1[,1], 
                                          group_no = aquatic_fluctuation_group_no_1[,1], 
                                          row.names = aquatic_fluctuation_rnames_1)
aquatic_fluctuation_table_1$name <- row.names(aquatic_fluctuation_table_1)

aquatic_fluctuation_raw_mean_1 <- c( unlist(unname(Aquatic_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                   select("InCVR"))))

aquatic_fluctuation_raw_name_1 <- c(replicate(267, "Sinusoidal (Sine Curve)"))

aquatic_fluctuation_raw_df_1 <- data.frame("Model" = aquatic_fluctuation_raw_name_1, 
                                           "Effect" = aquatic_fluctuation_raw_mean_1)

# Graph code - Part 1

Aquatic_Fluctuation_Order_1 <- c("Sinusoidal (Sine Curve)")

density_aquatic_fluctuation_CVR_1 <- aquatic_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Fluctuation_Order_1)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = aquatic_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Fluctuation_Order_1)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.25, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table_1)[1], 1)), xmin = min(aquatic_fluctuation_raw_df_1$Effect)-0.1, xmax = -1.5, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                     size = 1, fatten = 2) +
                                     theme_bw() +
                                     guides(fill = "none", colour = "none") +
                                     labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                           vjust = c(-0.8))) +
                                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                     theme(axis.ticks = element_blank()) +
                                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                     scale_colour_manual(values = c("#2B4E7A")) +
                                     scale_fill_manual(values = c("#2B4E7A")) +
                                     coord_cartesian(xlim = c(-1, 1)) +
                                     annotate('text',  x = 1, y = (seq(1, dim(aquatic_fluctuation_table_1)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(aquatic_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                 c(aquatic_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.75, y = (seq(1, dim(aquatic_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_aquatic_fluctuation_CVR_1 #(400x160)

# Preparing Graph - Part 2

aquatic_fluctuation_rnames_2 <- c("Alternating")

aquatic_fluctuation_k_2 <- data.frame("k" = c(Aquatic_Fluctuation_Exploration["Alternating", "Freq"]), 
                                      row.names = aquatic_fluctuation_rnames_2)

aquatic_fluctuation_group_no_2 <- data.frame("Spp No." = c(Aquatic_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                             row.names = aquatic_fluctuation_rnames_2)

aquatic_fluctuation_study_2 <- data.frame("Study" = c(Aquatic_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                          row.names = aquatic_fluctuation_rnames_2)

Aquatic_Fluctuation_Model_CVR_Estimates_Reorder_2 <- Aquatic_Fluctuation_Model_CVR_Estimates[c("Alternating"), ]

aquatic_fluctuation_table_2 <- data.frame(estimate = Aquatic_Fluctuation_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                          lowerCL = Aquatic_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                          upperCL = Aquatic_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                          K = aquatic_fluctuation_k_2[,1], 
                                          group_no = aquatic_fluctuation_group_no_2[,1], 
                                          row.names = aquatic_fluctuation_rnames_2)
aquatic_fluctuation_table_2$name <- row.names(aquatic_fluctuation_table_2)

aquatic_fluctuation_raw_mean_2 <- c(unlist(unname(Aquatic_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                   select("InCVR"))))

aquatic_fluctuation_raw_name_2 <- c(replicate(98, "Alternating"))

aquatic_fluctuation_raw_df_2 <- data.frame("Model" = aquatic_fluctuation_raw_name_2, 
                                           "Effect" = aquatic_fluctuation_raw_mean_2)

# Graph code - Part 2

Aquatic_Fluctuation_Order_2 <- c("Alternating")

density_aquatic_fluctuation_CVR_2 <- aquatic_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Fluctuation_Order_2)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = aquatic_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Fluctuation_Order_2)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.25, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(aquatic_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                     size = 1, fatten = 2) +
                                     theme_bw() +
                                     guides(fill = "none", colour = "none") +
                                     labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                           vjust = c(-2.7))) +
                                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                     theme(axis.ticks = element_blank()) +
                                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                     scale_colour_manual(values = c("#5D7AA1")) +
                                     scale_fill_manual(values = c("#5D7AA1")) +
                                     coord_cartesian(xlim = c(-1, 1)) +
                                     annotate('text',  x = 1, y = (seq(1, dim(aquatic_fluctuation_table_2)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(aquatic_fluctuation_table_2["Alternating", "K"]), "~","(", 
                                                                 c(aquatic_fluctuation_table_2["Alternating", "group_no"]), 
                                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.75, y = (seq(1, dim(aquatic_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_aquatic_fluctuation_CVR_2 #(400x160)

##### Aquatic Subset Model - Trait Meta-Regression - CVR #####
Aquatic_Trait_Exploration <- Aquatic_Subset_Data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Aquatic_Trait_Exploration) <- Aquatic_Trait_Exploration$Trait_Category

Aquatic_Trait_Species_Count <- Aquatic_Subset_Data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>%
                               filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Aquatic_Trait_Species_Count) <- Aquatic_Trait_Species_Count$Trait_Category

Aquatic_Trait_Study_Count <- Aquatic_Subset_Data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>%
                             filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Aquatic_Trait_Study_Count) <- Aquatic_Trait_Study_Count$Trait_Category

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Trait_Model_CVR <- metafor::rma.mv(InCVR, V = Aquatic_VCV_InCVR, test = "t", dfs = "contain",
                                               mods = ~ Trait_Category - 1,
                                               random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                             ~1|Shared_Animal_Number, ~1|Measurement), 
                                               R = list(phylo=Aquatic_A_cor), data = Aquatic_Subset_Data, method = "REML", sparse = TRUE, 
                                               control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Trait_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Aquatic_Trait_Model_CVR.rds")
  } else {
            Aquatic_Trait_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Aquatic_Trait_Model_CVR.rds")})

Aquatic_Trait_Model_CVR_rob <- robust(Aquatic_Trait_Model_CVR, cluster = Aquatic_Subset_Data$Study_ID, adjust = TRUE)

Aquatic_Trait_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Aquatic_Trait_Model_CVR$b), 15, 100),
                                                estimate = Aquatic_Trait_Model_CVR$b, 
                                                ci.lb = Aquatic_Trait_Model_CVR$ci.lb, 
                                                ci.ub = Aquatic_Trait_Model_CVR$ci.ub)
rownames(Aquatic_Trait_Model_CVR_Estimates) <- Aquatic_Trait_Model_CVR_Estimates$Category
Aquatic_Trait_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Trait_Model_CVR), 2))

# Preparing Graph - Combined

aquatic_trait_rnames <- c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits", "Morphology", "Physiological")

aquatic_trait_k <- data.frame("k" = c(Aquatic_Trait_Exploration["Behavioural", "Freq"], 
                                      Aquatic_Trait_Exploration["Biochemical Assay", "Freq"], 
                                      Aquatic_Trait_Exploration["Gene Expression", "Freq"], 
                                      Aquatic_Trait_Exploration["Life-History Traits", "Freq"], 
                                      Aquatic_Trait_Exploration["Morphology", "Freq"], 
                                      Aquatic_Trait_Exploration["Physiological", "Freq"]), 
                              row.names = aquatic_trait_rnames)

aquatic_trait_group_no <- data.frame("Spp No." = c(Aquatic_Trait_Species_Count["Behavioural", "Freq"], 
                                                   Aquatic_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                   Aquatic_Trait_Species_Count["Gene Expression", "Freq"], 
                                                   Aquatic_Trait_Species_Count["Life-History Traits", "Freq"],
                                                   Aquatic_Trait_Species_Count["Morphology", "Freq"],
                                                   Aquatic_Trait_Species_Count["Physiological", "Freq"]), 
                                     row.names = aquatic_trait_rnames)

aquatic_trait_study <- data.frame("Study" = c(Aquatic_Trait_Study_Count["Behavioural", "Freq"], 
                                              Aquatic_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                              Aquatic_Trait_Study_Count["Gene Expression", "Freq"], 
                                              Aquatic_Trait_Study_Count["Life-History Traits", "Freq"],
                                              Aquatic_Trait_Study_Count["Morphology", "Freq"],
                                              Aquatic_Trait_Study_Count["Physiological", "Freq"]), 
                                  row.names = aquatic_trait_rnames)

Aquatic_Trait_Model_CVR_Estimates_Reorder <- Aquatic_Trait_Model_CVR_Estimates[c("Behavioural", "Biochemical Assay", "Gene Expression", 
                                                                                 "Life-History Traits", "Morphology", "Physiological"), ]

aquatic_trait_table <- data.frame(estimate = Aquatic_Trait_Model_CVR_Estimates_Reorder[,"estimate"], 
                                  lowerCL = Aquatic_Trait_Model_CVR_Estimates_Reorder[,"ci.lb"], 
                                  upperCL = Aquatic_Trait_Model_CVR_Estimates_Reorder[,"ci.ub"], 
                                  K = aquatic_trait_k[,1], 
                                  group_no = aquatic_trait_group_no[,1], 
                                  row.names = aquatic_trait_rnames)
aquatic_trait_table$name <- row.names(aquatic_trait_table)

aquatic_trait_raw_mean <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Behavioural") %>% 
                                            select("InCVR"))), 
                            unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                            select("InCVR"))), 
                            unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Gene Expression") %>% 
                                            select("InCVR"))), 
                            unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                            select("InCVR"))), 
                            unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                            select("InCVR"))),
                            unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                            select("InCVR"))))

aquatic_trait_raw_name <- c(replicate(13, "Behavioural"), 
                            replicate(75, "Biochemical Assay"), 
                            replicate(38, "Gene Expression"), 
                            replicate(67, "Life-history Traits"), 
                            replicate(90, "Morphology"),
                            replicate(112, "Physiological"))

aquatic_trait_raw_df <- data.frame("Model" = aquatic_trait_raw_name, 
                                   "Effect" = aquatic_trait_raw_mean)

# Graph code - Combined

Aquatic_Trait_Order <- c("Physiological", "Morphology", "Life-history Traits", 
                         "Gene Expression", "Biochemical Assay", "Behavioural")

density_aquatic_trait_CVR <- aquatic_trait_table %>% mutate(name = fct_relevel(name, Aquatic_Trait_Order)) %>%
                             ggplot() +
                             geom_density_ridges(data = aquatic_trait_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Trait_Order)), 
                                                 aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                             geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                            size = 1) +
                             geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                             size = 1, fatten = 2) +
                             theme_bw() +
                             guides(fill = "none", colour = "none") +
                             labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                             theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                   vjust = c(-2.7, -2.7, -0.8, -0.8, -0.8, -2.7))) +
                             theme(axis.text.x = element_text(margin = margin(b = 5))) +
                             theme(axis.ticks = element_blank()) +
                             theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                             scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                             scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                             scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                             coord_cartesian(xlim = c(-1, 1)) +
                             annotate('text',  x = 1, y = (seq(1, dim(aquatic_trait_table)[1], 1)+0.4),
                             label= paste("italic(k)==", c(aquatic_trait_table["Physiological", "K"], 
                                                           aquatic_trait_table["Morphology", "K"], 
                                                           aquatic_trait_table["Life-history Traits", "K"], 
                                                           aquatic_trait_table["Gene Expression", "K"],
                                                           aquatic_trait_table["Biochemical Assay", "K"],
                                                           aquatic_trait_table["Behavioural", "K"]), "~","(", 
                                                         c(aquatic_trait_table["Physiological", "group_no"], 
                                                           aquatic_trait_table["Morphology", "group_no"], 
                                                           aquatic_trait_table["Life-history Traits", "group_no"], 
                                                           aquatic_trait_table["Gene Expression", "group_no"],
                                                           aquatic_trait_table["Biochemical Assay", "group_no"],
                                                           aquatic_trait_table["Behavioural", "group_no"]), 
                                          ")"), parse = TRUE, hjust = "right", size = 3.5) +
                             geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                    paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                    paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "% *"), 
                                                    paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Gene Expression", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                    paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                    paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Behavioural", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.75, y = (seq(1, dim(aquatic_trait_table)[1], 1)+0.4)), size = 3.5, 
                             fontface = c("plain", "plain", "bold", "plain", "plain", "plain"))

density_aquatic_trait_CVR #(400x540)

# Preparing Graph - Part 1

aquatic_trait_rnames_1 <- c("Behavioural", "Biochemical Assay", "Gene Expression")

aquatic_trait_k_1 <- data.frame("k" = c(Aquatic_Trait_Exploration["Behavioural", "Freq"], 
                                        Aquatic_Trait_Exploration["Biochemical Assay", "Freq"], 
                                        Aquatic_Trait_Exploration["Gene Expression", "Freq"]), 
                                row.names = aquatic_trait_rnames_1)

aquatic_trait_group_no_1 <- data.frame("Spp No." = c(Aquatic_Trait_Species_Count["Behavioural", "Freq"], 
                                                     Aquatic_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                     Aquatic_Trait_Species_Count["Gene Expression", "Freq"]), 
                                       row.names = aquatic_trait_rnames_1)

aquatic_trait_study_1 <- data.frame("Study" = c(Aquatic_Trait_Study_Count["Behavioural", "Freq"], 
                                                Aquatic_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                Aquatic_Trait_Study_Count["Gene Expression", "Freq"]), 
                                       row.names = aquatic_trait_rnames_1)

Aquatic_Trait_Model_CVR_Estimates_Reorder_1 <- Aquatic_Trait_Model_CVR_Estimates[c("Behavioural", "Biochemical Assay", "Gene Expression"), ]

aquatic_trait_table_1 <- data.frame(estimate = Aquatic_Trait_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                    lowerCL = Aquatic_Trait_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                    upperCL = Aquatic_Trait_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                    K = aquatic_trait_k_1[,1], 
                                    group_no = aquatic_trait_group_no_1[,1], 
                                    row.names = aquatic_trait_rnames_1)
aquatic_trait_table_1$name <- row.names(aquatic_trait_table_1)

aquatic_trait_raw_mean_1 <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Behavioural") %>% 
                                            select("InCVR"))), 
                              unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                            select("InCVR"))), 
                              unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Gene Expression") %>% 
                                            select("InCVR"))))

aquatic_trait_raw_name_1 <- c(replicate(13, "Behavioural"), 
                              replicate(75, "Biochemical Assay"), 
                              replicate(38, "Gene Expression"))

aquatic_trait_raw_df_1 <- data.frame("Model" = aquatic_trait_raw_name_1, 
                                     "Effect" = aquatic_trait_raw_mean_1)

# Graph code - Part 1

Aquatic_Trait_Order_1 <- c("Gene Expression", "Biochemical Assay", "Behavioural")

density_aquatic_trait_CVR_1 <- aquatic_trait_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Trait_Order_1)) %>%
                               ggplot() +
                               geom_density_ridges(data = aquatic_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Trait_Order_1)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                               size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                     vjust = c(-0.8, -0.8, -2.7))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                               scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                               coord_cartesian(xlim = c(-1, 1)) +
                               annotate('text',  x = 1, y = (seq(1, dim(aquatic_trait_table_1)[1], 1)+0.4),
                               label= paste("italic(k)==", c(aquatic_trait_table_1["Gene Expression", "K"],
                                                             aquatic_trait_table_1["Biochemical Assay", "K"],
                                                             aquatic_trait_table_1["Behavioural", "K"]), "~","(", 
                                                           c(aquatic_trait_table_1["Gene Expression", "group_no"],
                                                             aquatic_trait_table_1["Biochemical Assay", "group_no"],
                                                             aquatic_trait_table_1["Behavioural", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                               geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Gene Expression", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                      paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                      paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Behavioural", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                          x = -0.75, y = (seq(1, dim(aquatic_trait_table_1)[1], 1)+0.4)), size = 3.5, 
                               fontface = c("plain", "plain", "plain"))

density_aquatic_trait_CVR_1 #(400x320)

# Preparing Graph - Part 2

aquatic_trait_rnames_2 <- c("Life-history Traits", "Morphology", "Physiological")

aquatic_trait_k_2 <- data.frame("k" = c(Aquatic_Trait_Exploration["Life-History Traits", "Freq"], 
                                        Aquatic_Trait_Exploration["Morphology", "Freq"], 
                                        Aquatic_Trait_Exploration["Physiological", "Freq"]), 
                                row.names = aquatic_trait_rnames_2)

aquatic_trait_group_no_2 <- data.frame("Spp No." = c(Aquatic_Trait_Species_Count["Life-History Traits", "Freq"],
                                                     Aquatic_Trait_Species_Count["Morphology", "Freq"],
                                                     Aquatic_Trait_Species_Count["Physiological", "Freq"]), 
                                       row.names = aquatic_trait_rnames_2)

aquatic_trait_study_2 <- data.frame("Study" = c(Aquatic_Trait_Study_Count["Life-History Traits", "Freq"],
                                                Aquatic_Trait_Study_Count["Morphology", "Freq"],
                                                Aquatic_Trait_Study_Count["Physiological", "Freq"]), 
                                    row.names = aquatic_trait_rnames_2)

Aquatic_Trait_Model_CVR_Estimates_Reorder_2 <- Aquatic_Trait_Model_CVR_Estimates[c("Life-History Traits", "Morphology", "Physiological"), ]

aquatic_trait_table_2 <- data.frame(estimate = Aquatic_Trait_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                    lowerCL = Aquatic_Trait_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                    upperCL = Aquatic_Trait_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                    K = aquatic_trait_k_2[,1], 
                                    group_no = aquatic_trait_group_no_2[,1], 
                                    row.names = aquatic_trait_rnames_2)
aquatic_trait_table_2$name <- row.names(aquatic_trait_table_2)

aquatic_trait_raw_mean_2 <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                            select("InCVR"))), 
                              unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                            select("InCVR"))),
                              unlist(unname(Aquatic_Subset_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                            select("InCVR"))))

aquatic_trait_raw_name_2 <- c(replicate(67, "Life-history Traits"), 
                              replicate(90, "Morphology"),
                              replicate(112, "Physiological"))

aquatic_trait_raw_df_2 <- data.frame("Model" = aquatic_trait_raw_name_2, 
                                     "Effect" = aquatic_trait_raw_mean_2)

# Graph code - Part 2

Aquatic_Trait_Order_2 <- c("Physiological", "Morphology", "Life-history Traits")

density_aquatic_trait_CVR_2 <- aquatic_trait_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Trait_Order_2)) %>%
                               ggplot() +
                               geom_density_ridges(data = aquatic_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Trait_Order_2)), 
                                                   aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                   scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                               geom_linerange(aes(y = rev(seq(1, dim(aquatic_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                               geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                               theme_bw() +
                               guides(fill = "none", colour = "none") +
                               labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                               theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                     vjust = c(-2.7, -2.7, -0.8))) +
                               theme(axis.text.x = element_text(margin = margin(b = 5))) +
                               theme(axis.ticks = element_blank()) +
                               theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                               scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                               scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                               scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                               coord_cartesian(xlim = c(-1, 1)) +
                               annotate('text',  x = 1, y = (seq(1, dim(aquatic_trait_table_2)[1], 1)+0.4),
                               label= paste("italic(k)==", c(aquatic_trait_table_2["Physiological", "K"], 
                                                             aquatic_trait_table_2["Morphology", "K"], 
                                                             aquatic_trait_table_2["Life-history Traits", "K"]), "~","(", 
                                                           c(aquatic_trait_table_2["Physiological", "group_no"], 
                                                             aquatic_trait_table_2["Morphology", "group_no"], 
                                                             aquatic_trait_table_2["Life-history Traits", "group_no"]), 
                                            ")"), parse = TRUE, hjust = "right", size = 3.5) +
                               geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                      paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                      paste(format(round(mean(exp(Aquatic_Trait_Model_CVR_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "% *")), 
                                          x = -0.75, y = (seq(1, dim(aquatic_trait_table_2)[1], 1)+0.4)), size = 3.5, 
                               fontface = c("plain", "plain", "bold"))

density_aquatic_trait_CVR_2 #(400x320)

##### Aquatic Subset Model - Plasticity Mechanism Meta-Regression - CVR #####
Aquatic_Plasticity_Exploration <- Aquatic_Subset_Data %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Aquatic_Plasticity_Exploration) <- Aquatic_Plasticity_Exploration$Plasticity_Mechanism

Aquatic_Plasticity_Species_Count <- Aquatic_Subset_Data %>% select("Scientific_Name", "Plasticity_Mechanism") %>% table() %>% data.frame() %>%
                                    filter(`Freq` != 0) %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Aquatic_Plasticity_Species_Count) <- Aquatic_Plasticity_Species_Count$Plasticity_Mechanism

Aquatic_Plasticity_Study_Count <- Aquatic_Subset_Data %>% select("Study_ID", "Plasticity_Mechanism") %>% table() %>% data.frame() %>%
                                  filter(`Freq` != 0) %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Aquatic_Plasticity_Study_Count) <- Aquatic_Plasticity_Study_Count$Plasticity_Mechanism

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Plasticity_Model_CVR <- metafor::rma.mv(InCVR, V = Aquatic_VCV_InCVR, test = "t", dfs = "contain",
                                                    mods = ~ Plasticity_Mechanism - 1,
                                                    random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                  ~1|Shared_Animal_Number, ~1|Measurement), 
                                                    R = list(phylo=Aquatic_A_cor), data = Aquatic_Subset_Data, method = "REML", sparse = TRUE, 
                                                    control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Plasticity_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Aquatic_Plasticity_Model_CVR.rds")
  } else {
            Aquatic_Plasticity_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Aquatic_Plasticity_Model_CVR.rds")})

Aquatic_Plasticity_Model_CVR_rob <- robust(Aquatic_Plasticity_Model_CVR, cluster = Aquatic_Subset_Data$Study_ID, adjust = TRUE)

Aquatic_Plasticity_Model_CVR_Estimates <- data.frame(Plasticity_Mechanism = substr(row.names(Aquatic_Plasticity_Model_CVR$b), 21, 100),
                                                     estimate = Aquatic_Plasticity_Model_CVR$b, 
                                                     ci.lb = Aquatic_Plasticity_Model_CVR$ci.lb, 
                                                     ci.ub = Aquatic_Plasticity_Model_CVR$ci.ub)
rownames(Aquatic_Plasticity_Model_CVR_Estimates) <- Aquatic_Plasticity_Model_CVR_Estimates$Plasticity_Mechanism
Aquatic_Plasticity_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Plasticity_Model_CVR), 2))

# Preparing Graph - Combined

aquatic_plasticity_rnames <- c("Acclimation", "Development")

aquatic_plasticity_k <- data.frame("k" = c(Aquatic_Plasticity_Exploration["Acclimation", "Freq"], 
                                           Aquatic_Plasticity_Exploration["Developmental Plasticity", "Freq"]), 
                                   row.names = aquatic_plasticity_rnames)

aquatic_plasticity_group_no <- data.frame("Spp No." = c(Aquatic_Plasticity_Species_Count["Acclimation", "Freq"], 
                                                        Aquatic_Plasticity_Species_Count["Developmental Plasticity", "Freq"]), 
                                          row.names = aquatic_plasticity_rnames)

aquatic_plasticity_study <- data.frame("Study" = c(Aquatic_Plasticity_Study_Count["Acclimation", "Freq"], 
                                                   Aquatic_Plasticity_Study_Count["Developmental Plasticity", "Freq"]), 
                                          row.names = aquatic_plasticity_rnames)

aquatic_plasticity_table <- data.frame(estimate = Aquatic_Plasticity_Model_CVR_Estimates[,"estimate"], 
                                       lowerCL = Aquatic_Plasticity_Model_CVR_Estimates[,"ci.lb"], 
                                       upperCL = Aquatic_Plasticity_Model_CVR_Estimates[,"ci.ub"], 
                                       K = aquatic_plasticity_k[,1], 
                                       group_no = aquatic_plasticity_group_no[,1], 
                                       row.names = aquatic_plasticity_rnames)
aquatic_plasticity_table$name <- row.names(aquatic_plasticity_table)

aquatic_plasticity_raw_mean <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Plasticity_Mechanism` == "Acclimation") %>% 
                                                 select("InCVR"))), 
                                 unlist(unname(Aquatic_Subset_Data %>% filter(`Plasticity_Mechanism` == "Developmental Plasticity") %>% 
                                                 select("InCVR"))))

aquatic_plasticity_raw_name <- c(replicate(197, "Acclimation"), 
                                 replicate(198, "Development"))

aquatic_plasticity_raw_df <- data.frame("Model" = aquatic_plasticity_raw_name, 
                                        "Effect" = aquatic_plasticity_raw_mean)

# Graph code - Combined

Aquatic_Plasticity_Order <- c("Development", "Acclimation")

density_aquatic_plasticity_CVR <- aquatic_plasticity_table %>% mutate(name = fct_relevel(name, Aquatic_Plasticity_Order)) %>%
                                  ggplot() +
                                  geom_density_ridges(data = aquatic_plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Plasticity_Order)), 
                                                      aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                      scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                  geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                 size = 1) +
                                  geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                  size = 1, fatten = 2) +
                                  theme_bw() +
                                  guides(fill = "none", colour = "none") +
                                  labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                   vjust = c(-2.7, -2.7))) +
                                  theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                  theme(axis.ticks = element_blank()) +
                                  theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                  scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                  scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                  scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                  coord_cartesian(xlim = c(-1, 1)) +
                                  annotate('text',  x = 1, y = (seq(1, dim(aquatic_plasticity_table)[1], 1)+0.4),
                                  label= paste("italic(k)==", c(aquatic_plasticity_table["Development", "K"], 
                                                                aquatic_plasticity_table["Acclimation", "K"]), "~","(", 
                                                              c(aquatic_plasticity_table["Development", "group_no"], 
                                                                aquatic_plasticity_table["Acclimation", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                  geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Plasticity_Model_CVR_Estimates["Development", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                         paste(format(round(mean(exp(Aquatic_Plasticity_Model_CVR_Estimates["Acclimation", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                             x = -0.75, y = (seq(1, dim(aquatic_plasticity_table)[1], 1)+0.4)), size = 3.5)

density_aquatic_plasticity_CVR #(400x240)

# Preparing Graph - Part 1

aquatic_plasticity_rnames_1 <- c("Acclimation")

aquatic_plasticity_k_1 <- data.frame("k" = c(Aquatic_Plasticity_Exploration["Acclimation", "Freq"]), 
                                     row.names = aquatic_plasticity_rnames_1)

aquatic_plasticity_group_no_1 <- data.frame("Spp No." = c(Aquatic_Plasticity_Species_Count["Acclimation", "Freq"]), 
                                            row.names = aquatic_plasticity_rnames_1)

aquatic_plasticity_study_1 <- data.frame("Study" = c(Aquatic_Plasticity_Study_Count["Acclimation", "Freq"]), 
                                         row.names = aquatic_plasticity_rnames_1)

Aquatic_Plasticity_Model_CVR_Estimates_Reorder_1 <- Aquatic_Plasticity_Model_CVR_Estimates[c("Acclimation"), ]

aquatic_plasticity_table_1 <- data.frame(estimate = Aquatic_Plasticity_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                         lowerCL = Aquatic_Plasticity_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                         upperCL = Aquatic_Plasticity_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                         K = aquatic_plasticity_k_1[,1], 
                                         group_no = aquatic_plasticity_group_no_1[,1], 
                                         row.names = aquatic_plasticity_rnames_1)
aquatic_plasticity_table_1$name <- row.names(aquatic_plasticity_table_1)

aquatic_plasticity_raw_mean_1 <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Plasticity_Mechanism` == "Acclimation") %>% 
                                                 select("InCVR"))))

aquatic_plasticity_raw_name_1 <- c(replicate(197, "Acclimation"))

aquatic_plasticity_raw_df_1 <- data.frame("Model" = aquatic_plasticity_raw_name_1, 
                                          "Effect" = aquatic_plasticity_raw_mean_1)

# Graph code - Part 1

Aquatic_Plasticity_Order_1 <- c("Acclimation")

density_aquatic_plasticity_CVR_1 <- aquatic_plasticity_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Plasticity_Order_1)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = aquatic_plasticity_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Plasticity_Order_1)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                        scale = 0.25, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table_1)[1], 1)), xmin = min(aquatic_plasticity_raw_df_1$Effect)-0.1, xmax = -1.5, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                    size = 1, fatten = 2) +
                                    theme_bw() +
                                    guides(fill = "none", colour = "none") +
                                    labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                    theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                          vjust = c(-2.7))) +
                                    theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                    theme(axis.ticks = element_blank()) +
                                    theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                    scale_colour_manual(values = c("#2B4E7A")) +
                                    scale_fill_manual(values = c("#2B4E7A")) +
                                    coord_cartesian(xlim = c(-1, 1)) +
                                    annotate('text',  x = 1, y = (seq(1, dim(aquatic_plasticity_table_1)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(aquatic_plasticity_table_1["Acclimation", "K"]), "~","(", 
                                                                c(aquatic_plasticity_table_1["Acclimation", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                    geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Plasticity_Model_CVR_Estimates["Acclimation", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                               x = -0.75, y = (seq(1, dim(aquatic_plasticity_table_1)[1], 1)+0.4)), size = 3.5)

density_aquatic_plasticity_CVR_1 #(400x160)

# Preparing Graph - Part 2

aquatic_plasticity_rnames_2 <- c("Development")

aquatic_plasticity_k_2 <- data.frame("k" = c(Aquatic_Plasticity_Exploration["Developmental Plasticity", "Freq"]), 
                                     row.names = aquatic_plasticity_rnames_2)

aquatic_plasticity_group_no_2 <- data.frame("Spp No." = c(Aquatic_Plasticity_Species_Count["Developmental Plasticity", "Freq"]), 
                                            row.names = aquatic_plasticity_rnames_2)

aquatic_plasticity_study_2 <- data.frame("Study" = c(Aquatic_Plasticity_Study_Count["Developmental Plasticity", "Freq"]), 
                                         row.names = aquatic_plasticity_rnames_2)

Aquatic_Plasticity_Model_CVR_Estimates_Reorder_2 <- Aquatic_Plasticity_Model_CVR_Estimates[c("Development"), ]

aquatic_plasticity_table_2 <- data.frame(estimate = Aquatic_Plasticity_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                         lowerCL = Aquatic_Plasticity_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                         upperCL = Aquatic_Plasticity_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                         K = aquatic_plasticity_k_2[,1], 
                                         group_no = aquatic_plasticity_group_no_2[,1], 
                                         row.names = aquatic_plasticity_rnames_2)
aquatic_plasticity_table_2$name <- row.names(aquatic_plasticity_table_2)

aquatic_plasticity_raw_mean_2 <- c(unlist(unname(Aquatic_Subset_Data %>% filter(`Plasticity_Mechanism` == "Developmental Plasticity") %>% 
                                                 select("InCVR"))))

aquatic_plasticity_raw_name_2 <- c(replicate(198, "Development"))

aquatic_plasticity_raw_df_2 <- data.frame("Model" = aquatic_plasticity_raw_name_2, 
                                          "Effect" = aquatic_plasticity_raw_mean_2)

# Graph code - Part 2

Aquatic_Plasticity_Order_2 <- c("Development")

density_aquatic_plasticity_CVR_2 <- aquatic_plasticity_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Plasticity_Order_2)) %>%
                                    ggplot() +
                                    geom_density_ridges(data = aquatic_plasticity_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Plasticity_Order_2)), 
                                                        aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                        scale = 0.27, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                    geom_linerange(aes(y = rev(seq(1, dim(aquatic_plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                   size = 1) +
                                    geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                    size = 1, fatten = 2) +
                                    theme_bw() +
                                    guides(fill = "none", colour = "none") +
                                    labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                    theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                          vjust = c(-2.7))) +
                                    theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                    theme(axis.ticks = element_blank()) +
                                    theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                    scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                    scale_colour_manual(values = c("#5D7AA1")) +
                                    scale_fill_manual(values = c("#5D7AA1")) +
                                    coord_cartesian(xlim = c(-1, 1)) +
                                    annotate('text',  x = 1, y = (seq(1, dim(aquatic_plasticity_table_2)[1], 1)+0.4),
                                    label= paste("italic(k)==", c(aquatic_plasticity_table_2["Development", "K"]), "~","(", 
                                                                c(aquatic_plasticity_table_2["Development", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                    geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Plasticity_Model_CVR_Estimates["Development", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                               x = -0.75, y = (seq(1, dim(aquatic_plasticity_table_2)[1], 1)+0.4)), size = 3.5)

density_aquatic_plasticity_CVR_2 #(400x160)

##### Aquatic Subset Model - Specific Trait Meta-Regression - CVR #####
Aquatic_Specific_Trait_Exploration <- Aquatic_Subset_Data %>% select("Measurement") %>% table() %>% data.frame()
Aquatic_Specific_Trait_Exploration <- Aquatic_Specific_Trait_Exploration %>% filter(Freq > 10)
rownames(Aquatic_Specific_Trait_Exploration) <- Aquatic_Specific_Trait_Exploration$Measurement

Aquatic_Specific_Trait_Data <- Aquatic_Subset_Data %>% filter(Measurement == "Apparent Digestability Coefficient"| 
                                                              Measurement == "Cortisol"|
                                                              Measurement == "Development Time"|
                                                              Measurement == "Food Consumption"|
                                                              Measurement == "Immune Defense"|
                                                              Measurement == "Length"|
                                                              Measurement == "Locomotor Performance"|
                                                              Measurement == "Mass"|
                                                              Measurement == "Metabolic Rate")

Aquatic_Specific_Trait_Species_Count <- Aquatic_Specific_Trait_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>%
                                        filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Aquatic_Specific_Trait_Species_Count) <- Aquatic_Specific_Trait_Species_Count$Measurement

Aquatic_Specific_Trait_Study_Count <- Aquatic_Specific_Trait_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>%
                                      filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Aquatic_Specific_Trait_Study_Count) <- Aquatic_Specific_Trait_Study_Count$Measurement

Aquatic_Specific_Trait_Species <- Aquatic_Specific_Trait_Data %>% select("phylo") %>% unique()

Aquatic_Specific_Trait_Fluctuation_A_cor <- as.data.frame(A_cor)
Aquatic_Specific_Trait_Fluctuation_A_cor <- Aquatic_Specific_Trait_Fluctuation_A_cor[c(Aquatic_Specific_Trait_Species$phylo), c(Aquatic_Specific_Trait_Species$phylo)]
Aquatic_Specific_Trait_Fluctuation_A_cor <- as.matrix(Aquatic_Specific_Trait_Fluctuation_A_cor)

Aquatic_Specific_Trait_Fluctuation_VCV_InCVR <- make_VCV_matrix(Aquatic_Specific_Trait_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Aquatic_Specific_Trait_Model_CVR <- metafor::rma.mv(InCVR, V = Aquatic_Specific_Trait_Fluctuation_VCV_InCVR, test = "t", dfs = "contain",
                                                        mods = ~ Measurement - 1,
                                                        random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                      ~1|Shared_Animal_Number), 
                                                        R = list(phylo=Aquatic_Specific_Trait_Fluctuation_A_cor), data = Aquatic_Specific_Trait_Data, method = "REML", sparse = TRUE, 
                                                        control=list(rel.tol=1e-9))
    saveRDS(Aquatic_Specific_Trait_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Aquatic_Specific_Trait_Model_CVR.rds")
  } else {
            Aquatic_Specific_Trait_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Aquatic_Specific_Trait_Model_CVR.rds")})

Aquatic_Specific_Trait_Model_CVR_rob <- robust(Aquatic_Specific_Trait_Model_CVR, cluster = Aquatic_Specific_Trait_Data$Study_ID, adjust = TRUE)

Aquatic_Specific_Trait_Model_CVR_Estimates <- data.frame(Trait = substr(row.names(Aquatic_Specific_Trait_Model_CVR$b), 12, 100),
                                                         estimate = Aquatic_Specific_Trait_Model_CVR$b, 
                                                         ci.lb = Aquatic_Specific_Trait_Model_CVR$ci.lb, 
                                                         ci.ub = Aquatic_Specific_Trait_Model_CVR$ci.ub)
rownames(Aquatic_Specific_Trait_Model_CVR_Estimates) <- Aquatic_Specific_Trait_Model_CVR_Estimates$Trait
Aquatic_Specific_Trait_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Aquatic_Specific_Trait_Model_CVR), 2))

# Preparing Graph - Combined

aquatic_specific_trait_rnames <- c("Apparent Digestibility Coefficient","Cortisol Levels", "Development Time", 
                                   "Food Consumption", "Immune Defense", "Length", "Locomotor Performance", "Mass", "Metabolic Rate")

aquatic_specific_trait_k <- data.frame("k" = c(Aquatic_Specific_Trait_Exploration["Apparent Digestability Coefficient", "Freq"], 
                                               Aquatic_Specific_Trait_Exploration["Cortisol", "Freq"], 
                                               Aquatic_Specific_Trait_Exploration["Development Time", "Freq"], 
                                               Aquatic_Specific_Trait_Exploration["Food Consumption", "Freq"],
                                               Aquatic_Specific_Trait_Exploration["Immune Defense", "Freq"], 
                                               Aquatic_Specific_Trait_Exploration["Length", "Freq"], 
                                               Aquatic_Specific_Trait_Exploration["Locomotor Performance", "Freq"], 
                                               Aquatic_Specific_Trait_Exploration["Mass", "Freq"], 
                                               Aquatic_Specific_Trait_Exploration["Metabolic Rate", "Freq"]), 
                                       row.names = aquatic_specific_trait_rnames)

aquatic_specific_trait_group_no <- data.frame("Spp No." = c(Aquatic_Specific_Trait_Species_Count["Apparent Digestability Coefficient", "Freq"],
                                                            Aquatic_Specific_Trait_Species_Count["Cortisol", "Freq"], 
                                                            Aquatic_Specific_Trait_Species_Count["Development Time", "Freq"],
                                                            Aquatic_Specific_Trait_Species_Count["Food Consumption", "Freq"],
                                                            Aquatic_Specific_Trait_Species_Count["Immune Defense", "Freq"], 
                                                            Aquatic_Specific_Trait_Species_Count["Length", "Freq"], 
                                                            Aquatic_Specific_Trait_Species_Count["Locomotor Performance", "Freq"],  
                                                            Aquatic_Specific_Trait_Species_Count["Mass", "Freq"], 
                                                            Aquatic_Specific_Trait_Species_Count["Metabolic Rate", "Freq"]), 
                                              row.names = aquatic_specific_trait_rnames)

aquatic_specific_trait_study <- data.frame("Study" = c(Aquatic_Specific_Trait_Study_Count["Apparent Digestability Coefficient", "Freq"],
                                                       Aquatic_Specific_Trait_Study_Count["Cortisol", "Freq"], 
                                                       Aquatic_Specific_Trait_Study_Count["Development Time", "Freq"],
                                                       Aquatic_Specific_Trait_Study_Count["Food Consumption", "Freq"],
                                                       Aquatic_Specific_Trait_Study_Count["Immune Defense", "Freq"], 
                                                       Aquatic_Specific_Trait_Study_Count["Length", "Freq"], 
                                                       Aquatic_Specific_Trait_Study_Count["Locomotor Performance", "Freq"],  
                                                       Aquatic_Specific_Trait_Study_Count["Mass", "Freq"], 
                                                       Aquatic_Specific_Trait_Study_Count["Metabolic Rate", "Freq"]), 
                                           row.names = aquatic_specific_trait_rnames)

Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder <- Aquatic_Specific_Trait_Model_CVR_Estimates[c("Apparent Digestability Coefficient","Cortisol", "Development Time", 
                                                                                                   "Food Consumption", "Immune Defense", "Length", "Locomotor Performance", 
                                                                                                   "Mass", "Metabolic Rate"), ]

aquatic_specific_trait_table <- data.frame(estimate = Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder[,"estimate"], 
                                           lowerCL = Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder[,"ci.lb"], 
                                           upperCL = Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder[,"ci.ub"], 
                                           K = aquatic_specific_trait_k[,1], 
                                           group_no = aquatic_specific_trait_group_no[,1], 
                                           row.names = aquatic_specific_trait_rnames)
aquatic_specific_trait_table$name <- row.names(aquatic_specific_trait_table)

aquatic_specific_trait_raw_mean <- c(unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Apparent Digestability Coefficient") %>% 
                                                     select("InCVR"))),
                                     unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Cortisol") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                                     select("InCVR"))),
                                     unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Food Consumption") %>% 
                                                     select("InCVR"))),
                                     unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Immune Defense") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Locomotor Performance") %>% 
                                                     select("InCVR"))),
                                     unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Metabolic Rate") %>% 
                                                     select("InCVR"))))

aquatic_specific_trait_raw_name <- c(replicate(16, "Apparent Digestibility Coefficient"), 
                                     replicate(13, "Cortisol Levels"), 
                                     replicate(35, "Development Time"),
                                     replicate(13, "Food Consumption"),
                                     replicate(12, "Immune Defense"), 
                                     replicate(43, "Length"), 
                                     replicate(11, "Locomotor Performance"),  
                                     replicate(20, "Mass"), 
                                     replicate(24, "Metabolic Rate"))

aquatic_specific_trait_raw_df <- data.frame("Model" = aquatic_specific_trait_raw_name, 
                                            "Effect" = aquatic_specific_trait_raw_mean)

# Graph code - Combined

Aquatic_Specific_Trait_Order <- c("Metabolic Rate", "Mass", "Locomotor Performance", "Length", 
                                  "Immune Defense", "Food Consumption", "Development Time", "Cortisol Levels", "Apparent Digestibility Coefficient")

density_aquatic_specific_trait_CVR <- aquatic_specific_trait_table %>% mutate(name = fct_relevel(name, Aquatic_Specific_Trait_Order)) %>%
                                      ggplot() +
                                      geom_density_ridges(data = aquatic_specific_trait_raw_df %>% mutate(Model = fct_relevel(Model, Aquatic_Specific_Trait_Order)), 
                                                          aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                              scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                      geom_linerange(aes(y = rev(seq(1, dim(aquatic_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                     size = 1) +
                                      geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                      size = 1, fatten = 2) +
                                      theme_bw() +
                                      guides(fill = "none", colour = "none") +
                                      labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                      theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                            vjust = c(-0.8, -2.7, -0.8, -2.7, -0.8, 
                                                      -0.8, -0.8, -0.8, -0.4))) +
                                      theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                      theme(axis.ticks = element_blank()) +
                                      theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                      scale_colour_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692", "#3D5E89", 
                                                                     "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                      scale_fill_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692", "#3D5E89", 
                                                                   "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                      coord_cartesian(xlim = c(-1, 1)) +
                                      annotate('text',  x = 1, y = (seq(1, dim(aquatic_specific_trait_table)[1], 1)+0.4),
                                      label= paste("italic(k)==", c(aquatic_specific_trait_table["Metabolic Rate", "K"],
                                                                    aquatic_specific_trait_table["Mass", "K"], 
                                                                    aquatic_specific_trait_table["Locomotor Performance", "K"], 
                                                                    aquatic_specific_trait_table["Length", "K"], 
                                                                    aquatic_specific_trait_table["Immune Defense", "K"],
                                                                    aquatic_specific_trait_table["Food Consumption", "K"],
                                                                    aquatic_specific_trait_table["Development Time", "K"],
                                                                    aquatic_specific_trait_table["Cortisol Levels", "K"],
                                                                    aquatic_specific_trait_table["Apparent Digestibility Coefficient", "K"]), "~","(", 
                                                                  c(aquatic_specific_trait_table["Metabolic Rate", "group_no"],
                                                                    aquatic_specific_trait_table["Mass", "group_no"],
                                                                    aquatic_specific_trait_table["Locomotor Performance", "group_no"], 
                                                                    aquatic_specific_trait_table["Length", "group_no"], 
                                                                    aquatic_specific_trait_table["Immune Defense", "group_no"],
                                                                    aquatic_specific_trait_table["Food Consumption", "group_no"], 
                                                                    aquatic_specific_trait_table["Development Time", "group_no"],
                                                                    aquatic_specific_trait_table["Cortisol Levels", "group_no"], 
                                                                    aquatic_specific_trait_table["Apparent Digestibility Coefficient", "group_no"]), 
                                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                       geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                              paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                              paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Locomotor Performance", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                              paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                              paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Immune Defense", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                              paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Food Consumption", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                              paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "% *"),
                                                              paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Cortisol", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                              paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Apparent Digestability Coefficient", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                      x = -0.75, y = (seq(1, dim(aquatic_specific_trait_table)[1], 1)+0.4)), size = 3.5, 
                                       fontface = c("plain", "plain", "plain", "plain", "plain", 
                                                    "plain", "bold", "plain", "plain"))

density_aquatic_specific_trait_CVR #(400x800)

# Preparing Graph - Part 1

aquatic_specific_trait_rnames_1 <- c("Apparent Digestibility Coefficient","Cortisol Levels", "Development Time", 
                                     "Food Consumption", "Immune Defense")

aquatic_specific_trait_k_1 <- data.frame("k" = c(Aquatic_Specific_Trait_Exploration["Apparent Digestability Coefficient", "Freq"], 
                                                 Aquatic_Specific_Trait_Exploration["Cortisol", "Freq"], 
                                                 Aquatic_Specific_Trait_Exploration["Development Time", "Freq"], 
                                                 Aquatic_Specific_Trait_Exploration["Food Consumption", "Freq"],
                                                 Aquatic_Specific_Trait_Exploration["Immune Defense", "Freq"]), 
                                         row.names = aquatic_specific_trait_rnames_1)

aquatic_specific_trait_group_no_1 <- data.frame("Spp No." = c(Aquatic_Specific_Trait_Species_Count["Apparent Digestability Coefficient", "Freq"],
                                                              Aquatic_Specific_Trait_Species_Count["Cortisol", "Freq"], 
                                                              Aquatic_Specific_Trait_Species_Count["Development Time", "Freq"],
                                                              Aquatic_Specific_Trait_Species_Count["Food Consumption", "Freq"],
                                                              Aquatic_Specific_Trait_Species_Count["Immune Defense", "Freq"]), 
                                                row.names = aquatic_specific_trait_rnames_1)

aquatic_specific_trait_study_1 <- data.frame("Study" = c(Aquatic_Specific_Trait_Study_Count["Apparent Digestability Coefficient", "Freq"],
                                                         Aquatic_Specific_Trait_Study_Count["Cortisol", "Freq"], 
                                                         Aquatic_Specific_Trait_Study_Count["Development Time", "Freq"],
                                                         Aquatic_Specific_Trait_Study_Count["Food Consumption", "Freq"],
                                                         Aquatic_Specific_Trait_Study_Count["Immune Defense", "Freq"]), 
                                                row.names = aquatic_specific_trait_rnames_1)

Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder_1 <- Aquatic_Specific_Trait_Model_CVR_Estimates[c("Apparent Digestability Coefficient","Cortisol", "Development Time", 
                                                                                                     "Food Consumption", "Immune Defense"), ]

aquatic_specific_trait_table_1 <- data.frame(estimate = Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                             lowerCL = Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                             upperCL = Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                             K = aquatic_specific_trait_k_1[,1], 
                                             group_no = aquatic_specific_trait_group_no_1[,1], 
                                             row.names = aquatic_specific_trait_rnames_1)
aquatic_specific_trait_table_1$name <- row.names(aquatic_specific_trait_table_1)

aquatic_specific_trait_raw_mean_1 <- c(unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Apparent Digestability Coefficient") %>% 
                                                     select("InCVR"))),
                                       unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Cortisol") %>% 
                                                     select("InCVR"))), 
                                       unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                                     select("InCVR"))),
                                       unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Food Consumption") %>% 
                                                     select("InCVR"))),
                                       unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Immune Defense") %>% 
                                                     select("InCVR"))))

aquatic_specific_trait_raw_name_1 <- c(replicate(16, "Apparent Digestibility Coefficient"), 
                                       replicate(13, "Cortisol Levels"), 
                                       replicate(35, "Development Time"),
                                       replicate(13, "Food Consumption"),
                                       replicate(12, "Immune Defense"))

aquatic_specific_trait_raw_df_1 <- data.frame("Model" = aquatic_specific_trait_raw_name_1, 
                                              "Effect" = aquatic_specific_trait_raw_mean_1)

# Graph code - Part 1

Aquatic_Specific_Trait_Order_1 <- c("Immune Defense", "Food Consumption", "Development Time", "Cortisol Levels", "Apparent Digestibility Coefficient")

density_aquatic_specific_trait_CVR_1 <- aquatic_specific_trait_table_1 %>% mutate(name = fct_relevel(name, Aquatic_Specific_Trait_Order_1)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = aquatic_specific_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Aquatic_Specific_Trait_Order_1)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(aquatic_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                        size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                              vjust = c(-0.8, -0.8, -0.8, -0.8, -0.4))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#3D5E89", "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                        scale_fill_manual(values = c("#3D5E89", "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                        coord_cartesian(xlim = c(-1, 1)) +
                                        annotate('text',  x = 1, y = (seq(1, dim(aquatic_specific_trait_table_1)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(aquatic_specific_trait_table_1["Immune Defense", "K"],
                                                                      aquatic_specific_trait_table_1["Food Consumption", "K"],
                                                                      aquatic_specific_trait_table_1["Development Time", "K"],
                                                                      aquatic_specific_trait_table_1["Cortisol Levels", "K"],
                                                                      aquatic_specific_trait_table_1["Apparent Digestibility Coefficient", "K"]), "~","(", 
                                                                    c(aquatic_specific_trait_table_1["Immune Defense", "group_no"],
                                                                      aquatic_specific_trait_table_1["Food Consumption", "group_no"], 
                                                                      aquatic_specific_trait_table_1["Development Time", "group_no"],
                                                                      aquatic_specific_trait_table_1["Cortisol Levels", "group_no"], 
                                                                      aquatic_specific_trait_table_1["Apparent Digestibility Coefficient", "group_no"]), 
                                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Immune Defense", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                               paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Food Consumption", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                               paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "% *"),
                                                               paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Cortisol", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                               paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Apparent Digestability Coefficient", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                  x = -0.75, y = (seq(1, dim(aquatic_specific_trait_table_1)[1], 1)+0.4)), size = 3.5, 
                                                  fontface = c("plain", "plain", "bold", "plain", "plain"))

density_aquatic_specific_trait_CVR_1 #(400x480)

# Preparing Graph - Part 2

aquatic_specific_trait_rnames_2 <- c("Length", "Locomotor Performance", "Mass", "Metabolic Rate")

aquatic_specific_trait_k_2 <- data.frame("k" = c(Aquatic_Specific_Trait_Exploration["Length", "Freq"], 
                                                 Aquatic_Specific_Trait_Exploration["Locomotor Performance", "Freq"], 
                                                 Aquatic_Specific_Trait_Exploration["Mass", "Freq"], 
                                                 Aquatic_Specific_Trait_Exploration["Metabolic Rate", "Freq"]), 
                                         row.names = aquatic_specific_trait_rnames_2)

aquatic_specific_trait_group_no_2 <- data.frame("Spp No." = c(Aquatic_Specific_Trait_Species_Count["Length", "Freq"], 
                                                              Aquatic_Specific_Trait_Species_Count["Locomotor Performance", "Freq"],  
                                                              Aquatic_Specific_Trait_Species_Count["Mass", "Freq"], 
                                                              Aquatic_Specific_Trait_Species_Count["Metabolic Rate", "Freq"]), 
                                                row.names = aquatic_specific_trait_rnames_2)

aquatic_specific_trait_study_2 <- data.frame("Study" = c(Aquatic_Specific_Trait_Study_Count["Length", "Freq"], 
                                                         Aquatic_Specific_Trait_Study_Count["Locomotor Performance", "Freq"],  
                                                         Aquatic_Specific_Trait_Study_Count["Mass", "Freq"], 
                                                         Aquatic_Specific_Trait_Study_Count["Metabolic Rate", "Freq"]), 
                                             row.names = aquatic_specific_trait_rnames_2)

Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder_2 <- Aquatic_Specific_Trait_Model_CVR_Estimates[c("Length", "Locomotor Performance", "Mass", "Metabolic Rate"), ]

aquatic_specific_trait_table_2 <- data.frame(estimate = Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                             lowerCL = Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                             upperCL = Aquatic_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                             K = aquatic_specific_trait_k_2[,1], 
                                             group_no = aquatic_specific_trait_group_no_2[,1], 
                                             row.names = aquatic_specific_trait_rnames_2)
aquatic_specific_trait_table_2$name <- row.names(aquatic_specific_trait_table_2)

aquatic_specific_trait_raw_mean_2 <- c(unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                                     select("InCVR"))), 
                                       unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Locomotor Performance") %>% 
                                                     select("InCVR"))),
                                       unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                                     select("InCVR"))), 
                                       unlist(unname(Aquatic_Specific_Trait_Data %>% filter(`Measurement` == "Metabolic Rate") %>% 
                                                     select("InCVR"))))

aquatic_specific_trait_raw_name_2 <- c(replicate(43, "Length"), 
                                       replicate(11, "Locomotor Performance"),  
                                       replicate(20, "Mass"), 
                                       replicate(24, "Metabolic Rate"))

aquatic_specific_trait_raw_df_2 <- data.frame("Model" = aquatic_specific_trait_raw_name_2, 
                                              "Effect" = aquatic_specific_trait_raw_mean_2)

# Graph code - Part 2

Aquatic_Specific_Trait_Order_2 <- c("Metabolic Rate", "Mass", "Locomotor Performance", "Length")

density_aquatic_specific_trait_CVR_2 <- aquatic_specific_trait_table_2 %>% mutate(name = fct_relevel(name, Aquatic_Specific_Trait_Order_2)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = aquatic_specific_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Aquatic_Specific_Trait_Order_2)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(aquatic_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_linerange(aes(y = rev(seq(1, dim(aquatic_specific_trait_table_2)[1], 1)), xmin = min(aquatic_specific_trait_raw_df_2$Effect)-0.1, xmax = -1.5, colour = name),
                                                           size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(aquatic_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                            size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                              vjust = c(-0.8, -2.7, -0.8, -2.7))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692")) +
                                        scale_fill_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692")) +
                                        coord_cartesian(xlim = c(-1, 1)) +
                                        annotate('text',  x = 1, y = (seq(1, dim(aquatic_specific_trait_table_2)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(aquatic_specific_trait_table_2["Metabolic Rate", "K"],
                                                                      aquatic_specific_trait_table_2["Mass", "K"], 
                                                                      aquatic_specific_trait_table_2["Locomotor Performance", "K"], 
                                                                      aquatic_specific_trait_table_2["Length", "K"]), "~","(", 
                                                                    c(aquatic_specific_trait_table_2["Metabolic Rate", "group_no"],
                                                                      aquatic_specific_trait_table_2["Mass", "group_no"],
                                                                      aquatic_specific_trait_table_2["Locomotor Performance", "group_no"], 
                                                                      aquatic_specific_trait_table_2["Length", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                               paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                               paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Locomotor Performance", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                               paste(format(round(mean(exp(Aquatic_Specific_Trait_Model_CVR_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.75, y = (seq(1, dim(aquatic_specific_trait_table_2)[1], 1)+0.4)), size = 3.5, 
                                        fontface = c("plain", "plain", "plain", "plain"))

density_aquatic_specific_trait_CVR_2 #(400x400)


##### Summary of Aquatic Plots #####

Aquatic_Layout <- rbind(c(1, 2), 
                        c(1, 2),
                        c(1, 2),
                        c(1, 2),
                        c(1, 2),
                        c(1, 2),
                        c(1, 3),
                        c(1, 3),
                        c(1, 3),
                        c(1, 4),
                        c(1, 4),
                        c(1, 4))

Aquatic_Combined_CVR <- grid.arrange(density_aquatic_specific_trait_CVR, density_aquatic_trait_CVR, 
                                     density_aquatic_fluctuation_CVR, density_aquatic_plasticity_CVR, 
                                     layout_matrix = Aquatic_Layout)

Aquatic_Combined_CVR #(850 x 900 - does not include amplitude plot)

##### Terrestrial Subset Model - CVR #####
Terrestrial_Subset_Data <- Individual_Subset_Data %>% filter(Ecosystem == "Terrestrial")
Terrestrial_Species <- Terrestrial_Subset_Data %>% select("phylo") %>% unique()

Terrestrial_A_cor <- as.data.frame(A_cor)
Terrestrial_A_cor <- Terrestrial_A_cor[c(Terrestrial_Species$phylo), c(Terrestrial_Species$phylo)]
Terrestrial_A_cor <- as.matrix(Terrestrial_A_cor)

Terrestrial_VCV_InCVR <- make_VCV_matrix(Terrestrial_Subset_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  8ish minutes
  if(run){
    Terrestrial_Model_CVR <- metafor::rma.mv(InCVR ~ 1, V = Terrestrial_VCV_InCVR, test = "t", dfs = "contain",
                                             random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                           ~1|Shared_Animal_Number, ~1|Measurement), 
                                             R = list(phylo=Terrestrial_A_cor), data = Terrestrial_Subset_Data, method = "REML", sparse = TRUE, 
                                             control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Terrestrial_Model_CVR.rds")
  } else {
            Terrestrial_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Terrestrial_Model_CVR.rds")})

Terrestrial_Model_CVR_rob <- robust(Terrestrial_Model_CVR, cluster = Terrestrial_Subset_Data$Study_ID, adjust = TRUE)

Terrestrial_Model_CVR_Estimates <- data.frame(estimate = Terrestrial_Model_CVR$b, 
                                              ci.lb = Terrestrial_Model_CVR$ci.lb, 
                                              ci.ub = Terrestrial_Model_CVR$ci.ub)
Terrestrial_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Model_CVR), 2))

#### Terrestrial Subset Model - Fluctuation Amplitude Meta-Regression - CVR ####
run <- FALSE
system.time( #  3ish minutes
  if(run){
    Terrestrial_Amplitude_Model_CVR <- metafor::rma.mv(InCVR, V = Terrestrial_VCV_InCVR, test = "t", dfs = "contain",
                                                       mods = ~ T2_Magnitude - 1,
                                                       random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                     ~1|Shared_Animal_Number, ~1|Measurement), 
                                                       R = list(phylo=Terrestrial_A_cor), data = Terrestrial_Subset_Data, method = "REML", sparse = TRUE, 
                                                       control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Amplitude_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Terrestrial_Amplitude_Model_CVR.rds")
  } else {
            Terrestrial_Amplitude_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Terrestrial_Amplitude_Model_CVR.rds")})

Terrestrial_Amplitude_Model_CVR_rob <- robust(Terrestrial_Amplitude_Model_CVR, cluster = Terrestrial_Subset_Data$Study_ID, adjust = TRUE)

Terrestrial_Amplitude_Model_CVR_Estimates <- data.frame(estimate = Terrestrial_Amplitude_Model_CVR$b, 
                                                        ci.lb = Terrestrial_Amplitude_Model_CVR$ci.lb, 
                                                        ci.ub = Terrestrial_Amplitude_Model_CVR$ci.ub)
Terrestrial_Amplitude_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Amplitude_Model_CVR), 2))

# Graph Preparing

Terrestrial_Plot_Data <- Terrestrial_Subset_Data
Terrestrial_Plot_Data <- Terrestrial_Plot_Data %>% mutate(n_category = ifelse(n_R1.1 <= 25, "25", 
                                                                       ifelse(n_R1.1 > 25 & n_R1.1 <= 50, "50", 
                                                                       ifelse(n_R1.1 > 50 & n_R1.1 <= 75, "75", "> 75"))))

# Graph Code

Terrestrial_Amplitude_Plot_CVR <- ggplot(Terrestrial_Plot_Data, aes(x = T2_Magnitude, y = InCVR)) + 
                                  geom_point(aes(x = T2_Magnitude, y = InCVR, 
                                             size = fct_relevel(n_category, c("25", "50", "75", "> 75"))), 
                                             shape = 21, fill = "#4292c6", alpha = 0.5) + 
                                  labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (lnCVR)", 
                                       size = "Sample Size", title = "Terrestrial Organisms") +
                                  theme_bw() +
                                  theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                                  theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                                  theme(legend.position = "bottom", legend.direction = "horizontal") + 
                                  geom_hline(yintercept = Terrestrial_Model_CVR_Estimates$estimate, lty = 2) + 
                                  geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                                  stat_poly_eq(formula = y ~ x, 
                                               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                               parse = TRUE) +
                                  coord_cartesian(xlim = c(0, 30), 
                                                  ylim = c(-2.5, 2.5))

Terrestrial_Amplitude_Plot_CVR #(400x400)

#### Terrestrial Subset Model - Type of Fluctuation Meta-Regression - CVR ####
Terrestrial_Fluctuation_Data <- Terrestrial_Subset_Data %>% filter(!is.na(Fluctuation_Category))

Terrestrial_Fluctuation_Exploration <- Terrestrial_Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Fluctuation_Exploration) <- Terrestrial_Fluctuation_Exploration$Fluctuation_Category

Terrestrial_Fluctuation_Species_Count <- Terrestrial_Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                         filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Fluctuation_Species_Count) <- Terrestrial_Fluctuation_Species_Count$Fluctuation_Category

Terrestrial_Fluctuation_Study_Count <- Terrestrial_Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                       filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Fluctuation_Study_Count) <- Terrestrial_Fluctuation_Study_Count$Fluctuation_Category

Terrestrial_Fluctuation_Species <- Terrestrial_Fluctuation_Data %>% select("phylo") %>% unique()

Terrestrial_Fluctuation_A_cor <- as.data.frame(A_cor)
Terrestrial_Fluctuation_A_cor <- Terrestrial_Fluctuation_A_cor[c(Terrestrial_Fluctuation_Species$phylo), c(Terrestrial_Fluctuation_Species$phylo)]
Terrestrial_Fluctuation_A_cor <- as.matrix(Terrestrial_Fluctuation_A_cor)

Terrestrial_Fluctuation_VCV_InCVR <- make_VCV_matrix(Terrestrial_Fluctuation_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  6ish minutes
  if(run){
    Terrestrial_Fluctuation_Model_CVR <- metafor::rma.mv(InCVR, V = Terrestrial_Fluctuation_VCV_InCVR, test = "t", dfs = "contain",
                                                         mods = ~ Fluctuation_Category - 1,
                                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                       ~1|Shared_Animal_Number, ~1|Measurement), 
                                                         R = list(phylo=Terrestrial_Fluctuation_A_cor), data = Terrestrial_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                         control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Fluctuation_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Terrestrial_Fluctuation_Model_CVR.rds")
  } else {
            Terrestrial_Fluctuation_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Terrestrial_Fluctuation_Model_CVR.rds")})

Terrestrial_Fluctuation_Model_CVR_rob <- robust(Terrestrial_Fluctuation_Model_CVR, cluster = Terrestrial_Fluctuation_Data$Study_ID, adjust = TRUE)

Terrestrial_Fluctuation_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Terrestrial_Fluctuation_Model_CVR$b), 21, 100),
                                                          estimate = Terrestrial_Fluctuation_Model_CVR$b, 
                                                          ci.lb = Terrestrial_Fluctuation_Model_CVR$ci.lb, 
                                                          ci.ub = Terrestrial_Fluctuation_Model_CVR$ci.ub)
rownames(Terrestrial_Fluctuation_Model_CVR_Estimates) <- Terrestrial_Fluctuation_Model_CVR_Estimates$Category
Terrestrial_Fluctuation_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Fluctuation_Model_CVR), 2))

# Preparing Graph - Combined

terrestrial_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise", "Stochastic")

terrestrial_fluctuation_k <- data.frame("k" = c(Terrestrial_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"],
                                                Terrestrial_Fluctuation_Exploration["Alternating", "Freq"], 
                                                Terrestrial_Fluctuation_Exploration["Stepwise", "Freq"], 
                                                Terrestrial_Fluctuation_Exploration["Stochastic", "Freq"]), 
                                        row.names = terrestrial_fluctuation_rnames)

terrestrial_fluctuation_group_no <- data.frame("Spp No." = c(Terrestrial_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                             Terrestrial_Fluctuation_Species_Count["Alternating", "Freq"], 
                                                             Terrestrial_Fluctuation_Species_Count["Stepwise", "Freq"], 
                                                             Terrestrial_Fluctuation_Species_Count["Stochastic", "Freq"]), 
                                               row.names = terrestrial_fluctuation_rnames)

terrestrial_fluctuation_study <- data.frame("Study" = c(Terrestrial_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                        Terrestrial_Fluctuation_Study_Count["Alternating", "Freq"], 
                                                        Terrestrial_Fluctuation_Study_Count["Stepwise", "Freq"], 
                                                        Terrestrial_Fluctuation_Study_Count["Stochastic", "Freq"]), 
                                               row.names = terrestrial_fluctuation_rnames)

Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder <- Terrestrial_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise", "Stochastic"), ]

terrestrial_fluctuation_table <- data.frame(estimate = Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder[,"estimate"], 
                                            lowerCL = Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.lb"], 
                                            upperCL = Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.ub"], 
                                            K = terrestrial_fluctuation_k[,1], 
                                            group_no = terrestrial_fluctuation_group_no[,1], 
                                            row.names = terrestrial_fluctuation_rnames)
terrestrial_fluctuation_table$name <- row.names(terrestrial_fluctuation_table)

terrestrial_fluctuation_raw_mean <- c(unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                      select("InCVR"))), 
                                      unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                      select("InCVR"))), 
                                      unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                      select("InCVR"))), 
                                      unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stochastic") %>% 
                                                      select("InCVR"))))

terrestrial_fluctuation_raw_name <- c(replicate(214, "Sinusoidal (Sine Curve)"), 
                                      replicate(435, "Alternating"), 
                                      replicate(186, "Stepwise"), 
                                      replicate(16, "Stochastic"))

terrestrial_fluctuation_raw_df <- data.frame("Model" = terrestrial_fluctuation_raw_name, 
                                             "Effect" = terrestrial_fluctuation_raw_mean)

# Graph code - Combined

Terrestrial_Fluctuation_Order <- c("Stochastic", "Stepwise", 
                                   "Alternating", "Sinusoidal (Sine Curve)")

density_terrestrial_fluctuation_CVR <- terrestrial_fluctuation_table %>% mutate(name = fct_relevel(name, Terrestrial_Fluctuation_Order)) %>%
                                       ggplot() +
                                       geom_density_ridges(data = terrestrial_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Fluctuation_Order)), 
                                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                       geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                      size = 1) +
                                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                       size = 1, fatten = 2) +
                                       theme_bw() +
                                       guides(fill = "none", colour = "none") +
                                       labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                        vjust = c(-2.7, -2.7, -2.7, -0.8))) +
                                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                       theme(axis.ticks = element_blank()) +
                                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                       scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                       scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                       coord_cartesian(xlim = c(-1, 1)) +
                                       annotate('text',  x = 1, y = (seq(1, dim(terrestrial_fluctuation_table)[1], 1)+0.4),
                                       label= paste("italic(k)==", c(terrestrial_fluctuation_table["Stochastic", "K"], 
                                                                     terrestrial_fluctuation_table["Stepwise", "K"], 
                                                                     terrestrial_fluctuation_table["Alternating", "K"], 
                                                                     terrestrial_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                   c(terrestrial_fluctuation_table["Stochastic", "group_no"], 
                                                                     terrestrial_fluctuation_table["Stepwise", "group_no"], 
                                                                     terrestrial_fluctuation_table["Alternating", "group_no"], 
                                                                     terrestrial_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                       geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_CVR_Estimates["Stochastic", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                              paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                              paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                              paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                  x = -0.75, y = (seq(1, dim(terrestrial_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_fluctuation_CVR #(400x400)

# Preparing Graph - Part 1

terrestrial_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)", "Alternating")

terrestrial_fluctuation_k_1 <- data.frame("k" = c(Terrestrial_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"],
                                                  Terrestrial_Fluctuation_Exploration["Alternating", "Freq"]), 
                                          row.names = terrestrial_fluctuation_rnames_1)

terrestrial_fluctuation_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                               Terrestrial_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                                 row.names = terrestrial_fluctuation_rnames_1)

terrestrial_fluctuation_study_1 <- data.frame("Study" = c(Terrestrial_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                          Terrestrial_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                              row.names = terrestrial_fluctuation_rnames_1)

Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder_1 <- Terrestrial_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

terrestrial_fluctuation_table_1 <- data.frame(estimate = Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                              lowerCL = Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                              upperCL = Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                              K = terrestrial_fluctuation_k_1[,1], 
                                              group_no = terrestrial_fluctuation_group_no_1[,1], 
                                              row.names = terrestrial_fluctuation_rnames_1)
terrestrial_fluctuation_table_1$name <- row.names(terrestrial_fluctuation_table_1)

terrestrial_fluctuation_raw_mean_1 <- c(unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                      select("InCVR"))), 
                                        unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                      select("InCVR"))))

terrestrial_fluctuation_raw_name_1 <- c(replicate(214, "Sinusoidal (Sine Curve)"), 
                                        replicate(435, "Alternating"))

terrestrial_fluctuation_raw_df_1 <- data.frame("Model" = terrestrial_fluctuation_raw_name_1, 
                                               "Effect" = terrestrial_fluctuation_raw_mean_1)

# Graph code - Part 1

Terrestrial_Fluctuation_Order_1 <- c("Alternating", "Sinusoidal (Sine Curve)")

density_terrestrial_fluctuation_CVR_1 <- terrestrial_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Fluctuation_Order_1)) %>%
                                         ggplot() +
                                         geom_density_ridges(data = terrestrial_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Fluctuation_Order_1)), 
                                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                 scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                         geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                        size = 1) +
                                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                         size = 1, fatten = 2) +
                                         theme_bw() +
                                         guides(fill = "none", colour = "none") +
                                         labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                         theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                               vjust = c(-2.7, -0.8))) +
                                         theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                         theme(axis.ticks = element_blank()) +
                                         theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                         scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                         scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                         coord_cartesian(xlim = c(-1, 1)) +
                                         annotate('text',  x = 1, y = (seq(1, dim(terrestrial_fluctuation_table_1)[1], 1)+0.4),
                                         label= paste("italic(k)==", c(terrestrial_fluctuation_table_1["Alternating", "K"], 
                                                                       terrestrial_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                     c(terrestrial_fluctuation_table_1["Alternating", "group_no"], 
                                                                       terrestrial_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                         geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                    x = -0.75, y = (seq(1, dim(terrestrial_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_fluctuation_CVR_1 #(400x240)

# Preparing Graph - Part 2

terrestrial_fluctuation_rnames_2 <- c("Stepwise", "Stochastic")

terrestrial_fluctuation_k_2 <- data.frame("k" = c(Terrestrial_Fluctuation_Exploration["Stepwise", "Freq"], 
                                                  Terrestrial_Fluctuation_Exploration["Stochastic", "Freq"]), 
                                          row.names = terrestrial_fluctuation_rnames_2)

terrestrial_fluctuation_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Fluctuation_Species_Count["Stepwise", "Freq"], 
                                                               Terrestrial_Fluctuation_Species_Count["Stochastic", "Freq"]), 
                                                 row.names = terrestrial_fluctuation_rnames_2)

terrestrial_fluctuation_study_2 <- data.frame("Study" = c(Terrestrial_Fluctuation_Study_Count["Stepwise", "Freq"], 
                                                          Terrestrial_Fluctuation_Study_Count["Stochastic", "Freq"]), 
                                              row.names = terrestrial_fluctuation_rnames_2)

Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder_2 <- Terrestrial_Fluctuation_Model_CVR_Estimates[c("Stepwise", "Stochastic"), ]

terrestrial_fluctuation_table_2 <- data.frame(estimate = Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                              lowerCL = Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                              upperCL = Terrestrial_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                              K = terrestrial_fluctuation_k_2[,1], 
                                              group_no = terrestrial_fluctuation_group_no_2[,1], 
                                              row.names = terrestrial_fluctuation_rnames_2)
terrestrial_fluctuation_table_2$name <- row.names(terrestrial_fluctuation_table_2)

terrestrial_fluctuation_raw_mean_2 <- c(unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                      select("InCVR"))), 
                                        unlist(unname(Terrestrial_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stochastic") %>% 
                                                      select("InCVR"))))

terrestrial_fluctuation_raw_name_2 <- c(replicate(186, "Stepwise"), 
                                        replicate(16, "Stochastic"))

terrestrial_fluctuation_raw_df_2 <- data.frame("Model" = terrestrial_fluctuation_raw_name_2, 
                                               "Effect" = terrestrial_fluctuation_raw_mean_2)

# Graph code - Part 2

Terrestrial_Fluctuation_Order_2 <- c("Stochastic", "Stepwise")

density_terrestrial_fluctuation_CVR_2 <- terrestrial_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Fluctuation_Order_2)) %>%
                                         ggplot() +
                                         geom_density_ridges(data = terrestrial_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Fluctuation_Order_2)), 
                                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                 scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                         geom_linerange(aes(y = rev(seq(1, dim(terrestrial_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                        size = 1) +
                                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                         size = 1, fatten = 2) +
                                         theme_bw() +
                                         guides(fill = "none", colour = "none") +
                                         labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                         theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                               vjust = c(-2.7, -2.7))) +
                                         theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                         theme(axis.ticks = element_blank()) +
                                         theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                         scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                         scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                         coord_cartesian(xlim = c(-1, 1)) +
                                         annotate('text',  x = 1, y = (seq(1, dim(terrestrial_fluctuation_table_2)[1], 1)+0.4),
                                         label= paste("italic(k)==", c(terrestrial_fluctuation_table["Stochastic", "K"], 
                                                                       terrestrial_fluctuation_table["Stepwise", "K"]), "~","(", 
                                                                     c(terrestrial_fluctuation_table["Stochastic", "group_no"], 
                                                                       terrestrial_fluctuation_table["Stepwise", "group_no"]), 
                                                      ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                         geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_CVR_Estimates["Stochastic", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                paste(format(round(mean(exp(Terrestrial_Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                    x = -0.75, y = (seq(1, dim(terrestrial_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_fluctuation_CVR_2 #(400x240)

##### Terrestrial Subset Model - Trait Meta-Regression - CVR #####
Terrestrial_Trait_Exploration <- Terrestrial_Subset_Data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Trait_Exploration) <- Terrestrial_Trait_Exploration$Trait_Category

Terrestrial_Trait_Species_Count <- Terrestrial_Subset_Data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>%
                                   filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Trait_Species_Count) <- Terrestrial_Trait_Species_Count$Trait_Category

Terrestrial_Trait_Study_Count <- Terrestrial_Subset_Data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>%
                                 filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Terrestrial_Trait_Study_Count) <- Terrestrial_Trait_Study_Count$Trait_Category

run <- FALSE
system.time( #  8ish minutes
  if(run){
    Terrestrial_Trait_Model_CVR <- metafor::rma.mv(InCVR, V = Terrestrial_VCV_InCVR, test = "t", dfs = "contain",
                                                   mods = ~ Trait_Category - 1,
                                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                 ~1|Shared_Animal_Number, ~1|Measurement), 
                                                   R = list(phylo=Terrestrial_A_cor), data = Terrestrial_Subset_Data, method = "REML", sparse = TRUE, 
                                                   control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Trait_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Terrestrial_Trait_Model_CVR.rds")
  } else {
            Terrestrial_Trait_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Terrestrial_Trait_Model_CVR.rds")})

Terrestrial_Trait_Model_CVR_rob <- robust(Terrestrial_Trait_Model_CVR, cluster = Terrestrial_Subset_Data$Study_ID, adjust = TRUE)

Terrestrial_Trait_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Terrestrial_Trait_Model_CVR$b), 15, 100),
                                                    estimate = Terrestrial_Trait_Model_CVR$b, 
                                                    ci.lb = Terrestrial_Trait_Model_CVR$ci.lb, 
                                                    ci.ub = Terrestrial_Trait_Model_CVR$ci.ub)
rownames(Terrestrial_Trait_Model_CVR_Estimates) <- Terrestrial_Trait_Model_CVR_Estimates$Category
Terrestrial_Trait_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Trait_Model_CVR), 2))

# Preparing Graph - Combined

terrestrial_trait_rnames <- c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits", 
                              "Morphology", "Physiological")

terrestrial_trait_k <- data.frame("k" = c(Terrestrial_Trait_Exploration["Behavioural", "Freq"], 
                                          Terrestrial_Trait_Exploration["Biochemical Assay", "Freq"], 
                                          Terrestrial_Trait_Exploration["Gene Expression", "Freq"], 
                                          Terrestrial_Trait_Exploration["Life-History Traits", "Freq"], 
                                          Terrestrial_Trait_Exploration["Morphology", "Freq"], 
                                          Terrestrial_Trait_Exploration["Physiological", "Freq"]), 
                                  row.names = terrestrial_trait_rnames)

terrestrial_trait_group_no <- data.frame("Spp No." = c(Terrestrial_Trait_Species_Count["Behavioural", "Freq"], 
                                                       Terrestrial_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                       Terrestrial_Trait_Species_Count["Gene Expression", "Freq"], 
                                                       Terrestrial_Trait_Species_Count["Life-History Traits", "Freq"],
                                                       Terrestrial_Trait_Species_Count["Morphology", "Freq"],
                                                       Terrestrial_Trait_Species_Count["Physiological", "Freq"]), 
                                         row.names = terrestrial_trait_rnames)

terrestrial_trait_study <- data.frame("Study" = c(Terrestrial_Trait_Study_Count["Behavioural", "Freq"], 
                                                  Terrestrial_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                  Terrestrial_Trait_Study_Count["Gene Expression", "Freq"], 
                                                  Terrestrial_Trait_Study_Count["Life-History Traits", "Freq"],
                                                  Terrestrial_Trait_Study_Count["Morphology", "Freq"],
                                                  Terrestrial_Trait_Study_Count["Physiological", "Freq"]), 
                                         row.names = terrestrial_trait_rnames)

terrestrial_trait_table <- data.frame(estimate = Terrestrial_Trait_Model_CVR_Estimates[,"estimate"], 
                                      lowerCL = Terrestrial_Trait_Model_CVR_Estimates[,"ci.lb"], 
                                      upperCL = Terrestrial_Trait_Model_CVR_Estimates[,"ci.ub"], 
                                      K = terrestrial_trait_k[,1], 
                                      group_no = terrestrial_trait_group_no[,1], 
                                      row.names = terrestrial_trait_rnames)
terrestrial_trait_table$name <- row.names(terrestrial_trait_table)

terrestrial_trait_raw_mean <- c(unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Behavioural") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Gene Expression") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                                select("InCVR"))),
                                unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                                select("InCVR"))))

terrestrial_trait_raw_name <- c(replicate(25, "Behavioural"), 
                                replicate(79, "Biochemical Assay"), 
                                replicate(12, "Gene Expression"), 
                                replicate(441, "Life-history Traits"), 
                                replicate(286, "Morphology"),
                                replicate(86, "Physiological"))

terrestrial_trait_raw_df <- data.frame("Model" = terrestrial_trait_raw_name, 
                                       "Effect" = terrestrial_trait_raw_mean)

# Graph code - Combined

Terrestrial_Trait_Order <- c("Physiological", "Morphology", "Life-history Traits",  
                             "Gene Expression", "Biochemical Assay", "Behavioural")

density_terrestrial_trait_CVR <- terrestrial_trait_table %>% mutate(name = fct_relevel(name, Terrestrial_Trait_Order)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = terrestrial_trait_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Trait_Order)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                 size = 1, fatten = 2) +
                                 theme_bw() +
                                 guides(fill = "none", colour = "none") +
                                 labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                 theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                  vjust = c(-2.7, -2.7, -0.8, -0.8, -0.8, -2.7))) +
                                 theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                 theme(axis.ticks = element_blank()) +
                                 theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                 theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                 scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                 scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                 scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                 coord_cartesian(xlim = c(-1, 1)) +
                                 annotate('text',  x = 1, y = (seq(1, dim(terrestrial_trait_table)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(terrestrial_trait_table["Physiological", "K"], 
                                                               terrestrial_trait_table["Morphology", "K"], 
                                                               terrestrial_trait_table["Life-history Traits", "K"],
                                                               terrestrial_trait_table["Gene Expression", "K"],
                                                               terrestrial_trait_table["Biochemical Assay", "K"],
                                                               terrestrial_trait_table["Behavioural", "K"]), "~","(", 
                                                             c(terrestrial_trait_table["Physiological", "group_no"], 
                                                               terrestrial_trait_table["Morphology", "group_no"], 
                                                               terrestrial_trait_table["Life-history Traits", "group_no"],
                                                               terrestrial_trait_table["Gene Expression", "group_no"],
                                                               terrestrial_trait_table["Biochemical Assay", "group_no"],
                                                               terrestrial_trait_table["Behavioural", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                 geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                        paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                        paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                        paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Gene Expression", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                        paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                        paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Behavioural", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.75, y = (seq(1, dim(terrestrial_trait_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_trait_CVR #(400x560)

# Preparing Graph - Part 1

terrestrial_trait_rnames_1 <- c("Behavioural", "Biochemical Assay", "Gene Expression")

terrestrial_trait_k_1 <- data.frame("k" = c(Terrestrial_Trait_Exploration["Behavioural", "Freq"], 
                                            Terrestrial_Trait_Exploration["Biochemical Assay", "Freq"], 
                                            Terrestrial_Trait_Exploration["Gene Expression", "Freq"]), 
                                    row.names = terrestrial_trait_rnames_1)

terrestrial_trait_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Trait_Species_Count["Behavioural", "Freq"], 
                                                         Terrestrial_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                         Terrestrial_Trait_Species_Count["Gene Expression", "Freq"]), 
                                           row.names = terrestrial_trait_rnames_1)

terrestrial_trait_study_1 <- data.frame("Study" = c(Terrestrial_Trait_Study_Count["Behavioural", "Freq"], 
                                                    Terrestrial_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                    Terrestrial_Trait_Study_Count["Gene Expression", "Freq"]), 
                                        row.names = terrestrial_trait_rnames_1)

Terrestrial_Trait_Model_CVR_Estimates_Reorder_1 <- Terrestrial_Trait_Model_CVR_Estimates[c("Behavioural", "Biochemical Assay", "Gene Expression"), ]

terrestrial_trait_table_1 <- data.frame(estimate = Terrestrial_Trait_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                      lowerCL = Terrestrial_Trait_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                      upperCL = Terrestrial_Trait_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                      K = terrestrial_trait_k_1[,1], 
                                      group_no = terrestrial_trait_group_no_1[,1], 
                                      row.names = terrestrial_trait_rnames_1)
terrestrial_trait_table_1$name <- row.names(terrestrial_trait_table_1)

terrestrial_trait_raw_mean_1 <- c(unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Behavioural") %>% 
                                                select("InCVR"))), 
                                  unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                                select("InCVR"))), 
                                  unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Gene Expression") %>% 
                                                select("InCVR"))))

terrestrial_trait_raw_name_1 <- c(replicate(25, "Behavioural"), 
                                  replicate(79, "Biochemical Assay"), 
                                  replicate(12, "Gene Expression"))

terrestrial_trait_raw_df_1 <- data.frame("Model" = terrestrial_trait_raw_name_1, 
                                         "Effect" = terrestrial_trait_raw_mean_1)

# Graph code - Part 1

Terrestrial_Trait_Order_1 <- c("Gene Expression", "Biochemical Assay", "Behavioural")

density_terrestrial_trait_CVR_1 <- terrestrial_trait_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Trait_Order_1)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = terrestrial_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Trait_Order_1)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                         vjust = c(-0.8, -0.8, -2.7))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                   scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(terrestrial_trait_table_1)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(terrestrial_trait_table_1["Gene Expression", "K"],
                                                                 terrestrial_trait_table_1["Biochemical Assay", "K"],
                                                                 terrestrial_trait_table_1["Behavioural", "K"]), "~","(", 
                                                               c(terrestrial_trait_table_1["Gene Expression", "group_no"],
                                                                 terrestrial_trait_table_1["Biochemical Assay", "group_no"],
                                                                 terrestrial_trait_table_1["Behavioural", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Gene Expression", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Behavioural", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(terrestrial_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_trait_CVR_1 #(400x320)

# Preparing Graph - Part 2

terrestrial_trait_rnames_2 <- c("Life-history Traits", "Morphology", "Physiological")

terrestrial_trait_k_2 <- data.frame("k" = c(Terrestrial_Trait_Exploration["Life-History Traits", "Freq"], 
                                            Terrestrial_Trait_Exploration["Morphology", "Freq"], 
                                            Terrestrial_Trait_Exploration["Physiological", "Freq"]), 
                                    row.names = terrestrial_trait_rnames_2)

terrestrial_trait_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Trait_Species_Count["Life-History Traits", "Freq"],
                                                         Terrestrial_Trait_Species_Count["Morphology", "Freq"],
                                                         Terrestrial_Trait_Species_Count["Physiological", "Freq"]), 
                                           row.names = terrestrial_trait_rnames_2)

terrestrial_trait_study_2 <- data.frame("Study" = c(Terrestrial_Trait_Study_Count["Life-History Traits", "Freq"],
                                                    Terrestrial_Trait_Study_Count["Morphology", "Freq"],
                                                    Terrestrial_Trait_Study_Count["Physiological", "Freq"]), 
                                        row.names = terrestrial_trait_rnames_2)

Terrestrial_Trait_Model_CVR_Estimates_Reorder_2 <- Terrestrial_Trait_Model_CVR_Estimates[c("Life-History Traits", "Morphology", "Physiological"), ]

terrestrial_trait_table_2 <- data.frame(estimate = Terrestrial_Trait_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                        lowerCL = Terrestrial_Trait_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                        upperCL = Terrestrial_Trait_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                        K = terrestrial_trait_k_2[,1], 
                                        group_no = terrestrial_trait_group_no_2[,1], 
                                        row.names = terrestrial_trait_rnames_2)
terrestrial_trait_table_2$name <- row.names(terrestrial_trait_table_2)

terrestrial_trait_raw_mean_2 <- c(unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                                select("InCVR"))), 
                                  unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                                select("InCVR"))),
                                  unlist(unname(Terrestrial_Subset_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                                select("InCVR"))))

terrestrial_trait_raw_name_2 <- c(replicate(441, "Life-history Traits"), 
                                  replicate(286, "Morphology"),
                                  replicate(86, "Physiological"))

terrestrial_trait_raw_df_2 <- data.frame("Model" = terrestrial_trait_raw_name_2, 
                                         "Effect" = terrestrial_trait_raw_mean_2)

# Graph code - Part 2

Terrestrial_Trait_Order_2 <- c("Physiological", "Morphology", "Life-history Traits")

density_terrestrial_trait_CVR_2 <- terrestrial_trait_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Trait_Order_2)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = terrestrial_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Trait_Order_2)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(terrestrial_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                         vjust = c(-2.7, -2.7, -0.8))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                                   scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(terrestrial_trait_table_2)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(terrestrial_trait_table_2["Physiological", "K"], 
                                                                 terrestrial_trait_table_2["Morphology", "K"], 
                                                                 terrestrial_trait_table_2["Life-history Traits", "K"]), "~","(", 
                                                               c(terrestrial_trait_table_2["Physiological", "group_no"], 
                                                                 terrestrial_trait_table_2["Morphology", "group_no"], 
                                                                 terrestrial_trait_table_2["Life-history Traits", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                          paste(format(round(mean(exp(Terrestrial_Trait_Model_CVR_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(terrestrial_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_trait_CVR_2 #(400x320)

##### Terrestrial Subset Model - Plasticity Mechanism Meta-Regression - CVR #####
Terrestrial_Plasticity_Exploration <- Terrestrial_Subset_Data %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Terrestrial_Plasticity_Exploration) <- Terrestrial_Plasticity_Exploration$Plasticity_Mechanism

Terrestrial_Plasticity_Species_Count <- Terrestrial_Subset_Data %>% select("Scientific_Name", "Plasticity_Mechanism") %>% table() %>% data.frame() %>%
                                        filter(`Freq` != 0) %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Terrestrial_Plasticity_Species_Count) <- Terrestrial_Plasticity_Species_Count$Plasticity_Mechanism

Terrestrial_Plasticity_Study_Count <- Terrestrial_Subset_Data %>% select("Study_ID", "Plasticity_Mechanism") %>% table() %>% data.frame() %>%
                                      filter(`Freq` != 0) %>% select("Plasticity_Mechanism") %>% table() %>% data.frame()
rownames(Terrestrial_Plasticity_Study_Count) <- Terrestrial_Plasticity_Study_Count$Plasticity_Mechanism

run <- FALSE
system.time( #  8ish minutes
  if(run){
    Terrestrial_Plasticity_Model_CVR <- metafor::rma.mv(InCVR, V = Terrestrial_VCV_InCVR, test = "t", dfs = "contain",
                                                        mods = ~ Plasticity_Mechanism - 1,
                                                        random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                      ~1|Shared_Animal_Number, ~1|Measurement), 
                                                        R = list(phylo=Terrestrial_A_cor), data = Terrestrial_Subset_Data, method = "REML", sparse = TRUE, 
                                                        control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Plasticity_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Terrestrial_Plasticity_Model_CVR.rds")
  } else {
            Terrestrial_Plasticity_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Terrestrial_Plasticity_Model_CVR.rds")})

Terrestrial_Plasticity_Model_CVR_rob <- robust(Terrestrial_Plasticity_Model_CVR, cluster = Terrestrial_Subset_Data$Study_ID, adjust = TRUE)

Terrestrial_Plasticity_Model_CVR_Estimates <- data.frame(Plasticity_Mechanism = substr(row.names(Terrestrial_Plasticity_Model_CVR$b), 21, 100),
                                                         estimate = Terrestrial_Plasticity_Model_CVR$b, 
                                                         ci.lb = Terrestrial_Plasticity_Model_CVR$ci.lb, 
                                                         ci.ub = Terrestrial_Plasticity_Model_CVR$ci.ub)
rownames(Terrestrial_Plasticity_Model_CVR_Estimates) <- Terrestrial_Plasticity_Model_CVR_Estimates$Plasticity_Mechanism
Terrestrial_Plasticity_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Plasticity_Model_CVR), 2))

# Preparing Graph - Combined

terrestrial_plasticity_rnames <- c("Acclimation", "Development")

terrestrial_plasticity_k <- data.frame("k" = c(Terrestrial_Plasticity_Exploration["Acclimation", "Freq"], 
                                               Terrestrial_Plasticity_Exploration["Developmental Plasticity", "Freq"]), 
                                       row.names = terrestrial_plasticity_rnames)

terrestrial_plasticity_group_no <- data.frame("Spp No." = c(Terrestrial_Plasticity_Species_Count["Acclimation", "Freq"], 
                                                            Terrestrial_Plasticity_Species_Count["Developmental Plasticity", "Freq"]), 
                                              row.names = terrestrial_plasticity_rnames)

terrestrial_plasticity_study <- data.frame("Study" = c(Terrestrial_Plasticity_Study_Count["Acclimation", "Freq"], 
                                                       Terrestrial_Plasticity_Study_Count["Developmental Plasticity", "Freq"]), 
                                              row.names = terrestrial_plasticity_rnames)

terrestrial_plasticity_table <- data.frame(estimate = Terrestrial_Plasticity_Model_CVR_Estimates[,"estimate"], 
                                           lowerCL = Terrestrial_Plasticity_Model_CVR_Estimates[,"ci.lb"], 
                                           upperCL = Terrestrial_Plasticity_Model_CVR_Estimates[,"ci.ub"], 
                                           K = terrestrial_plasticity_k[,1], 
                                           group_no = terrestrial_plasticity_group_no[,1], 
                                           row.names = terrestrial_plasticity_rnames)
terrestrial_plasticity_table$name <- row.names(terrestrial_plasticity_table)

terrestrial_plasticity_raw_mean <- c(unlist(unname(Terrestrial_Subset_Data %>% filter(`Plasticity_Mechanism` == "Acclimation") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Terrestrial_Subset_Data %>% filter(`Plasticity_Mechanism` == "Developmental Plasticity") %>% 
                                                     select("InCVR"))))

terrestrial_plasticity_raw_name <- c(replicate(139, "Acclimation"), 
                                     replicate(790, "Development"))

terrestrial_plasticity_raw_df <- data.frame("Model" = terrestrial_plasticity_raw_name, 
                                            "Effect" = terrestrial_plasticity_raw_mean)

# Graph code - Combined

Terrestrial_Plasticity_Order <- c("Development", "Acclimation")

density_terrestrial_plasticity_CVR <- terrestrial_plasticity_table %>% mutate(name = fct_relevel(name, Terrestrial_Plasticity_Order)) %>%
                                      ggplot() +
                                      geom_density_ridges(data = terrestrial_plasticity_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Plasticity_Order)), 
                                                          aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                          scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                      geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                     size = 1) +
                                      geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_plasticity_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                      size = 1, fatten = 2) +
                                      theme_bw() +
                                      guides(fill = "none", colour = "none") +
                                      labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                      theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                       vjust = c(-2.7, -2.7))) +
                                      theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                      theme(axis.ticks = element_blank()) +
                                      theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                      scale_colour_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                      scale_fill_manual(values = c("#5D7AA1", "#2B4E7A")) +
                                      coord_cartesian(xlim = c(-1, 1)) +
                                      annotate('text',  x = 1, y = (seq(1, dim(terrestrial_plasticity_table)[1], 1)+0.4),
                                      label= paste("italic(k)==", c(terrestrial_plasticity_table["Development", "K"], 
                                                                    terrestrial_plasticity_table["Acclimation", "K"]), "~","(", 
                                                                  c(terrestrial_plasticity_table["Development", "group_no"], 
                                                                    terrestrial_plasticity_table["Acclimation", "group_no"]), 
                                                   ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                      geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Plasticity_Model_CVR_Estimates["Development", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                             paste(format(round(mean(exp(Terrestrial_Plasticity_Model_CVR_Estimates["Acclimation", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                 x = -0.75, y = (seq(1, dim(terrestrial_plasticity_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_plasticity_CVR #(400x240)

# Preparing Graph - Part 1

terrestrial_plasticity_rnames_1 <- c("Acclimation")

terrestrial_plasticity_k_1 <- data.frame("k" = c(Terrestrial_Plasticity_Exploration["Acclimation", "Freq"]), 
                                         row.names = terrestrial_plasticity_rnames_1)

terrestrial_plasticity_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Plasticity_Species_Count["Acclimation", "Freq"]), 
                                                row.names = terrestrial_plasticity_rnames_1)

terrestrial_plasticity_study_1 <- data.frame("Study" = c(Terrestrial_Plasticity_Study_Count["Acclimation", "Freq"]), 
                                             row.names = terrestrial_plasticity_rnames_1)

Terrestrial_Plasticity_Model_CVR_Estimates_Reorder_1 <- Terrestrial_Plasticity_Model_CVR_Estimates[c("Acclimation"), ]

terrestrial_plasticity_table_1 <- data.frame(estimate = Terrestrial_Plasticity_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                             lowerCL = Terrestrial_Plasticity_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                             upperCL = Terrestrial_Plasticity_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                             K = terrestrial_plasticity_k_1[,1], 
                                             group_no = terrestrial_plasticity_group_no_1[,1], 
                                             row.names = terrestrial_plasticity_rnames_1)
terrestrial_plasticity_table_1$name <- row.names(terrestrial_plasticity_table_1)

terrestrial_plasticity_raw_mean_1 <- c(unlist(unname(Terrestrial_Subset_Data %>% filter(`Plasticity_Mechanism` == "Acclimation") %>% 
                                                     select("InCVR"))))

terrestrial_plasticity_raw_name_1 <- c(replicate(139, "Acclimation"))

terrestrial_plasticity_raw_df_1 <- data.frame("Model" = terrestrial_plasticity_raw_name_1, 
                                              "Effect" = terrestrial_plasticity_raw_mean_1)

# Graph code - Part 1

Terrestrial_Plasticity_Order_1 <- c("Acclimation")

density_terrestrial_plasticity_CVR_1 <- terrestrial_plasticity_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Plasticity_Order_1)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = terrestrial_plasticity_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Plasticity_Order_1)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.13, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_plasticity_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                        size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                              vjust = c(-2.7))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#2B4E7A")) +
                                        scale_fill_manual(values = c("#2B4E7A")) +
                                        coord_cartesian(xlim = c(-1, 1)) +
                                        annotate('text',  x = 1, y = (seq(1, dim(terrestrial_plasticity_table_1)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(terrestrial_plasticity_table_1["Acclimation", "K"]), "~","(", 
                                                                    c(terrestrial_plasticity_table_1["Acclimation", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Plasticity_Model_CVR_Estimates["Acclimation", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.75, y = (seq(1, dim(terrestrial_plasticity_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_plasticity_CVR_1 #(400x160)

# Preparing Graph - Part 2

terrestrial_plasticity_rnames_2 <- c("Development")

terrestrial_plasticity_k_2 <- data.frame("k" = c(Terrestrial_Plasticity_Exploration["Developmental Plasticity", "Freq"]), 
                                         row.names = terrestrial_plasticity_rnames_2)

terrestrial_plasticity_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Plasticity_Species_Count["Developmental Plasticity", "Freq"]), 
                                                row.names = terrestrial_plasticity_rnames_2)

terrestrial_plasticity_study_2 <- data.frame("Study" = c(Terrestrial_Plasticity_Study_Count["Developmental Plasticity", "Freq"]), 
                                             row.names = terrestrial_plasticity_rnames_2)

Terrestrial_Plasticity_Model_CVR_Estimates_Reorder_2 <- Terrestrial_Plasticity_Model_CVR_Estimates[c("Development"), ]

terrestrial_plasticity_table_2 <- data.frame(estimate = Terrestrial_Plasticity_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                             lowerCL = Terrestrial_Plasticity_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                             upperCL = Terrestrial_Plasticity_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                             K = terrestrial_plasticity_k_2[,1], 
                                             group_no = terrestrial_plasticity_group_no_2[,1], 
                                             row.names = terrestrial_plasticity_rnames_2)
terrestrial_plasticity_table_2$name <- row.names(terrestrial_plasticity_table_2)

terrestrial_plasticity_raw_mean_2 <- c(unlist(unname(Terrestrial_Subset_Data %>% filter(`Plasticity_Mechanism` == "Developmental Plasticity") %>% 
                                                     select("InCVR"))))

terrestrial_plasticity_raw_name_2 <- c(replicate(790, "Development"))

terrestrial_plasticity_raw_df_2 <- data.frame("Model" = terrestrial_plasticity_raw_name_2, 
                                              "Effect" = terrestrial_plasticity_raw_mean_2)

# Graph code - Part 2

Terrestrial_Plasticity_Order_2 <- c("Development")

density_terrestrial_plasticity_CVR_2 <- terrestrial_plasticity_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Plasticity_Order_2)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = terrestrial_plasticity_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Plasticity_Order_2)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.13, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(terrestrial_plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_plasticity_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                        size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                              vjust = c(-2.7))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#5D7AA1")) +
                                        scale_fill_manual(values = c("#5D7AA1")) +
                                        coord_cartesian(xlim = c(-1, 1)) +
                                        annotate('text',  x = 1, y = (seq(1, dim(terrestrial_plasticity_table_2)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(terrestrial_plasticity_table_2["Development", "K"]), "~","(", 
                                                                    c(terrestrial_plasticity_table_2["Development", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Plasticity_Model_CVR_Estimates["Development", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.75, y = (seq(1, dim(terrestrial_plasticity_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_plasticity_CVR_2 #(400x160)

##### Terrestrial Subset Model - Specific Trait Meta-Regression - CVR #####
Terrestrial_Specific_Trait_Exploration <- Terrestrial_Subset_Data %>% select("Measurement") %>% table() %>% data.frame()
Terrestrial_Specific_Trait_Exploration <- Terrestrial_Specific_Trait_Exploration %>% filter(Freq > 10)
rownames(Terrestrial_Specific_Trait_Exploration) <- Terrestrial_Specific_Trait_Exploration$Measurement

Terrestrial_Specific_Trait_Data <- Terrestrial_Subset_Data %>% filter(Measurement == "Development Time"| 
                                                                      Measurement == "Fecundity"|
                                                                      Measurement == "Food Consumption"|
                                                                      Measurement == "Head Width"|
                                                                      Measurement == "Length"|
                                                                      Measurement == "Locomotor Performance"|
                                                                      Measurement == "Longevity"|
                                                                      Measurement == "Mass"|
                                                                      Measurement == "Metabolic Rate"|
                                                                      Measurement == "PO Activity"|
                                                                      Measurement == "Reproductive Rate"|
                                                                      Measurement == "Tail Length")

Terrestrial_Specific_Trait_Species_Count <- Terrestrial_Specific_Trait_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>%
                                            filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Terrestrial_Specific_Trait_Species_Count) <- Terrestrial_Specific_Trait_Species_Count$Measurement

Terrestrial_Specific_Trait_Study_Count <- Terrestrial_Specific_Trait_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>%
                                          filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Terrestrial_Specific_Trait_Study_Count) <- Terrestrial_Specific_Trait_Study_Count$Measurement

Terrestrial_Specific_Trait_Species <- Terrestrial_Specific_Trait_Data %>% select("phylo") %>% unique()

Terrestrial_Specific_Trait_A_cor <- as.data.frame(A_cor)
Terrestrial_Specific_Trait_A_cor <- Terrestrial_Specific_Trait_A_cor[c(Terrestrial_Specific_Trait_Species$phylo), c(Terrestrial_Specific_Trait_Species$phylo)]
Terrestrial_Specific_Trait_A_cor <- as.matrix(Terrestrial_Specific_Trait_A_cor)

Terrestrial_Specific_Trait_VCV_InCVR <- make_VCV_matrix(Terrestrial_Specific_Trait_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  3ish minutes
  if(run){
    Terrestrial_Specific_Trait_Model_CVR <- metafor::rma.mv(InCVR, V = Terrestrial_Specific_Trait_VCV_InCVR, test = "t", dfs = "contain",
                                                            mods = ~ Measurement - 1,
                                                            random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                          ~1|Shared_Animal_Number), 
                                                            R = list(phylo=Terrestrial_Specific_Trait_A_cor), data = Terrestrial_Specific_Trait_Data, method = "REML", sparse = TRUE, 
                                                            control=list(rel.tol=1e-9))
    saveRDS(Terrestrial_Specific_Trait_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Terrestrial_Specific_Trait_Model_CVR.rds")
  } else {
            Terrestrial_Specific_Trait_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Terrestrial_Specific_Trait_Model_CVR.rds")})

Terrestrial_Specific_Trait_Model_CVR_rob <- robust(Terrestrial_Specific_Trait_Model_CVR, cluster = Terrestrial_Specific_Trait_Data$Study_ID, adjust = TRUE)

Terrestrial_Specific_Trait_Model_CVR_Estimates <- data.frame(Trait = substr(row.names(Terrestrial_Specific_Trait_Model_CVR$b), 12, 100),
                                                             estimate = Terrestrial_Specific_Trait_Model_CVR$b, 
                                                             ci.lb = Terrestrial_Specific_Trait_Model_CVR$ci.lb, 
                                                             ci.ub = Terrestrial_Specific_Trait_Model_CVR$ci.ub)
rownames(Terrestrial_Specific_Trait_Model_CVR_Estimates) <- Terrestrial_Specific_Trait_Model_CVR_Estimates$Trait
Terrestrial_Specific_Trait_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Terrestrial_Specific_Trait_Model_CVR), 2))

# Preparing Graph - Combined

terrestrial_specific_trait_rnames <- c("Development Time", "Fecundity", "Food Consumption", "Head Width", 
                                       "Length", "Locomotor Performance", "Longevity", "Mass", 
                                       "Metabolic Rate", "PO Activity", "Reproductive Rate", "Tail Length")

terrestrial_specific_trait_k <- data.frame("k" = c(Terrestrial_Specific_Trait_Exploration["Development Time", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Fecundity", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Food Consumption", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Head Width", "Freq"],
                                                   Terrestrial_Specific_Trait_Exploration["Length", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Locomotor Performance", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Longevity", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Mass", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Metabolic Rate", "Freq"],  
                                                   Terrestrial_Specific_Trait_Exploration["PO Activity", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Reproductive Rate", "Freq"], 
                                                   Terrestrial_Specific_Trait_Exploration["Tail Length", "Freq"]), 
                                           row.names = terrestrial_specific_trait_rnames)

terrestrial_specific_trait_group_no <- data.frame("Spp No." = c(Terrestrial_Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                                Terrestrial_Specific_Trait_Species_Count["Fecundity", "Freq"], 
                                                                Terrestrial_Specific_Trait_Species_Count["Food Consumption", "Freq"],
                                                                Terrestrial_Specific_Trait_Species_Count["Head Width", "Freq"],
                                                                Terrestrial_Specific_Trait_Species_Count["Length", "Freq"], 
                                                                Terrestrial_Specific_Trait_Species_Count["Locomotor Performance", "Freq"], 
                                                                Terrestrial_Specific_Trait_Species_Count["Longevity", "Freq"], 
                                                                Terrestrial_Specific_Trait_Species_Count["Mass", "Freq"], 
                                                                Terrestrial_Specific_Trait_Species_Count["Metabolic Rate", "Freq"],  
                                                                Terrestrial_Specific_Trait_Species_Count["PO Activity", "Freq"], 
                                                                Terrestrial_Specific_Trait_Species_Count["Reproductive Rate", "Freq"], 
                                                                Terrestrial_Specific_Trait_Species_Count["Tail Length", "Freq"]), 
                                                  row.names = terrestrial_specific_trait_rnames)

terrestrial_specific_trait_study <- data.frame("Study" = c(Terrestrial_Specific_Trait_Study_Count["Development Time", "Freq"], 
                                                           Terrestrial_Specific_Trait_Study_Count["Fecundity", "Freq"], 
                                                           Terrestrial_Specific_Trait_Study_Count["Food Consumption", "Freq"],
                                                           Terrestrial_Specific_Trait_Study_Count["Head Width", "Freq"],
                                                           Terrestrial_Specific_Trait_Study_Count["Length", "Freq"], 
                                                           Terrestrial_Specific_Trait_Study_Count["Locomotor Performance", "Freq"], 
                                                           Terrestrial_Specific_Trait_Study_Count["Longevity", "Freq"], 
                                                           Terrestrial_Specific_Trait_Study_Count["Mass", "Freq"], 
                                                           Terrestrial_Specific_Trait_Study_Count["Metabolic Rate", "Freq"],  
                                                           Terrestrial_Specific_Trait_Study_Count["PO Activity", "Freq"], 
                                                           Terrestrial_Specific_Trait_Study_Count["Reproductive Rate", "Freq"], 
                                                           Terrestrial_Specific_Trait_Study_Count["Tail Length", "Freq"]), 
                                                  row.names = terrestrial_specific_trait_rnames)

terrestrial_specific_trait_table <- data.frame(estimate = Terrestrial_Specific_Trait_Model_CVR_Estimates[,"estimate"], 
                                               lowerCL = Terrestrial_Specific_Trait_Model_CVR_Estimates[,"ci.lb"], 
                                               upperCL = Terrestrial_Specific_Trait_Model_CVR_Estimates[,"ci.ub"], 
                                               K = terrestrial_specific_trait_k[,1], 
                                               group_no = terrestrial_specific_trait_group_no[,1], 
                                               row.names = terrestrial_specific_trait_rnames)
terrestrial_specific_trait_table$name <- row.names(terrestrial_specific_trait_table)

terrestrial_specific_trait_raw_mean <- c(unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Fecundity") %>% 
                                                         select("InCVR"))),
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Food Consumption") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Head Width") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Locomotor Performance") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Longevity") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Metabolic Rate") %>% 
                                                         select("InCVR"))),
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "PO Activity") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Reproductive Rate") %>% 
                                                         select("InCVR"))),
                                         unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Tail Length") %>% 
                                                         select("InCVR"))))

terrestrial_specific_trait_raw_name <- c(replicate(259, "Development Time"), 
                                         replicate(67, "Fecundity"),
                                         replicate(11, "Food Consumption"), 
                                         replicate(14, "Head Width"), 
                                         replicate(46, "Length"), 
                                         replicate(24, "Locomotor Performance"), 
                                         replicate(86, "Longevity"), 
                                         replicate(96, "Mass"), 
                                         replicate(22, "Metabolic Rate"),
                                         replicate(14, "PO Activity"), 
                                         replicate(11, "Reproductive Rate"), 
                                         replicate(21, "Tail Length"))

terrestrial_specific_trait_raw_df <- data.frame("Model" = terrestrial_specific_trait_raw_name, 
                                                "Effect" = terrestrial_specific_trait_raw_mean)

# Graph code - Combined

Terrestrial_Specific_Trait_Order <- c("Tail Length", "Reproductive Rate", "PO Activity", "Metabolic Rate", 
                                      "Mass", "Longevity", "Locomotor Performance", "Length", "Head Width", "Food Consumption", 
                                      "Fecundity", "Development Time")

density_terrestrial_specific_trait_CVR <- terrestrial_specific_trait_table %>% mutate(name = fct_relevel(name, Terrestrial_Specific_Trait_Order)) %>%
                                          ggplot() +
                                          geom_density_ridges(data = terrestrial_specific_trait_raw_df %>% mutate(Model = fct_relevel(Model, Terrestrial_Specific_Trait_Order)), 
                                                              aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                              scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                          geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                         size = 1) +
                                          geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                          size = 1, fatten = 2) +
                                          theme_bw() +
                                          guides(fill = "none", colour = "none") +
                                          labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                          theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                           vjust = c(-2.7, -0.8, -2.7, -0.8, -2.7, -2.7, 
                                                                                     -0.8, -2.7, -2.7, -0.8, -2.7, -0.8))) +
                                          theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                          theme(axis.ticks = element_blank()) +
                                          theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                          theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                          scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                          scale_colour_manual(values = c("#6C87AB", "#6582A9", "#607C9F", "#5777A1", "#4A6E9C", "#446692",
                                                                         "#3D5E89", "#355784", "#234373","#1C375F", "#0D2A51","#0F2643")) +
                                          scale_fill_manual(values = c("#6C87AB", "#6582A9", "#607C9F", "#5777A1", "#4A6E9C", "#446692",
                                                                       "#3D5E89", "#355784", "#234373","#1C375F", "#0D2A51","#0F2643")) +
                                          coord_cartesian(xlim = c(-1, 1)) +
                                          annotate('text',  x = 1, y = (seq(1, dim(terrestrial_specific_trait_table)[1], 1)+0.4),
                                          label= paste("italic(k)==", c(terrestrial_specific_trait_table["Tail Length", "K"],
                                                                        terrestrial_specific_trait_table["Reproductive Rate", "K"],
                                                                        terrestrial_specific_trait_table["PO Activity", "K"],
                                                                        terrestrial_specific_trait_table["Metabolic Rate", "K"],
                                                                        terrestrial_specific_trait_table["Mass", "K"],
                                                                        terrestrial_specific_trait_table["Longevity", "K"], 
                                                                        terrestrial_specific_trait_table["Locomotor Performance", "K"], 
                                                                        terrestrial_specific_trait_table["Length", "K"],
                                                                        terrestrial_specific_trait_table["Head Width", "K"],
                                                                        terrestrial_specific_trait_table["Food Consumption", "K"], 
                                                                        terrestrial_specific_trait_table["Fecundity", "K"],
                                                                        terrestrial_specific_trait_table["Development Time", "K"]), "~","(", 
                                                                      c(terrestrial_specific_trait_table["Tail Length", "group_no"],
                                                                        terrestrial_specific_trait_table["Reproductive Rate", "group_no"],
                                                                        terrestrial_specific_trait_table["PO Activity", "group_no"], 
                                                                        terrestrial_specific_trait_table["Metabolic Rate", "group_no"],
                                                                        terrestrial_specific_trait_table["Mass", "group_no"],
                                                                        terrestrial_specific_trait_table["Longevity", "group_no"], 
                                                                        terrestrial_specific_trait_table["Locomotor Performance", "group_no"], 
                                                                        terrestrial_specific_trait_table["Length", "group_no"], 
                                                                        terrestrial_specific_trait_table["Head Width", "group_no"],
                                                                        terrestrial_specific_trait_table["Food Consumption", "group_no"], 
                                                                        terrestrial_specific_trait_table["Fecundity", "group_no"],
                                                                        terrestrial_specific_trait_table["Development Time", "group_no"]), 
                                                       ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                          geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Tail Length", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Reproductive Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["PO Activity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Longevity", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Locomotor Performance", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Head Width", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Food Consumption", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Fecundity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                     x = -0.75, y = (seq(1, dim(terrestrial_specific_trait_table)[1], 1)+0.4)), size = 3.5)

density_terrestrial_specific_trait_CVR #(400x1040)

# Preparing Graph - Part 1

terrestrial_specific_trait_rnames_1 <- c("Development Time", "Fecundity", "Food Consumption", "Head Width", "Length", "Locomotor Performance")

terrestrial_specific_trait_k_1 <- data.frame("k" = c(Terrestrial_Specific_Trait_Exploration["Development Time", "Freq"], 
                                                     Terrestrial_Specific_Trait_Exploration["Fecundity", "Freq"], 
                                                     Terrestrial_Specific_Trait_Exploration["Food Consumption", "Freq"], 
                                                     Terrestrial_Specific_Trait_Exploration["Head Width", "Freq"],
                                                     Terrestrial_Specific_Trait_Exploration["Length", "Freq"], 
                                                     Terrestrial_Specific_Trait_Exploration["Locomotor Performance", "Freq"]), 
                                             row.names = terrestrial_specific_trait_rnames_1)

terrestrial_specific_trait_group_no_1 <- data.frame("Spp No." = c(Terrestrial_Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                                  Terrestrial_Specific_Trait_Species_Count["Fecundity", "Freq"], 
                                                                  Terrestrial_Specific_Trait_Species_Count["Food Consumption", "Freq"],
                                                                  Terrestrial_Specific_Trait_Species_Count["Head Width", "Freq"],
                                                                  Terrestrial_Specific_Trait_Species_Count["Length", "Freq"], 
                                                                  Terrestrial_Specific_Trait_Species_Count["Locomotor Performance", "Freq"]), 
                                                    row.names = terrestrial_specific_trait_rnames_1)

terrestrial_specific_trait_study_1 <- data.frame("Study" = c(Terrestrial_Specific_Trait_Study_Count["Development Time", "Freq"], 
                                                             Terrestrial_Specific_Trait_Study_Count["Fecundity", "Freq"], 
                                                             Terrestrial_Specific_Trait_Study_Count["Food Consumption", "Freq"],
                                                             Terrestrial_Specific_Trait_Study_Count["Head Width", "Freq"],
                                                             Terrestrial_Specific_Trait_Study_Count["Length", "Freq"], 
                                                             Terrestrial_Specific_Trait_Study_Count["Locomotor Performance", "Freq"]), 
                                                 row.names = terrestrial_specific_trait_rnames_1)

Terrestrial_Specific_Trait_Model_CVR_Estimates_Reorder_1 <- Terrestrial_Specific_Trait_Model_CVR_Estimates[c("Development Time", "Fecundity", "Food Consumption", "Head Width", "Length", "Locomotor Performance"), ]

terrestrial_specific_trait_table_1 <- data.frame(estimate = Terrestrial_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                                 lowerCL = Terrestrial_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                                 upperCL = Terrestrial_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                                 K = terrestrial_specific_trait_k_1[,1], 
                                                 group_no = terrestrial_specific_trait_group_no_1[,1], 
                                                 row.names = terrestrial_specific_trait_rnames_1)
terrestrial_specific_trait_table_1$name <- row.names(terrestrial_specific_trait_table_1)

terrestrial_specific_trait_raw_mean_1 <- c(unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                                         select("InCVR"))), 
                                          unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Fecundity") %>% 
                                                         select("InCVR"))),
                                          unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Food Consumption") %>% 
                                                         select("InCVR"))), 
                                          unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Head Width") %>% 
                                                         select("InCVR"))), 
                                          unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                                         select("InCVR"))), 
                                          unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Locomotor Performance") %>% 
                                                         select("InCVR"))))

terrestrial_specific_trait_raw_name_1 <- c(replicate(259, "Development Time"), 
                                           replicate(67, "Fecundity"),
                                           replicate(11, "Food Consumption"), 
                                           replicate(14, "Head Width"), 
                                           replicate(46, "Length"), 
                                           replicate(24, "Locomotor Performance"))

terrestrial_specific_trait_raw_df_1 <- data.frame("Model" = terrestrial_specific_trait_raw_name_1, 
                                                  "Effect" = terrestrial_specific_trait_raw_mean_1)

# Graph code - Part 1

Terrestrial_Specific_Trait_Order_1 <- c("Locomotor Performance", "Length", "Head Width", "Food Consumption", "Fecundity", "Development Time")

density_terrestrial_specific_trait_CVR_1 <- terrestrial_specific_trait_table_1 %>% mutate(name = fct_relevel(name, Terrestrial_Specific_Trait_Order_1)) %>%
                                            ggplot() +
                                            geom_density_ridges(data = terrestrial_specific_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Terrestrial_Specific_Trait_Order_1)), 
                                                                aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                            geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                           size = 1) +
                                            geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                            size = 1, fatten = 2) +
                                            theme_bw() +
                                            guides(fill = "none", colour = "none") +
                                            labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                            theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                  vjust = c(-0.8, -2.7, -2.7, -0.8, -2.7, -0.8))) +
                                            theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                            theme(axis.ticks = element_blank()) +
                                            theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                            theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                            scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                            scale_colour_manual(values = c("#3D5E89", "#355784", "#234373","#1C375F", "#0D2A51","#0F2643")) +
                                            scale_fill_manual(values = c("#3D5E89", "#355784", "#234373","#1C375F", "#0D2A51","#0F2643")) +
                                            coord_cartesian(xlim = c(-1, 1)) +
                                            annotate('text',  x = 1, y = (seq(1, dim(terrestrial_specific_trait_table_1)[1], 1)+0.4),
                                            label= paste("italic(k)==", c(terrestrial_specific_trait_table_1["Locomotor Performance", "K"], 
                                                                          terrestrial_specific_trait_table_1["Length", "K"],
                                                                          terrestrial_specific_trait_table_1["Head Width", "K"],
                                                                          terrestrial_specific_trait_table_1["Food Consumption", "K"], 
                                                                          terrestrial_specific_trait_table_1["Fecundity", "K"],
                                                                          terrestrial_specific_trait_table_1["Development Time", "K"]), "~","(", 
                                                                        c(terrestrial_specific_trait_table_1["Locomotor Performance", "group_no"], 
                                                                          terrestrial_specific_trait_table_1["Length", "group_no"], 
                                                                          terrestrial_specific_trait_table_1["Head Width", "group_no"],
                                                                          terrestrial_specific_trait_table_1["Food Consumption", "group_no"], 
                                                                          terrestrial_specific_trait_table_1["Fecundity", "group_no"],
                                                                          terrestrial_specific_trait_table_1["Development Time", "group_no"]), 
                                                       ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                                geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Locomotor Performance", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                       paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                       paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Head Width", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                       paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Food Consumption", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                       paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Fecundity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                       paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                          x = -0.75, y = (seq(1, dim(terrestrial_specific_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_terrestrial_specific_trait_CVR_1 #(400x540)

# Preparing Graph - Part 2

terrestrial_specific_trait_rnames_2 <- c("Longevity", "Mass", "Metabolic Rate", "PO Activity", "Reproductive Rate", "Tail Length")

terrestrial_specific_trait_k_2 <- data.frame("k" = c(Terrestrial_Specific_Trait_Exploration["Longevity", "Freq"], 
                                                     Terrestrial_Specific_Trait_Exploration["Mass", "Freq"], 
                                                     Terrestrial_Specific_Trait_Exploration["Metabolic Rate", "Freq"],  
                                                     Terrestrial_Specific_Trait_Exploration["PO Activity", "Freq"], 
                                                     Terrestrial_Specific_Trait_Exploration["Reproductive Rate", "Freq"], 
                                                     Terrestrial_Specific_Trait_Exploration["Tail Length", "Freq"]), 
                                             row.names = terrestrial_specific_trait_rnames_2)

terrestrial_specific_trait_group_no_2 <- data.frame("Spp No." = c(Terrestrial_Specific_Trait_Species_Count["Longevity", "Freq"], 
                                                                  Terrestrial_Specific_Trait_Species_Count["Mass", "Freq"], 
                                                                  Terrestrial_Specific_Trait_Species_Count["Metabolic Rate", "Freq"],  
                                                                  Terrestrial_Specific_Trait_Species_Count["PO Activity", "Freq"], 
                                                                  Terrestrial_Specific_Trait_Species_Count["Reproductive Rate", "Freq"], 
                                                                  Terrestrial_Specific_Trait_Species_Count["Tail Length", "Freq"]), 
                                                    row.names = terrestrial_specific_trait_rnames_2)

terrestrial_specific_trait_study_2 <- data.frame("Study" = c(Terrestrial_Specific_Trait_Study_Count["Longevity", "Freq"], 
                                                             Terrestrial_Specific_Trait_Study_Count["Mass", "Freq"], 
                                                             Terrestrial_Specific_Trait_Study_Count["Metabolic Rate", "Freq"],  
                                                             Terrestrial_Specific_Trait_Study_Count["PO Activity", "Freq"], 
                                                             Terrestrial_Specific_Trait_Study_Count["Reproductive Rate", "Freq"], 
                                                             Terrestrial_Specific_Trait_Study_Count["Tail Length", "Freq"]), 
                                                 row.names = terrestrial_specific_trait_rnames_2)

Terrestrial_Specific_Trait_Model_CVR_Estimates_Reorder_2 <- Terrestrial_Specific_Trait_Model_CVR_Estimates[c("Longevity", "Mass", "Metabolic Rate", "PO Activity", "Reproductive Rate", "Tail Length"), ]

terrestrial_specific_trait_table_2 <- data.frame(estimate = Terrestrial_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                                 lowerCL = Terrestrial_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                                 upperCL = Terrestrial_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                                 K = terrestrial_specific_trait_k_2[,1], 
                                                 group_no = terrestrial_specific_trait_group_no_2[,1], 
                                                 row.names = terrestrial_specific_trait_rnames_2)
terrestrial_specific_trait_table_2$name <- row.names(terrestrial_specific_trait_table_2)

terrestrial_specific_trait_raw_mean_2 <- c(unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Longevity") %>% 
                                                         select("InCVR"))), 
                                           unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                                         select("InCVR"))), 
                                           unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Metabolic Rate") %>% 
                                                         select("InCVR"))),
                                           unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "PO Activity") %>% 
                                                         select("InCVR"))), 
                                           unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Reproductive Rate") %>% 
                                                         select("InCVR"))),
                                           unlist(unname(Terrestrial_Specific_Trait_Data %>% filter(`Measurement` == "Tail Length") %>% 
                                                         select("InCVR"))))

terrestrial_specific_trait_raw_name_2 <- c(replicate(86, "Longevity"), 
                                           replicate(96, "Mass"), 
                                           replicate(22, "Metabolic Rate"),
                                           replicate(14, "PO Activity"), 
                                           replicate(11, "Reproductive Rate"), 
                                           replicate(21, "Tail Length"))

terrestrial_specific_trait_raw_df_2 <- data.frame("Model" = terrestrial_specific_trait_raw_name_2, 
                                                  "Effect" = terrestrial_specific_trait_raw_mean_2)

# Graph code - Part 2

Terrestrial_Specific_Trait_Order_2 <- c("Tail Length", "Reproductive Rate", "PO Activity", "Metabolic Rate", 
                                        "Mass", "Longevity")

density_terrestrial_specific_trait_CVR_2 <- terrestrial_specific_trait_table_2 %>% mutate(name = fct_relevel(name, Terrestrial_Specific_Trait_Order_2)) %>%
                                            ggplot() +
                                            geom_density_ridges(data = terrestrial_specific_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Terrestrial_Specific_Trait_Order_2)), 
                                                                aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                            geom_linerange(aes(y = rev(seq(1, dim(terrestrial_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                           size = 1) +
                                            geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(terrestrial_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                            size = 1, fatten = 2) +
                                            theme_bw() +
                                            guides(fill = "none", colour = "none") +
                                            labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                            theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                  vjust = c(-2.7, -0.8, -2.7, -0.8, -2.7, -2.7))) +
                                            theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                            theme(axis.ticks = element_blank()) +
                                            theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                            theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                            scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                            scale_colour_manual(values = c("#6C87AB", "#6582A9", "#607C9F", "#5777A1", "#4A6E9C", "#446692")) +
                                            scale_fill_manual(values = c("#6C87AB", "#6582A9", "#607C9F", "#5777A1", "#4A6E9C", "#446692")) +
                                            coord_cartesian(xlim = c(-1, 1)) +
                                            annotate('text',  x = 1, y = (seq(1, dim(terrestrial_specific_trait_table_2)[1], 1)+0.4),
                                            label= paste("italic(k)==", c(terrestrial_specific_trait_table_2["Tail Length", "K"],
                                                                          terrestrial_specific_trait_table_2["Reproductive Rate", "K"],
                                                                          terrestrial_specific_trait_table_2["PO Activity", "K"],
                                                                          terrestrial_specific_trait_table_2["Metabolic Rate", "K"],
                                                                          terrestrial_specific_trait_table_2["Mass", "K"],
                                                                          terrestrial_specific_trait_table_2["Longevity", "K"]), "~","(", 
                                                                        c(terrestrial_specific_trait_table_2["Tail Length", "group_no"],
                                                                          terrestrial_specific_trait_table_2["Reproductive Rate", "group_no"],
                                                                          terrestrial_specific_trait_table_2["PO Activity", "group_no"], 
                                                                          terrestrial_specific_trait_table_2["Metabolic Rate", "group_no"],
                                                                          terrestrial_specific_trait_table_2["Mass", "group_no"],
                                                                          terrestrial_specific_trait_table_2["Longevity", "group_no"]), 
                                                         ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                            geom_label(aes(label=c(paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Tail Length", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                   paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Reproductive Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["PO Activity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Terrestrial_Specific_Trait_Model_CVR_Estimates["Longevity", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                       x = -0.75, y = (seq(1, dim(terrestrial_specific_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_terrestrial_specific_trait_CVR_2 #(400x540)

##### Summary of Terrestrial Plots #####

Terrestrial_Layout <- rbind(c(1, 2), 
                            c(1, 2),
                            c(1, 2),
                            c(1, 2),
                            c(1, 2),
                            c(1, 2),
                            c(1, 3),
                            c(1, 3),
                            c(1, 3),
                            c(1, 3),
                            c(1, 4),
                            c(1, 4))

Terrestrial_Combined_CVR <- grid.arrange(density_terrestrial_specific_trait_CVR, density_terrestrial_trait_CVR, 
                                         density_terrestrial_fluctuation_CVR, density_terrestrial_plasticity_CVR, 
                                         layout_matrix = Terrestrial_Layout)

Terrestrial_Combined_CVR #(850 x 900 - does not include amplitude plot)

##### Acclimation Subset Model - CVR #####
Acclimation_Subset_Data <- Individual_Subset_Data %>% filter(Plasticity_Mechanism == "Acclimation")
Acclimation_Species <- Acclimation_Subset_Data %>% select("phylo") %>% unique()

Acclimation_A_cor <- as.data.frame(A_cor)
Acclimation_A_cor <- Acclimation_A_cor[c(Acclimation_Species$phylo), c(Acclimation_Species$phylo)]
Acclimation_A_cor <- as.matrix(Acclimation_A_cor)

Acclimation_VCV_InCVR <- make_VCV_matrix(Acclimation_Subset_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Model_CVR <- metafor::rma.mv(InCVR ~ 1, V = Acclimation_VCV_InCVR, test = "t", dfs = "contain",
                                             random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                           ~1|Shared_Animal_Number, ~1|Measurement), 
                                             R = list(phylo=Acclimation_A_cor), data = Acclimation_Subset_Data, method = "REML", sparse = TRUE, 
                                             control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Acclimation_Model_CVR.rds")
  } else {
            Acclimation_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Acclimation_Model_CVR.rds")})

Acclimation_Model_CVR_rob <- robust(Acclimation_Model_CVR, cluster = Acclimation_Subset_Data$Study_ID, adjust = TRUE)

Acclimation_Model_CVR_Estimates <- data.frame(estimate = Acclimation_Model_CVR$b, 
                                              ci.lb = Acclimation_Model_CVR$ci.lb, 
                                              ci.ub = Acclimation_Model_CVR$ci.ub)
Acclimation_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Model_CVR), 2))

#### Acclimation Subset Model - Fluctuation Amplitude Meta-Regression - CVR ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Amplitude_Model_CVR <- metafor::rma.mv(InCVR, V = Acclimation_VCV_InCVR, test = "t", dfs = "contain",
                                                       mods = ~ T2_Magnitude - 1,
                                                       random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                     ~1|Shared_Animal_Number, ~1|Measurement), 
                                                       R = list(phylo=Acclimation_A_cor), data = Acclimation_Subset_Data, method = "REML", sparse = TRUE, 
                                                       control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Amplitude_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Acclimation_Amplitude_Model_CVR.rds")
  } else {
            Acclimation_Amplitude_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Acclimation_Amplitude_Model_CVR.rds")})

Acclimation_Amplitude_Model_CVR_rob <- robust(Acclimation_Amplitude_Model_CVR, cluster = Acclimation_Subset_Data$Study_ID, adjust = TRUE)

Acclimation_Amplitude_Model_CVR_Estimates <- data.frame(estimate = Acclimation_Amplitude_Model_CVR$b, 
                                                        ci.lb = Acclimation_Amplitude_Model_CVR$ci.lb, 
                                                        ci.ub = Acclimation_Amplitude_Model_CVR$ci.ub)
Acclimation_Amplitude_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Amplitude_Model_CVR), 2))

# Graph Preparing

Acclimation_Plot_Data <- Acclimation_Subset_Data
Acclimation_Plot_Data <- Acclimation_Plot_Data %>% mutate(n_category = ifelse(n_R1.1 <= 25, "25", 
                                                                       ifelse(n_R1.1 > 25 & n_R1.1 <= 50, "50", 
                                                                       ifelse(n_R1.1 > 50 & n_R1.1 <= 75, "75", "> 75"))))

# Graph Code

Acclimation_Amplitude_Plot_CVR <- ggplot(Acclimation_Plot_Data, aes(x = T2_Magnitude, y = InCVR)) + 
                                  geom_point(aes(x = T2_Magnitude, y = InCVR, 
                                             size = fct_relevel(n_category, c("25", "50", "75", "> 75"))), 
                                             shape = 21, fill = "#4292c6", alpha = 0.5) + 
                                  labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (lnCVR)", 
                                       size = "Sample Size", title = "Acclimation Treatments") +
                                  theme_bw() +
                                  theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                                  theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                                  theme(legend.position = "bottom", legend.direction = "horizontal") + 
                                  geom_hline(yintercept = Acclimation_Model_CVR_Estimates$estimate, lty = 2) + 
                                  geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                                  stat_poly_eq(formula = y ~ x, 
                                            aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                            parse = TRUE) +
                                  coord_cartesian(xlim = c(0, 30), 
                                                  ylim = c(-2.5, 2.5))

Acclimation_Amplitude_Plot_CVR #(400x400)

#### Acclimation Subset Model - Exposure Time Meta-Regression - CVR ####
run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Exposure_Model_CVR <- metafor::rma.mv(InCVR, V = Acclimation_VCV_InCVR, test = "t", dfs = "contain",
                                                      mods = ~ Acclimation_Exposure_Time - 1,
                                                      random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                    ~1|Shared_Animal_Number, ~1|Measurement), 
                                                      R = list(phylo=Acclimation_A_cor), data = Acclimation_Subset_Data, method = "REML", sparse = TRUE, 
                                                      control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Exposure_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Acclimation_Exposure_Model_CVR.rds")
  } else {
            Acclimation_Exposure_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Acclimation_Exposure_Model_CVR.rds")})

Acclimation_Exposure_Model_CVR_rob <- robust(Acclimation_Exposure_Model_CVR, cluster = Acclimation_Subset_Data$Study_ID, adjust = TRUE)

Acclimation_Exposure_Model_CVR_Estimates <- data.frame(estimate = Acclimation_Exposure_Model_CVR$b, 
                                                       ci.lb = Acclimation_Exposure_Model_CVR$ci.lb, 
                                                       ci.ub = Acclimation_Exposure_Model_CVR$ci.ub)
Acclimation_Exposure_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Exposure_Model_CVR), 2))

# Graph Preparing

Acclimation_Exposure_Plot_Data <- Acclimation_Subset_Data
Acclimation_Exposure_Plot_Data <- Acclimation_Exposure_Plot_Data %>% mutate(n_category = ifelse(n_R1.1 <= 25, "25", 
                                                                                         ifelse(n_R1.1 > 25 & n_R1.1 <= 50, "50", 
                                                                                         ifelse(n_R1.1 > 50 & n_R1.1 <= 75, "75", "> 75"))))

# Graph Code

Acclimation_Exposure_Amplitude_Plot_CVR <- ggplot(Acclimation_Exposure_Plot_Data, aes(x = Acclimation_Exposure_Time, y = InCVR)) + 
                                           geom_point(aes(x = Acclimation_Exposure_Time, y = InCVR, 
                                                      size = fct_relevel(n_category, c("25", "50", "75", "> 75"))), 
                                                      shape = 21, fill = "#4292c6", alpha = 0.5) + 
                                           labs(x = "Exposure Time (Days)", y = "Effect Size (lnCVR)", 
                                                size = "Sample Size", title = "Acclimation Treatments") +
                                           theme_bw() +
                                           theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                                           theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                                           theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                                           theme(legend.position = "bottom", legend.direction = "horizontal") + 
                                           geom_hline(yintercept = Acclimation_Model_CVR_Estimates$estimate, lty = 2) + 
                                           geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                                           stat_poly_eq(formula = y ~ x, 
                                                        aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                                        parse = TRUE) +
                                           coord_cartesian(xlim = c(0, 100), 
                                                           ylim = c(-2.5, 2.5))

Acclimation_Exposure_Amplitude_Plot_CVR #(400x400)

#### Acclimation Subset Model - Fluctuation Frequency Meta-Regression - CVR ####
Acclimation_Frequency_Data <- Acclimation_Subset_Data %>% filter(!is.na(Number_Of_Fluctuations))
Acclimation_Frequency_Species <- Acclimation_Frequency_Data %>% select("phylo") %>% unique()

Acclimation_Frequency_A_cor <- as.data.frame(A_cor)
Acclimation_Frequency_A_cor <- Acclimation_Frequency_A_cor[c(Acclimation_Frequency_Species$phylo), c(Acclimation_Frequency_Species$phylo)]
Acclimation_Frequency_A_cor <- as.matrix(Acclimation_Frequency_A_cor)

Acclimation_Frequency_VCV_InCVR <- make_VCV_matrix(Acclimation_Frequency_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Frequency_Model_CVR <- metafor::rma.mv(InCVR, V = Acclimation_Frequency_VCV_InCVR, test = "t", dfs = "contain",
                                                       mods = ~ Number_Of_Fluctuations - 1,
                                                       random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                     ~1|Shared_Animal_Number, ~1|Measurement), 
                                                       R = list(phylo=Acclimation_Frequency_A_cor), data = Acclimation_Frequency_Data, method = "REML", sparse = TRUE, 
                                                       control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Frequency_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Acclimation_Frequency_Model_CVR.rds")
  } else {
            Acclimation_Frequency_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Acclimation_Frequency_Model_CVR.rds")})

Acclimation_Frequency_Model_CVR_rob <- robust(Acclimation_Frequency_Model_CVR, cluster = Acclimation_Frequency_Data$Study_ID, adjust = TRUE)

Acclimation_Frequency_Model_CVR_Estimates <- data.frame(estimate = Acclimation_Frequency_Model_CVR$b, 
                                                        ci.lb = Acclimation_Frequency_Model_CVR$ci.lb, 
                                                        ci.ub = Acclimation_Frequency_Model_CVR$ci.ub)
Acclimation_Frequency_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Frequency_Model_CVR), 2))

# Graph Preparing

Acclimation_Frequency_Plot_Data <- Acclimation_Frequency_Data
Acclimation_Frequency_Plot_Data <- Acclimation_Frequency_Plot_Data %>% mutate(n_category = ifelse(n_R1.1 <= 25, "25", 
                                                                                           ifelse(n_R1.1 > 25 & n_R1.1 <= 50, "50", 
                                                                                           ifelse(n_R1.1 > 50 & n_R1.1 <= 75, "75", "> 75"))))

# Graph Code

Acclimation_Frequency_Plot_CVR <- ggplot(Acclimation_Frequency_Plot_Data, aes(x = Number_Of_Fluctuations, y = InCVR)) + 
                                  geom_point(aes(x = Number_Of_Fluctuations, y = InCVR, 
                                             size = fct_relevel(n_category, c("25", "50", "75", "> 75"))), 
                                             shape = 21, fill = "#4292c6", alpha = 0.5) + 
                                  labs(x = "Number of Fluctuations", y = "Effect Size (lnCVR)", 
                                       size = "Sample Size", title = "Acclimation Treatments") +
                                  theme_bw() +
                                  theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                                  theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                                  theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                                  theme(legend.position = "bottom", legend.direction = "horizontal") + 
                                  geom_hline(yintercept = Acclimation_Model_CVR_Estimates$estimate, lty = 2) + 
                                  geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                                  stat_poly_eq(formula = y ~ x, 
                                               aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                               parse = TRUE) +
                                  coord_cartesian(xlim = c(0, 100), 
                                                  ylim = c(-2.5, 2.5))

Acclimation_Frequency_Plot_CVR #(400x400)

#### Acclimation Subset Model - Type of Fluctuation Meta-Regression - CVR ####
Acclimation_Fluctuation_Data <- Acclimation_Subset_Data %>% filter(!is.na(Fluctuation_Category))

Acclimation_Fluctuation_Exploration <- Acclimation_Fluctuation_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Acclimation_Fluctuation_Exploration) <- Acclimation_Fluctuation_Exploration$Fluctuation_Category

Acclimation_Fluctuation_Data <- Acclimation_Fluctuation_Data %>% filter(Fluctuation_Category != "Stochastic")

Acclimation_Fluctuation_Species_Count <- Acclimation_Fluctuation_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                         filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Acclimation_Fluctuation_Species_Count) <- Acclimation_Fluctuation_Species_Count$Fluctuation_Category

Acclimation_Fluctuation_Study_Count <- Acclimation_Fluctuation_Data %>% select("Study_ID", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                       filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Acclimation_Fluctuation_Study_Count) <- Acclimation_Fluctuation_Study_Count$Fluctuation_Category

Acclimation_Fluctuation_Species <- Acclimation_Fluctuation_Data %>% select("phylo") %>% unique()

Acclimation_Fluctuation_A_cor <- as.data.frame(A_cor)
Acclimation_Fluctuation_A_cor <- Acclimation_Fluctuation_A_cor[c(Acclimation_Fluctuation_Species$phylo), c(Acclimation_Fluctuation_Species$phylo)]
Acclimation_Fluctuation_A_cor <- as.matrix(Acclimation_Fluctuation_A_cor)

Acclimation_Fluctuation_VCV_InCVR <- make_VCV_matrix(Acclimation_Fluctuation_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Fluctuation_Model_CVR <- metafor::rma.mv(InCVR, V = Acclimation_Fluctuation_VCV_InCVR, test = "t", dfs = "contain",
                                                         mods = ~ Fluctuation_Category - 1,
                                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                       ~1|Shared_Animal_Number, ~1|Measurement), 
                                                         R = list(phylo=Acclimation_Fluctuation_A_cor), data = Acclimation_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                         control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Fluctuation_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Acclimation_Fluctuation_Model_CVR.rds")
  } else {
            Acclimation_Fluctuation_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Acclimation_Fluctuation_Model_CVR.rds")})

Acclimation_Fluctuation_Model_CVR_rob <- robust(Acclimation_Fluctuation_Model_CVR, cluster = Acclimation_Fluctuation_Data$Study_ID, adjust = TRUE)

Acclimation_Fluctuation_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Acclimation_Fluctuation_Model_CVR$b), 21, 100),
                                                          estimate = Acclimation_Fluctuation_Model_CVR$b, 
                                                          ci.lb = Acclimation_Fluctuation_Model_CVR$ci.lb, 
                                                          ci.ub = Acclimation_Fluctuation_Model_CVR$ci.ub)
rownames(Acclimation_Fluctuation_Model_CVR_Estimates) <- Acclimation_Fluctuation_Model_CVR_Estimates$Category
Acclimation_Fluctuation_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Fluctuation_Model_CVR), 2))

# Preparing Graph - Combined

acclimation_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise")

acclimation_fluctuation_k <- data.frame("k" = c(Acclimation_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                Acclimation_Fluctuation_Exploration["Alternating", "Freq"], 
                                                Acclimation_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                        row.names = acclimation_fluctuation_rnames)

acclimation_fluctuation_group_no <- data.frame("Spp No." = c(Acclimation_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                             Acclimation_Fluctuation_Species_Count["Alternating", "Freq"], 
                                                             Acclimation_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                               row.names = acclimation_fluctuation_rnames)

acclimation_fluctuation_study <- data.frame("Study" = c(Acclimation_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                        Acclimation_Fluctuation_Study_Count["Alternating", "Freq"], 
                                                        Acclimation_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                               row.names = acclimation_fluctuation_rnames)

Acclimation_Fluctuation_Model_CVR_Estimates_Reorder <- Acclimation_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise"), ]

acclimation_fluctuation_table <- data.frame(estimate = Acclimation_Fluctuation_Model_CVR_Estimates_Reorder[,"estimate"], 
                                            lowerCL = Acclimation_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.lb"], 
                                            upperCL = Acclimation_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.ub"], 
                                            K = acclimation_fluctuation_k[,1], 
                                            group_no = acclimation_fluctuation_group_no[,1], 
                                            row.names = acclimation_fluctuation_rnames)
acclimation_fluctuation_table$name <- row.names(acclimation_fluctuation_table)

acclimation_fluctuation_raw_mean <- c(unlist(unname(Acclimation_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                      select("InCVR"))), 
                                      unlist(unname(Acclimation_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                      select("InCVR"))), 
                                      unlist(unname(Acclimation_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                      select("InCVR"))))

acclimation_fluctuation_raw_name <- c(replicate(151, "Sinusoidal (Sine Curve)"), 
                                      replicate(75, "Alternating"), 
                                      replicate(76, "Stepwise"))

acclimation_fluctuation_raw_df <- data.frame("Model" = acclimation_fluctuation_raw_name, 
                                             "Effect" = acclimation_fluctuation_raw_mean)

# Graph code - Combined

Acclimation_Fluctuation_Order <- c("Stepwise", "Alternating", "Sinusoidal (Sine Curve)")

density_acclimation_fluctuation_CVR <- acclimation_fluctuation_table %>% mutate(name = fct_relevel(name, Acclimation_Fluctuation_Order)) %>%
                                       ggplot() +
                                       geom_density_ridges(data = acclimation_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Acclimation_Fluctuation_Order)), 
                                                           aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                       geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                      size = 1) +
                                       geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table)[1], 1)), xmin = min(acclimation_fluctuation_raw_df$Effect)-0.13, xmax = -1.5, colour = name),
                                                      size = 1) +
                                       geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                       size = 1, fatten = 2) +
                                       theme_bw() +
                                       guides(fill = "none", colour = "none") +
                                       labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                       theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                        vjust = c(-2.7, -2.7, -0.8))) +
                                       theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                       theme(axis.ticks = element_blank()) +
                                       theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                       scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                       scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                       scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#2B4E7A")) +
                                       coord_cartesian(xlim = c(-1, 1)) +
                                       annotate('text',  x = 1, y = (seq(1, dim(acclimation_fluctuation_table)[1], 1)+0.4),
                                       label= paste("italic(k)==", c(acclimation_fluctuation_table["Stepwise", "K"], 
                                                                     acclimation_fluctuation_table["Alternating", "K"], 
                                                                     acclimation_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                   c(acclimation_fluctuation_table["Stepwise", "group_no"], 
                                                                     acclimation_fluctuation_table["Alternating", "group_no"], 
                                                                     acclimation_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                                    ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                       geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                              paste(format(round(mean(exp(Acclimation_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                              paste(format(round(mean(exp(Acclimation_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                  x = -0.75, y = (seq(1, dim(acclimation_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_acclimation_fluctuation_CVR #(400x320)

# Preparing Graph - Part 1

acclimation_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)", "Alternating")

acclimation_fluctuation_k_1 <- data.frame("k" = c(Acclimation_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                  Acclimation_Fluctuation_Exploration["Alternating", "Freq"]), 
                                          row.names = acclimation_fluctuation_rnames_1)

acclimation_fluctuation_group_no_1 <- data.frame("Spp No." = c(Acclimation_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                               Acclimation_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                                 row.names = acclimation_fluctuation_rnames_1)

acclimation_fluctuation_study_1 <- data.frame("Study" = c(Acclimation_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                          Acclimation_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                              row.names = acclimation_fluctuation_rnames_1)

Acclimation_Fluctuation_Model_CVR_Estimates_Reorder_1 <- Acclimation_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

acclimation_fluctuation_table_1 <- data.frame(estimate = Acclimation_Fluctuation_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                              lowerCL = Acclimation_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                              upperCL = Acclimation_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                              K = acclimation_fluctuation_k_1[,1], 
                                              group_no = acclimation_fluctuation_group_no_1[,1], 
                                              row.names = acclimation_fluctuation_rnames_1)
acclimation_fluctuation_table_1$name <- row.names(acclimation_fluctuation_table_1)

acclimation_fluctuation_raw_mean_1 <- c(unlist(unname(Acclimation_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                      select("InCVR"))), 
                                        unlist(unname(Acclimation_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                      select("InCVR"))))

acclimation_fluctuation_raw_name_1 <- c(replicate(151, "Sinusoidal (Sine Curve)"), 
                                        replicate(75, "Alternating"))

acclimation_fluctuation_raw_df_1 <- data.frame("Model" = acclimation_fluctuation_raw_name_1, 
                                               "Effect" = acclimation_fluctuation_raw_mean_1)

# Graph code - Part 1

Acclimation_Fluctuation_Order_1 <- c("Alternating", "Sinusoidal (Sine Curve)")

density_acclimation_fluctuation_CVR_1 <- acclimation_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Acclimation_Fluctuation_Order_1)) %>%
                                         ggplot() +
                                         geom_density_ridges(data = acclimation_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Acclimation_Fluctuation_Order_1)), 
                                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                 scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                         geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                        size = 1) +
                                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                         size = 1, fatten = 2) +
                                         theme_bw() +
                                         guides(fill = "none", colour = "none") +
                                         labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                         theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                               vjust = c(-2.7, -0.8))) +
                                         theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                         theme(axis.ticks = element_blank()) +
                                         theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                         scale_colour_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                         scale_fill_manual(values = c("#4A6E9C", "#2B4E7A")) +
                                         coord_cartesian(xlim = c(-1, 1)) +
                                         annotate('text',  x = 1, y = (seq(1, dim(acclimation_fluctuation_table_1)[1], 1)+0.4),
                                         label= paste("italic(k)==", c(acclimation_fluctuation_table_1["Alternating", "K"], 
                                                                       acclimation_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                     c(acclimation_fluctuation_table_1["Alternating", "group_no"], 
                                                                       acclimation_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                                      ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                         geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                paste(format(round(mean(exp(Acclimation_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                    x = -0.75, y = (seq(1, dim(acclimation_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_acclimation_fluctuation_CVR_1 #(400x240)

# Preparing Graph - Part 2

acclimation_fluctuation_rnames_2 <- c("Stepwise")

acclimation_fluctuation_k_2 <- data.frame("k" = c(Acclimation_Fluctuation_Exploration["Stepwise", "Freq"]), 
                                          row.names = acclimation_fluctuation_rnames_2)

acclimation_fluctuation_group_no_2 <- data.frame("Spp No." = c(Acclimation_Fluctuation_Species_Count["Stepwise", "Freq"]), 
                                                 row.names = acclimation_fluctuation_rnames_2)

acclimation_fluctuation_study_2 <- data.frame("Study" = c(Acclimation_Fluctuation_Study_Count["Stepwise", "Freq"]), 
                                              row.names = acclimation_fluctuation_rnames_2)

Acclimation_Fluctuation_Model_CVR_Estimates_Reorder_2 <- Acclimation_Fluctuation_Model_CVR_Estimates[c("Stepwise"), ]

acclimation_fluctuation_table_2 <- data.frame(estimate = Acclimation_Fluctuation_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                              lowerCL = Acclimation_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                              upperCL = Acclimation_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                              K = acclimation_fluctuation_k_2[,1], 
                                              group_no = acclimation_fluctuation_group_no_2[,1], 
                                              row.names = acclimation_fluctuation_rnames_2)
acclimation_fluctuation_table_2$name <- row.names(acclimation_fluctuation_table_2)

acclimation_fluctuation_raw_mean_2 <- c(unlist(unname(Acclimation_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                      select("InCVR"))))

acclimation_fluctuation_raw_name_2 <- c(replicate(76, "Stepwise"))

acclimation_fluctuation_raw_df_2 <- data.frame("Model" = acclimation_fluctuation_raw_name_2, 
                                               "Effect" = acclimation_fluctuation_raw_mean_2)

# Graph code - Part 2

Acclimation_Fluctuation_Order_2 <- c("Stepwise")

density_acclimation_fluctuation_CVR_2 <- acclimation_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Acclimation_Fluctuation_Order_2)) %>%
                                         ggplot() +
                                         geom_density_ridges(data = acclimation_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Acclimation_Fluctuation_Order_2)), 
                                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                 scale = 0.07, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                         geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                        size = 1) +
                                         geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table_2)[1], 1)), xmin = max(acclimation_fluctuation_raw_df_2$Effect)+0.03, xmax = 1.5, colour = name),
                                                        size = 1) +
                                         geom_linerange(aes(y = rev(seq(1, dim(acclimation_fluctuation_table_2)[1], 1)), xmin = min(acclimation_fluctuation_raw_df_2$Effect)-0.03, xmax = -1.5, colour = name),
                                                        size = 1) +
                                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                         size = 1, fatten = 2) +
                                         theme_bw() +
                                         guides(fill = "none", colour = "none") +
                                         labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                         theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                               vjust = c(-2.7))) +
                                         theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                         theme(axis.ticks = element_blank()) +
                                         theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                         scale_colour_manual(values = c("#5D7AA1")) +
                                         scale_fill_manual(values = c("#5D7AA1")) +
                                         coord_cartesian(xlim = c(-1, 1)) +
                                         annotate('text',  x = 1, y = (seq(1, dim(acclimation_fluctuation_table_2)[1], 1)+0.4),
                                         label= paste("italic(k)==", c(acclimation_fluctuation_table_2["Stepwise", "K"]), "~","(", 
                                                                     c(acclimation_fluctuation_table_2["Stepwise", "group_no"]), 
                                                      ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                         geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                    x = -0.75, y = (seq(1, dim(acclimation_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_acclimation_fluctuation_CVR_2 #(400x160)

##### Acclimation Subset Model - Trait Meta-Regression - CVR #####
Acclimation_Trait_Exploration <- Acclimation_Subset_Data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Acclimation_Trait_Exploration) <- Acclimation_Trait_Exploration$Trait_Category

Acclimation_Trait_Data <- Acclimation_Subset_Data %>% filter(Trait_Category != "Gene Expression")

Acclimation_Trait_Species_Count <- Acclimation_Trait_Data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>%
                                   filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Acclimation_Trait_Species_Count) <- Acclimation_Trait_Species_Count$Trait_Category

Acclimation_Trait_Study_Count <- Acclimation_Trait_Data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>%
                                 filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Acclimation_Trait_Study_Count) <- Acclimation_Trait_Study_Count$Trait_Category

Acclimation_Trait_Species <- Acclimation_Trait_Data %>% select("phylo") %>% unique()

Acclimation_Trait_A_cor <- as.data.frame(A_cor)
Acclimation_Trait_A_cor <- Acclimation_Trait_A_cor[c(Acclimation_Trait_Species$phylo), c(Acclimation_Trait_Species$phylo)]
Acclimation_Trait_A_cor <- as.matrix(Acclimation_Trait_A_cor)

Acclimation_Trait_VCV_InCVR <- make_VCV_matrix(Acclimation_Trait_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Trait_Model_CVR <- metafor::rma.mv(InCVR, V = Acclimation_Trait_VCV_InCVR, test = "t", dfs = "contain",
                                                   mods = ~ Trait_Category - 1,
                                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                 ~1|Shared_Animal_Number, ~1|Measurement), 
                                                   R = list(phylo=Acclimation_Trait_A_cor), data = Acclimation_Trait_Data, method = "REML", sparse = TRUE, 
                                                   control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Trait_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Acclimation_Trait_Model_CVR.rds")
  } else {
            Acclimation_Trait_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Acclimation_Trait_Model_CVR.rds")})

Acclimation_Trait_Model_CVR_rob <- robust(Acclimation_Trait_Model_CVR, cluster = Acclimation_Trait_Data$Study_ID, adjust = TRUE)

Acclimation_Trait_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Acclimation_Trait_Model_CVR$b), 15, 100),
                                                    estimate = Acclimation_Trait_Model_CVR$b, 
                                                    ci.lb = Acclimation_Trait_Model_CVR$ci.lb, 
                                                    ci.ub = Acclimation_Trait_Model_CVR$ci.ub)
rownames(Acclimation_Trait_Model_CVR_Estimates) <- Acclimation_Trait_Model_CVR_Estimates$Category
Acclimation_Trait_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Trait_Model_CVR), 2))

# Preparing Graph - Combined

acclimation_trait_rnames <- c("Behavioural", "Biochemical Assay", "Life-history Traits", "Physiological")

acclimation_trait_k <- data.frame("k" = c(Acclimation_Trait_Exploration["Behavioural", "Freq"], 
                                          Acclimation_Trait_Exploration["Biochemical Assay", "Freq"],
                                          Acclimation_Trait_Exploration["Life-History Traits", "Freq"],
                                          Acclimation_Trait_Exploration["Physiological", "Freq"]), 
                                  row.names = acclimation_trait_rnames)

acclimation_trait_group_no <- data.frame("Spp No." = c(Acclimation_Trait_Species_Count["Behavioural", "Freq"], 
                                                       Acclimation_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                       Acclimation_Trait_Species_Count["Life-History Traits", "Freq"],
                                                       Acclimation_Trait_Species_Count["Physiological", "Freq"]), 
                                         row.names = acclimation_trait_rnames)

acclimation_trait_study <- data.frame("Study" = c(Acclimation_Trait_Study_Count["Behavioural", "Freq"], 
                                                  Acclimation_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                  Acclimation_Trait_Study_Count["Life-History Traits", "Freq"],
                                                  Acclimation_Trait_Study_Count["Physiological", "Freq"]), 
                                         row.names = acclimation_trait_rnames)

acclimation_trait_table <- data.frame(estimate = Acclimation_Trait_Model_CVR_Estimates[,"estimate"], 
                                      lowerCL = Acclimation_Trait_Model_CVR_Estimates[,"ci.lb"], 
                                      upperCL = Acclimation_Trait_Model_CVR_Estimates[,"ci.ub"], 
                                      K = acclimation_trait_k[,1], 
                                      group_no = acclimation_trait_group_no[,1], 
                                      row.names = acclimation_trait_rnames)
acclimation_trait_table$name <- row.names(acclimation_trait_table)

acclimation_trait_raw_mean <- c(unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Behavioural") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                                select("InCVR"))))

acclimation_trait_raw_name <- c(replicate(24, "Behavioural"), 
                                replicate(140, "Biochemical Assay"), 
                                replicate(28, "Life-history Traits"), 
                                replicate(134, "Physiological"))

acclimation_trait_raw_df <- data.frame("Model" = acclimation_trait_raw_name, 
                                       "Effect" = acclimation_trait_raw_mean)

# Graph code - Combined

Acclimation_Trait_Order <- c("Physiological", "Life-history Traits", "Biochemical Assay", "Behavioural")

density_acclimation_trait_CVR <- acclimation_trait_table %>% mutate(name = fct_relevel(name, Acclimation_Trait_Order)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = acclimation_trait_raw_df %>% mutate(Model = fct_relevel(Model, Acclimation_Trait_Order)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                 size = 1, fatten = 2) +
                                 theme_bw() +
                                 guides(fill = "none", colour = "none") +
                                 labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                 theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                  vjust = c(-2.7, -0.8, -0.8, -2.7))) +
                                 theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                 theme(axis.ticks = element_blank()) +
                                 theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                 theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                 scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                 scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                 scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                 coord_cartesian(xlim = c(-1, 1)) +
                                 annotate('text',  x = 1, y = (seq(1, dim(acclimation_trait_table)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(acclimation_trait_table["Physiological", "K"],
                                                               acclimation_trait_table["Life-history Traits", "K"],
                                                               acclimation_trait_table["Biochemical Assay", "K"],
                                                               acclimation_trait_table["Behavioural", "K"]), "~","(", 
                                                             c(acclimation_trait_table["Physiological", "group_no"],
                                                               acclimation_trait_table["Life-history Traits", "group_no"],
                                                               acclimation_trait_table["Biochemical Assay", "group_no"],
                                                               acclimation_trait_table["Behavioural", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                 geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Trait_Model_CVR_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                        paste(format(round(mean(exp(Acclimation_Trait_Model_CVR_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                        paste(format(round(mean(exp(Acclimation_Trait_Model_CVR_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                        paste(format(round(mean(exp(Acclimation_Trait_Model_CVR_Estimates["Behavioural", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.75, y = (seq(1, dim(acclimation_trait_table)[1], 1)+0.4)), size = 3.5)

density_acclimation_trait_CVR #(400x400)

# Preparing Graph - Part 1

acclimation_trait_rnames_1 <- c("Behavioural", "Biochemical Assay")

acclimation_trait_k_1 <- data.frame("k" = c(Acclimation_Trait_Exploration["Behavioural", "Freq"], 
                                            Acclimation_Trait_Exploration["Biochemical Assay", "Freq"]), 
                                    row.names = acclimation_trait_rnames_1)

acclimation_trait_group_no_1 <- data.frame("Spp No." = c(Acclimation_Trait_Species_Count["Behavioural", "Freq"], 
                                                         Acclimation_Trait_Species_Count["Biochemical Assay", "Freq"]), 
                                           row.names = acclimation_trait_rnames_1)

acclimation_trait_study_1 <- data.frame("Study" = c(Acclimation_Trait_Study_Count["Behavioural", "Freq"], 
                                                    Acclimation_Trait_Study_Count["Biochemical Assay", "Freq"]), 
                                        row.names = acclimation_trait_rnames_1)

Acclimation_Trait_Model_CVR_Estimates_Reorder_1 <- Acclimation_Trait_Model_CVR_Estimates[c("Behavioural", "Biochemical Assay"), ]

acclimation_trait_table_1 <- data.frame(estimate = Acclimation_Trait_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                        lowerCL = Acclimation_Trait_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                        upperCL = Acclimation_Trait_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                        K = acclimation_trait_k_1[,1], 
                                        group_no = acclimation_trait_group_no_1[,1], 
                                        row.names = acclimation_trait_rnames_1)
acclimation_trait_table_1$name <- row.names(acclimation_trait_table_1)

acclimation_trait_raw_mean_1 <- c(unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Behavioural") %>% 
                                                select("InCVR"))), 
                                  unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                                select("InCVR"))))

acclimation_trait_raw_name_1 <- c(replicate(24, "Behavioural"), 
                                  replicate(140, "Biochemical Assay"))

acclimation_trait_raw_df_1 <- data.frame("Model" = acclimation_trait_raw_name_1, 
                                         "Effect" = acclimation_trait_raw_mean_1)

# Graph code - Part 1

Acclimation_Trait_Order_1 <- c("Biochemical Assay", "Behavioural")

density_acclimation_trait_CVR_1 <- acclimation_trait_table_1 %>% mutate(name = fct_relevel(name, Acclimation_Trait_Order_1)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = acclimation_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Acclimation_Trait_Order_1)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table_1)[1], 1)), xmin = min(acclimation_trait_raw_df_1$Effect)-0.07, xmax = -1.5, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                         vjust = c(-0.8, -2.7))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                   scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(acclimation_trait_table_1)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(acclimation_trait_table["Biochemical Assay", "K"],
                                                                 acclimation_trait_table["Behavioural", "K"]), "~","(", 
                                                               c(acclimation_trait_table["Biochemical Assay", "group_no"],
                                                                 acclimation_trait_table["Behavioural", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Trait_Model_CVR_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Acclimation_Trait_Model_CVR_Estimates["Behavioural", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(acclimation_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_acclimation_trait_CVR_1 #(400x240)

# Preparing Graph - Part 2

acclimation_trait_rnames_2 <- c("Life-history Traits", "Physiological")

acclimation_trait_k_2 <- data.frame("k" = c(Acclimation_Trait_Exploration["Life-History Traits", "Freq"],
                                            Acclimation_Trait_Exploration["Physiological", "Freq"]), 
                                    row.names = acclimation_trait_rnames_2)

acclimation_trait_group_no_2 <- data.frame("Spp No." = c(Acclimation_Trait_Species_Count["Life-History Traits", "Freq"],
                                                         Acclimation_Trait_Species_Count["Physiological", "Freq"]), 
                                           row.names = acclimation_trait_rnames_2)

acclimation_trait_study_2 <- data.frame("Study" = c(Acclimation_Trait_Study_Count["Life-History Traits", "Freq"],
                                                    Acclimation_Trait_Study_Count["Physiological", "Freq"]), 
                                        row.names = acclimation_trait_rnames_2)

Acclimation_Trait_Model_CVR_Estimates_Reorder_2 <- Acclimation_Trait_Model_CVR_Estimates[c("Life-History Traits", "Physiological"), ]

acclimation_trait_table_2 <- data.frame(estimate = Acclimation_Trait_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                        lowerCL = Acclimation_Trait_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                        upperCL = Acclimation_Trait_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                        K = acclimation_trait_k_2[,1], 
                                        group_no = acclimation_trait_group_no_2[,1], 
                                        row.names = acclimation_trait_rnames_2)
acclimation_trait_table_2$name <- row.names(acclimation_trait_table_2)

acclimation_trait_raw_mean_2 <- c(unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                                select("InCVR"))), 
                                  unlist(unname(Acclimation_Trait_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                                select("InCVR"))))

acclimation_trait_raw_name_2 <- c(replicate(28, "Life-history Traits"), 
                                  replicate(134, "Physiological"))

acclimation_trait_raw_df_2 <- data.frame("Model" = acclimation_trait_raw_name_2, 
                                         "Effect" = acclimation_trait_raw_mean_2)

# Graph code - Part 2

Acclimation_Trait_Order_2 <- c("Physiological", "Life-history Traits")

density_acclimation_trait_CVR_2 <- acclimation_trait_table_2 %>% mutate(name = fct_relevel(name, Acclimation_Trait_Order_2)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = acclimation_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Acclimation_Trait_Order_2)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                         vjust = c(-2.7, -0.8))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                   scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(acclimation_trait_table_2)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(acclimation_trait_table_2["Physiological", "K"],
                                                                 acclimation_trait_table_2["Life-history Traits", "K"]), "~","(", 
                                                               c(acclimation_trait_table_2["Physiological", "group_no"],
                                                                 acclimation_trait_table_2["Life-history Traits", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Trait_Model_CVR_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Acclimation_Trait_Model_CVR_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(acclimation_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_acclimation_trait_CVR_2 #(400x240)

##### Acclimation Subset Model - Life-History Stage Meta-Regression - CVR #####
Acclimation_Stage_Exploration <- Acclimation_Subset_Data %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_Stage_Exploration) <- Acclimation_Stage_Exploration$Acclimation_Life.History_Stage_Category

Acclimation_Stage_Data <- Acclimation_Subset_Data %>% filter(Acclimation_Life.History_Stage_Category != "Pupae")

Acclimation_Stage_Species_Count <- Acclimation_Stage_Data %>% select("Scientific_Name", "Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame() %>%
                                   filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_Stage_Species_Count) <- Acclimation_Stage_Species_Count$Acclimation_Life.History_Stage_Category

Acclimation_Stage_Study_Count <- Acclimation_Stage_Data %>% select("Study_ID", "Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame() %>%
                                 filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_Stage_Study_Count) <- Acclimation_Stage_Study_Count$Acclimation_Life.History_Stage_Category

Acclimation_Stage_Species <- Acclimation_Stage_Data %>% select("phylo") %>% unique()

Acclimation_Stage_A_cor <- as.data.frame(A_cor)
Acclimation_Stage_A_cor <- Acclimation_Stage_A_cor[c(Acclimation_Stage_Species$phylo), c(Acclimation_Stage_Species$phylo)]
Acclimation_Stage_A_cor <- as.matrix(Acclimation_Stage_A_cor)

Acclimation_Stage_VCV_InCVR <- make_VCV_matrix(Acclimation_Stage_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Stage_Model_CVR <- metafor::rma.mv(InCVR, V = Acclimation_Stage_VCV_InCVR, test = "t", dfs = "contain",
                                                   mods = ~ Acclimation_Life.History_Stage_Category - 1,
                                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                 ~1|Shared_Animal_Number, ~1|Measurement), 
                                                   R = list(phylo=Acclimation_Stage_A_cor), data = Acclimation_Stage_Data, method = "REML", sparse = TRUE, 
                                                   control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Stage_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Acclimation_Stage_Model_CVR.rds")
  } else {
            Acclimation_Stage_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Acclimation_Stage_Model_CVR.rds")})

Acclimation_Stage_Model_CVR_rob <- robust(Acclimation_Stage_Model_CVR, cluster = Acclimation_Stage_Data$Study_ID, adjust = TRUE)

Acclimation_Stage_Model_CVR_Estimates <- data.frame(Stage = substr(row.names(Acclimation_Stage_Model_CVR$b), 40, 100),
                                                    estimate = Acclimation_Stage_Model_CVR$b, 
                                                    ci.lb = Acclimation_Stage_Model_CVR$ci.lb, 
                                                    ci.ub = Acclimation_Stage_Model_CVR$ci.ub)
rownames(Acclimation_Stage_Model_CVR_Estimates) <- Acclimation_Stage_Model_CVR_Estimates$Stage
Acclimation_Stage_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Stage_Model_CVR), 2))

# Preparing Graph - Combined

acclimation_stage_rnames <- c("Adult", "Embryo", "Juvenile", "Larva")

acclimation_stage_k <- data.frame("k" = c(Acclimation_Stage_Exploration["Adult", "Freq"], 
                                          Acclimation_Stage_Exploration["Embryo", "Freq"],
                                          Acclimation_Stage_Exploration["Juvenile", "Freq"],
                                          Acclimation_Stage_Exploration["Larvae", "Freq"]), 
                                  row.names = acclimation_stage_rnames)

acclimation_stage_group_no <- data.frame("Spp No." = c(Acclimation_Stage_Species_Count["Adult", "Freq"], 
                                                       Acclimation_Stage_Species_Count["Embryo", "Freq"], 
                                                       Acclimation_Stage_Species_Count["Juvenile", "Freq"],
                                                       Acclimation_Stage_Species_Count["Larvae", "Freq"]), 
                                         row.names = acclimation_stage_rnames)

acclimation_stage_study <- data.frame("Study" = c(Acclimation_Stage_Study_Count["Adult", "Freq"], 
                                                  Acclimation_Stage_Study_Count["Embryo", "Freq"], 
                                                  Acclimation_Stage_Study_Count["Juvenile", "Freq"],
                                                  Acclimation_Stage_Study_Count["Larvae", "Freq"]), 
                                         row.names = acclimation_stage_rnames)

acclimation_stage_table <- data.frame(estimate = Acclimation_Stage_Model_CVR_Estimates[,"estimate"], 
                                      lowerCL = Acclimation_Stage_Model_CVR_Estimates[,"ci.lb"], 
                                      upperCL = Acclimation_Stage_Model_CVR_Estimates[,"ci.ub"], 
                                      K = acclimation_stage_k[,1], 
                                      group_no = acclimation_stage_group_no[,1], 
                                      row.names = acclimation_stage_rnames)
acclimation_stage_table$name <- row.names(acclimation_stage_table)

acclimation_stage_raw_mean <- c(unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Adult") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Embryo") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Juvenile") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Larvae") %>% 
                                                select("InCVR"))))

acclimation_stage_raw_name <- c(replicate(125, "Adult"), 
                                replicate(11, "Embryo"), 
                                replicate(101, "Juvenile"), 
                                replicate(92, "Larva"))

acclimation_stage_raw_df <- data.frame("Model" = acclimation_stage_raw_name, 
                                       "Effect" = acclimation_stage_raw_mean)

# Graph code - Combined

Acclimation_Stage_Order <- c("Larva", "Juvenile", "Embryo", "Adult")

density_acclimation_stage_CVR <- acclimation_stage_table %>% mutate(name = fct_relevel(name, Acclimation_Stage_Order)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = acclimation_stage_raw_df %>% mutate(Model = fct_relevel(Model, Acclimation_Stage_Order)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_stage_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                 size = 1, fatten = 2) +
                                 theme_bw() +
                                 guides(fill = "none", colour = "none") +
                                 labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                 theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                  vjust = c(-2.7, -2.7, -2.7, -2.7))) +
                                 theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                 theme(axis.ticks = element_blank()) +
                                 theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                 theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                 scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                 scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                 scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                 coord_cartesian(xlim = c(-1, 1)) +
                                 annotate('text',  x = 1, y = (seq(1, dim(acclimation_stage_table)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(acclimation_stage_table["Larva", "K"],
                                                               acclimation_stage_table["Juvenile", "K"],
                                                               acclimation_stage_table["Embryo", "K"],
                                                               acclimation_stage_table["Adult", "K"]), "~","(", 
                                                             c(acclimation_stage_table["Larva", "group_no"],
                                                               acclimation_stage_table["Juvenile", "group_no"],
                                                               acclimation_stage_table["Embryo", "group_no"],
                                                               acclimation_stage_table["Adult", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                 geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Stage_Model_CVR_Estimates["Larvae", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                        paste(format(round(mean(exp(Acclimation_Stage_Model_CVR_Estimates["Juvenile", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                        paste(format(round(mean(exp(Acclimation_Stage_Model_CVR_Estimates["Embryo", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                        paste(format(round(mean(exp(Acclimation_Stage_Model_CVR_Estimates["Adult", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.75, y = (seq(1, dim(acclimation_stage_table)[1], 1)+0.4)), size = 3.5)

density_acclimation_stage_CVR #(400x400)

# Preparing Graph - Part 1

acclimation_stage_rnames_1 <- c("Adult", "Embryo")

acclimation_stage_k_1 <- data.frame("k" = c(Acclimation_Stage_Exploration["Adult", "Freq"], 
                                            Acclimation_Stage_Exploration["Embryo", "Freq"]), 
                                    row.names = acclimation_stage_rnames_1)

acclimation_stage_group_no_1 <- data.frame("Spp No." = c(Acclimation_Stage_Species_Count["Adult", "Freq"], 
                                                         Acclimation_Stage_Species_Count["Embryo", "Freq"]), 
                                           row.names = acclimation_stage_rnames_1)

acclimation_stage_study_1 <- data.frame("Study" = c(Acclimation_Stage_Study_Count["Adult", "Freq"], 
                                                    Acclimation_Stage_Study_Count["Embryo", "Freq"]), 
                                        row.names = acclimation_stage_rnames_1)

Acclimation_Stage_Model_CVR_Estimates_Reorder_1 <- Acclimation_Stage_Model_CVR_Estimates[c("Adult", "Embryo"), ]

acclimation_stage_table_1 <- data.frame(estimate = Acclimation_Stage_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                        lowerCL = Acclimation_Stage_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                        upperCL = Acclimation_Stage_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                        K = acclimation_stage_k_1[,1], 
                                        group_no = acclimation_stage_group_no_1[,1], 
                                        row.names = acclimation_stage_rnames_1)
acclimation_stage_table_1$name <- row.names(acclimation_stage_table_1)

acclimation_stage_raw_mean_1 <- c(unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Adult") %>% 
                                                select("InCVR"))), 
                                  unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Embryo") %>% 
                                                select("InCVR"))))

acclimation_stage_raw_name_1 <- c(replicate(125, "Adult"), 
                                  replicate(11, "Embryo"))

acclimation_stage_raw_df_1 <- data.frame("Model" = acclimation_stage_raw_name_1, 
                                         "Effect" = acclimation_stage_raw_mean_1)

# Graph code - Part 1

Acclimation_Stage_Order_1 <- c("Embryo", "Adult")

density_acclimation_stage_CVR_1 <- acclimation_stage_table_1 %>% mutate(name = fct_relevel(name, Acclimation_Stage_Order_1)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = acclimation_stage_raw_df_1 %>% mutate(Model = fct_relevel(Model, Acclimation_Stage_Order_1)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_stage_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                         vjust = c(-2.7, -2.7))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                   scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(acclimation_stage_table_1)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(acclimation_stage_table_1["Embryo", "K"],
                                                                 acclimation_stage_table_1["Adult", "K"]), "~","(", 
                                                               c(acclimation_stage_table_1["Embryo", "group_no"],
                                                                 acclimation_stage_table_1["Adult", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Stage_Model_CVR_Estimates["Embryo", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Acclimation_Stage_Model_CVR_Estimates["Adult", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(acclimation_stage_table_1)[1], 1)+0.4)), size = 3.5)

density_acclimation_stage_CVR_1 #(400x240)

# Preparing Graph - Part 2

acclimation_stage_rnames_2 <- c("Juvenile", "Larva")

acclimation_stage_k_2 <- data.frame("k" = c(Acclimation_Stage_Exploration["Juvenile", "Freq"],
                                            Acclimation_Stage_Exploration["Larvae", "Freq"]), 
                                    row.names = acclimation_stage_rnames_2)

acclimation_stage_group_no_2 <- data.frame("Spp No." = c(Acclimation_Stage_Species_Count["Juvenile", "Freq"],
                                                         Acclimation_Stage_Species_Count["Larvae", "Freq"]), 
                                           row.names = acclimation_stage_rnames_2)

acclimation_stage_study_2 <- data.frame("Study" = c(Acclimation_Stage_Study_Count["Juvenile", "Freq"],
                                                    Acclimation_Stage_Study_Count["Larvae", "Freq"]), 
                                        row.names = acclimation_stage_rnames_2)

Acclimation_Stage_Model_CVR_Estimates_Reorder_2 <- Acclimation_Stage_Model_CVR_Estimates[c("Juvenile", "Larva"), ]

acclimation_stage_table_2 <- data.frame(estimate = Acclimation_Stage_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                        lowerCL = Acclimation_Stage_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                        upperCL = Acclimation_Stage_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                        K = acclimation_stage_k_2[,1], 
                                        group_no = acclimation_stage_group_no_2[,1], 
                                        row.names = acclimation_stage_rnames_2)
acclimation_stage_table_2$name <- row.names(acclimation_stage_table_2)

acclimation_stage_raw_mean_2 <- c(unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Juvenile") %>% 
                                                select("InCVR"))), 
                                  unlist(unname(Acclimation_Stage_Data %>% filter(`Acclimation_Life.History_Stage_Category` == "Larvae") %>% 
                                                select("InCVR"))))

acclimation_stage_raw_name_2 <- c(replicate(101, "Juvenile"), 
                                  replicate(92, "Larva"))

acclimation_stage_raw_df_2 <- data.frame("Model" = acclimation_stage_raw_name_2, 
                                         "Effect" = acclimation_stage_raw_mean_2)

# Graph code - Part 2

Acclimation_Stage_Order_2 <- c("Larva", "Juvenile")

density_acclimation_stage_CVR_2 <- acclimation_stage_table_2 %>% mutate(name = fct_relevel(name, Acclimation_Stage_Order_2)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = acclimation_stage_raw_df_2 %>% mutate(Model = fct_relevel(Model, Acclimation_Stage_Order_2)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_stage_table_2)[1], 1)), xmin = min(acclimation_stage_raw_df_2$Effect)-0.05, xmax = -1.5, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_stage_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                         vjust = c(-2.7, -2.7))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                   scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(acclimation_stage_table_2)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(acclimation_stage_table_2["Larva", "K"],
                                                                 acclimation_stage_table_2["Juvenile", "K"]), "~","(", 
                                                               c(acclimation_stage_table_2["Larva", "group_no"],
                                                                 acclimation_stage_table_2["Juvenile", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Stage_Model_CVR_Estimates["Larvae", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Acclimation_Stage_Model_CVR_Estimates["Juvenile", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(acclimation_stage_table_2)[1], 1)+0.4)), size = 3.5)

density_acclimation_stage_CVR_2 #(400x240)

##### Acclimation Subset Model - Class Meta-Regression - CVR #####
Acclimation_Class_Exploration <- Acclimation_Subset_Data %>% select("Class") %>% table() %>% data.frame()
rownames(Acclimation_Class_Exploration) <- Acclimation_Class_Exploration$Class

Acclimation_Class_Data <- Acclimation_Subset_Data %>% filter(Class != "Amphibia" &
                                                             Class != "Clitellata" &
                                                             Class != "Collembola")

Acclimation_Class_Species_Count <- Acclimation_Class_Data %>% select("Scientific_Name", "Class") %>% table() %>% data.frame() %>% 
                                   filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Acclimation_Class_Species_Count) <- Acclimation_Class_Species_Count$Class

Acclimation_Class_Study_Count <- Acclimation_Class_Data %>% select("Study_ID", "Class") %>% table() %>% data.frame() %>% 
                                 filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Acclimation_Class_Study_Count) <- Acclimation_Class_Study_Count$Class

Acclimation_Class_Species <- Acclimation_Class_Data %>% select("phylo") %>% unique()

Acclimation_Class_A_cor <- as.data.frame(A_cor)
Acclimation_Class_A_cor <- Acclimation_Class_A_cor[c(Acclimation_Class_Species$phylo), c(Acclimation_Class_Species$phylo)]
Acclimation_Class_A_cor <- as.matrix(Acclimation_Class_A_cor)

Acclimation_Class_VCV_InCVR <- make_VCV_matrix(Acclimation_Class_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Class_Model_CVR <- metafor::rma.mv(InCVR, V = Acclimation_Class_VCV_InCVR, test = "t", dfs = "contain",
                                                   mods = ~ Class - 1,
                                                   random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                 ~1|Shared_Animal_Number, ~1|Measurement), 
                                                   R = list(phylo=Acclimation_Class_A_cor), data = Acclimation_Class_Data, method = "REML", sparse = TRUE, 
                                                   control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Class_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Acclimation_Class_Model_CVR.rds")
  } else {
            Acclimation_Class_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Acclimation_Class_Model_CVR.rds")})

Acclimation_Class_Model_CVR_rob <- robust(Acclimation_Class_Model_CVR, cluster = Acclimation_Class_Data$Study_ID, adjust = TRUE)

Acclimation_Class_Model_CVR_Estimates <- data.frame(Class = substr(row.names(Acclimation_Class_Model_CVR$b), 6, 100),
                                                    estimate = Acclimation_Class_Model_CVR$b, 
                                                    ci.lb = Acclimation_Class_Model_CVR$ci.lb, 
                                                    ci.ub = Acclimation_Class_Model_CVR$ci.ub)
rownames(Acclimation_Class_Model_CVR_Estimates) <- Acclimation_Class_Model_CVR_Estimates$Class
Acclimation_Class_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Class_Model_CVR), 2))

# Preparing Graph - Combined

acclimation_class_rnames <- c("Actinopteri", "Bivalvia", "Gastropoda", "Holothuroidea", 
                              "Insecta", "Malacostraca")

acclimation_class_k <- data.frame("k" = c(Acclimation_Class_Exploration["Actinopteri", "Freq"],  
                                          Acclimation_Class_Exploration["Bivalvia", "Freq"], 
                                          Acclimation_Class_Exploration["Gastropoda", "Freq"], 
                                          Acclimation_Class_Exploration["Holothuroidea", "Freq"],
                                          Acclimation_Class_Exploration["Insecta", "Freq"],
                                          Acclimation_Class_Exploration["Malacostraca", "Freq"]), 
                                  row.names = acclimation_class_rnames)

acclimation_class_group_no <- data.frame("Spp No." = c(Acclimation_Class_Species_Count["Actinopteri", "Freq"],
                                                       Acclimation_Class_Species_Count["Bivalvia", "Freq"],
                                                       Acclimation_Class_Species_Count["Gastropoda", "Freq"], 
                                                       Acclimation_Class_Species_Count["Holothuroidea", "Freq"],
                                                       Acclimation_Class_Species_Count["Insecta", "Freq"],
                                                       Acclimation_Class_Species_Count["Malacostraca", "Freq"]), 
                                         row.names = acclimation_class_rnames)

acclimation_class_study <- data.frame("Study" = c(Acclimation_Class_Study_Count["Actinopteri", "Freq"],
                                                  Acclimation_Class_Study_Count["Bivalvia", "Freq"],
                                                  Acclimation_Class_Study_Count["Gastropoda", "Freq"], 
                                                  Acclimation_Class_Study_Count["Holothuroidea", "Freq"],
                                                  Acclimation_Class_Study_Count["Insecta", "Freq"],
                                                  Acclimation_Class_Study_Count["Malacostraca", "Freq"]), 
                                         row.names = acclimation_class_rnames)

acclimation_class_table <- data.frame(estimate = Acclimation_Class_Model_CVR_Estimates[,"estimate"], 
                                      lowerCL = Acclimation_Class_Model_CVR_Estimates[,"ci.lb"], 
                                      upperCL = Acclimation_Class_Model_CVR_Estimates[,"ci.ub"], 
                                      K = acclimation_class_k[,1], 
                                      group_no = acclimation_class_group_no[,1], 
                                      row.names = acclimation_class_rnames)
acclimation_class_table$name <- row.names(acclimation_class_table)

acclimation_class_raw_mean <- c(unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Actinopteri") %>% 
                                                select("InCVR"))),
                                unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Bivalvia") %>% 
                                                select("InCVR"))),
                                unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Gastropoda") %>% 
                                                select("InCVR"))), 
                                unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Holothuroidea") %>% 
                                                select("InCVR"))),
                                unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                                select("InCVR"))),
                                unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Malacostraca") %>% 
                                                select("InCVR"))))

acclimation_class_raw_name <- c(replicate(86, "Actinopteri"), 
                                replicate(12, "Bivalvia"),
                                replicate(21, "Gastropoda"), 
                                replicate(39, "Holothuroidea"),
                                replicate(103, "Insecta"),
                                replicate(48, "Malacostraca"))

acclimation_class_raw_df <- data.frame("Model" = acclimation_class_raw_name, 
                                       "Effect" = acclimation_class_raw_mean)

# Graph code - Combined

Acclimation_Class_Order <- c("Malacostraca", "Insecta", "Holothuroidea", "Gastropoda",  
                             "Bivalvia", "Actinopteri")

density_acclimation_class_CVR <- acclimation_class_table %>% mutate(name = fct_relevel(name, Acclimation_Class_Order)) %>%
                                 ggplot() +
                                 geom_density_ridges(data = acclimation_class_raw_df %>% mutate(Model = fct_relevel(Model, Acclimation_Class_Order)), 
                                                     aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                     scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                 geom_linerange(aes(y = rev(seq(1, dim(acclimation_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                size = 1) +
                                 geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                 size = 1, fatten = 2) +
                                 theme_bw() +
                                 guides(fill = "none", colour = "none") +
                                 labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                 theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                  vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7, -2.7))) +
                                 theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                 theme(axis.ticks = element_blank()) +
                                 theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                 theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                 scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                 scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                 scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                 coord_cartesian(xlim = c(-1, 1)) +
                                 annotate('text',  x = 1, y = (seq(1, dim(acclimation_class_table)[1], 1)+0.4),
                                 label= paste("italic(k)==", c(acclimation_class_table["Malacostraca", "K"], 
                                                               acclimation_class_table["Insecta", "K"], 
                                                               acclimation_class_table["Holothuroidea", "K"], 
                                                               acclimation_class_table["Gastropoda", "K"],
                                                               acclimation_class_table["Bivalvia", "K"],
                                                               acclimation_class_table["Actinopteri", "K"]), "~","(", 
                                                             c(acclimation_class_table["Malacostraca", "group_no"], 
                                                               acclimation_class_table["Insecta", "group_no"], 
                                                               acclimation_class_table["Holothuroidea", "group_no"], 
                                                               acclimation_class_table["Gastropoda", "group_no"],
                                                               acclimation_class_table["Bivalvia", "group_no"],
                                                               acclimation_class_table["Actinopteri", "group_no"]), 
                                              ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                 geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Malacostraca", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                        paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                        paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Holothuroidea", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                        paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Gastropoda", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                        paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Bivalvia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                        paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Actinopteri", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                            x = -0.75, y = (seq(1, dim(acclimation_class_table)[1], 1)+0.4)), size = 3.5)

density_acclimation_class_CVR #(400x560)

# Preparing Graph - Part 1

acclimation_class_rnames_1 <- c("Actinopteri", "Bivalvia", "Gastropoda")

acclimation_class_k_1 <- data.frame("k" = c(Acclimation_Class_Exploration["Actinopteri", "Freq"],  
                                            Acclimation_Class_Exploration["Bivalvia", "Freq"], 
                                            Acclimation_Class_Exploration["Gastropoda", "Freq"]), 
                                    row.names = acclimation_class_rnames_1)

acclimation_class_group_no_1 <- data.frame("Spp No." = c(Acclimation_Class_Species_Count["Actinopteri", "Freq"],
                                                         Acclimation_Class_Species_Count["Bivalvia", "Freq"],
                                                         Acclimation_Class_Species_Count["Gastropoda", "Freq"]), 
                                           row.names = acclimation_class_rnames_1)

acclimation_class_study_1 <- data.frame("Study" = c(Acclimation_Class_Study_Count["Actinopteri", "Freq"],
                                                    Acclimation_Class_Study_Count["Bivalvia", "Freq"],
                                                    Acclimation_Class_Study_Count["Gastropoda", "Freq"]), 
                                        row.names = acclimation_class_rnames_1)

Acclimation_Class_Model_CVR_Estimates_Reorder_1 <- Acclimation_Class_Model_CVR_Estimates[c("Actinopteri", "Bivalvia", "Gastropoda"), ]

acclimation_class_table_1 <- data.frame(estimate = Acclimation_Class_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                        lowerCL = Acclimation_Class_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                        upperCL = Acclimation_Class_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                        K = acclimation_class_k_1[,1], 
                                        group_no = acclimation_class_group_no_1[,1], 
                                        row.names = acclimation_class_rnames_1)
acclimation_class_table_1$name <- row.names(acclimation_class_table_1)

acclimation_class_raw_mean_1 <- c(unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Actinopteri") %>% 
                                                select("InCVR"))),
                                  unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Bivalvia") %>% 
                                                select("InCVR"))),
                                  unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Gastropoda") %>% 
                                                select("InCVR"))))

acclimation_class_raw_name_1 <- c(replicate(86, "Actinopteri"), 
                                  replicate(12, "Bivalvia"),
                                  replicate(21, "Gastropoda"))

acclimation_class_raw_df_1 <- data.frame("Model" = acclimation_class_raw_name_1, 
                                         "Effect" = acclimation_class_raw_mean_1)

# Graph code - Part 1

Acclimation_Class_Order_1 <- c("Gastropoda", "Bivalvia", "Actinopteri")

density_acclimation_class_CVR_1 <- acclimation_class_table_1 %>% mutate(name = fct_relevel(name, Acclimation_Class_Order_1)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = acclimation_class_raw_df_1 %>% mutate(Model = fct_relevel(Model, Acclimation_Class_Order_1)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                         vjust = c(-2.7, -2.7, -2.7))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                   scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(acclimation_class_table_1)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(acclimation_class_table_1["Gastropoda", "K"],
                                                                 acclimation_class_table_1["Bivalvia", "K"],
                                                                 acclimation_class_table_1["Actinopteri", "K"]), "~","(", 
                                                               c(acclimation_class_table_1["Gastropoda", "group_no"],
                                                                 acclimation_class_table_1["Bivalvia", "group_no"],
                                                                 acclimation_class_table_1["Actinopteri", "group_no"]), 
                                               ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Gastropoda", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                          paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Bivalvia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Actinopteri", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(acclimation_class_table_1)[1], 1)+0.4)), size = 3.5)

density_acclimation_class_CVR_1 #(400x320)

# Preparing Graph - Part 2

acclimation_class_rnames_2 <- c("Holothuroidea", "Insecta", "Malacostraca")

acclimation_class_k_2 <- data.frame("k" = c(Acclimation_Class_Exploration["Holothuroidea", "Freq"],
                                            Acclimation_Class_Exploration["Insecta", "Freq"],
                                            Acclimation_Class_Exploration["Malacostraca", "Freq"]), 
                                    row.names = acclimation_class_rnames_2)

acclimation_class_group_no_2 <- data.frame("Spp No." = c(Acclimation_Class_Species_Count["Holothuroidea", "Freq"],
                                                         Acclimation_Class_Species_Count["Insecta", "Freq"],
                                                         Acclimation_Class_Species_Count["Malacostraca", "Freq"]), 
                                           row.names = acclimation_class_rnames_2)

acclimation_class_study_2 <- data.frame("Study" = c(Acclimation_Class_Study_Count["Holothuroidea", "Freq"],
                                                    Acclimation_Class_Study_Count["Insecta", "Freq"],
                                                    Acclimation_Class_Study_Count["Malacostraca", "Freq"]), 
                                        row.names = acclimation_class_rnames_2)

Acclimation_Class_Model_CVR_Estimates_Reorder_2 <- Acclimation_Class_Model_CVR_Estimates[c("Holothuroidea", "Insecta", "Malacostraca"), ]

acclimation_class_table_2 <- data.frame(estimate = Acclimation_Class_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                        lowerCL = Acclimation_Class_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                        upperCL = Acclimation_Class_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                        K = acclimation_class_k_2[,1], 
                                        group_no = acclimation_class_group_no_2[,1], 
                                        row.names = acclimation_class_rnames_2)
acclimation_class_table_2$name <- row.names(acclimation_class_table_2)

acclimation_class_raw_mean_2 <- c(unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Holothuroidea") %>% 
                                                select("InCVR"))),
                                  unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                                select("InCVR"))),
                                  unlist(unname(Acclimation_Class_Data %>% filter(`Class` == "Malacostraca") %>% 
                                                select("InCVR"))))

acclimation_class_raw_name_2 <- c(replicate(39, "Holothuroidea"),
                                  replicate(103, "Insecta"),
                                  replicate(48, "Malacostraca"))

acclimation_class_raw_df_2 <- data.frame("Model" = acclimation_class_raw_name_2, 
                                         "Effect" = acclimation_class_raw_mean_2)

# Graph code - Part 2

Acclimation_Class_Order_2 <- c("Malacostraca", "Insecta", "Holothuroidea")

density_acclimation_class_CVR_2 <- acclimation_class_table_2 %>% mutate(name = fct_relevel(name, Acclimation_Class_Order_2)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = acclimation_class_raw_df_2 %>% mutate(Model = fct_relevel(Model, Acclimation_Class_Order_2)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                           scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_linerange(aes(y = rev(seq(1, dim(acclimation_class_table_2)[1], 1)), xmin = min(acclimation_class_raw_df_2$Effect)-0.05, xmax = -1.5, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                         vjust = c(-2.7, -2.7, -2.7))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                                   scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(acclimation_class_table_2)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(acclimation_class_table_2["Malacostraca", "K"], 
                                                                 acclimation_class_table_2["Insecta", "K"], 
                                                                 acclimation_class_table_2["Holothuroidea", "K"]), "~","(", 
                                                               c(acclimation_class_table_2["Malacostraca", "group_no"], 
                                                                 acclimation_class_table_2["Insecta", "group_no"], 
                                                                 acclimation_class_table_2["Holothuroidea", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Malacostraca", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                          paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Acclimation_Class_Model_CVR_Estimates["Holothuroidea", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(acclimation_class_table_2)[1], 1)+0.4)), size = 3.5)

density_acclimation_class_CVR_2 #(400x320)

##### Acclimation Subset Model - Specific Trait Meta-Regression - CVR #####
Acclimation_Specific_Trait_Exploration <- Acclimation_Subset_Data %>% select("Measurement") %>% table() %>% data.frame()
Acclimation_Specific_Trait_Exploration <- Acclimation_Specific_Trait_Exploration %>% filter(Freq > 10)
rownames(Acclimation_Specific_Trait_Exploration) <- Acclimation_Specific_Trait_Exploration$Measurement

Acclimation_Specific_Trait_Data <- Acclimation_Subset_Data %>% filter(Measurement == "Apparent Digestability Coefficient"| 
                                                                      Measurement == "Catalase Activity"|
                                                                      Measurement == "Cortisol"|
                                                                      Measurement == "Food Consumption"|
                                                                      Measurement == "Immune Defense"|
                                                                      Measurement == "Metabolic Rate"|
                                                                      Measurement == "SOD Activity")

Acclimation_Specific_Trait_Species_Count <- Acclimation_Specific_Trait_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>% 
                                            filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Acclimation_Specific_Trait_Species_Count) <- Acclimation_Specific_Trait_Species_Count$Measurement

Acclimation_Specific_Trait_Study_Count <- Acclimation_Specific_Trait_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>% 
                                          filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Acclimation_Specific_Trait_Study_Count) <- Acclimation_Specific_Trait_Study_Count$Measurement

Acclimation_Specific_Trait_Species <- Acclimation_Specific_Trait_Data %>% select("phylo") %>% unique()

Acclimation_Specific_Trait_A_cor <- as.data.frame(A_cor)
Acclimation_Specific_Trait_A_cor <- Acclimation_Specific_Trait_A_cor[c(Acclimation_Specific_Trait_Species$phylo), c(Acclimation_Specific_Trait_Species$phylo)]
Acclimation_Specific_Trait_A_cor <- as.matrix(Acclimation_Specific_Trait_A_cor)

Acclimation_Specific_Trait_VCV_InCVR <- make_VCV_matrix(Acclimation_Specific_Trait_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  1ish minutes
  if(run){
    Acclimation_Specific_Trait_Model_CVR <- metafor::rma.mv(InCVR, V = Acclimation_Specific_Trait_VCV_InCVR, test = "t", dfs = "contain",
                                                            mods = ~ Measurement - 1,
                                                            random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                          ~1|Shared_Animal_Number), 
                                                            R = list(phylo=Acclimation_Specific_Trait_A_cor), data = Acclimation_Specific_Trait_Data, method = "REML", sparse = TRUE, 
                                                            control=list(rel.tol=1e-9))
    saveRDS(Acclimation_Specific_Trait_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Acclimation_Specific_Trait_Model_CVR.rds")
  } else {
            Acclimation_Specific_Trait_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Acclimation_Specific_Trait_Model_CVR.rds")})

Acclimation_Specific_Trait_Model_CVR_rob <- robust(Acclimation_Specific_Trait_Model_CVR, cluster = Acclimation_Specific_Trait_Data$Study_ID, adjust = TRUE)

Acclimation_Specific_Trait_Model_CVR_Estimates <- data.frame(Trait = substr(row.names(Acclimation_Specific_Trait_Model_CVR$b), 12, 100),
                                                             estimate = Acclimation_Specific_Trait_Model_CVR$b, 
                                                             ci.lb = Acclimation_Specific_Trait_Model_CVR$ci.lb, 
                                                             ci.ub = Acclimation_Specific_Trait_Model_CVR$ci.ub)
rownames(Acclimation_Specific_Trait_Model_CVR_Estimates) <- Acclimation_Specific_Trait_Model_CVR_Estimates$Trait
Acclimation_Specific_Trait_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Acclimation_Specific_Trait_Model_CVR), 2))

# Preparing Graph - Combined

acclimation_specific_trait_rnames <- c("Apparent Digestibility Coefficient", "Catalase Activity", "Cortisol Levels", 
                                       "Food Consumption", "Immune Defense", "Metabolic Rate", "SOD Activity")

acclimation_specific_trait_k <- data.frame("k" = c(Acclimation_Specific_Trait_Exploration["Apparent Digestability Coefficient", "Freq"], 
                                                   Acclimation_Specific_Trait_Exploration["Catalase Activity", "Freq"], 
                                                   Acclimation_Specific_Trait_Exploration["Cortisol", "Freq"], 
                                                   Acclimation_Specific_Trait_Exploration["Food Consumption", "Freq"],
                                                   Acclimation_Specific_Trait_Exploration["Immune Defense", "Freq"],
                                                   Acclimation_Specific_Trait_Exploration["Metabolic Rate", "Freq"], 
                                                   Acclimation_Specific_Trait_Exploration["SOD Activity", "Freq"]), 
                                           row.names = acclimation_specific_trait_rnames)

acclimation_specific_trait_group_no <- data.frame("Spp No." = c(Acclimation_Specific_Trait_Species_Count["Apparent Digestability Coefficient", "Freq"], 
                                                                Acclimation_Specific_Trait_Species_Count["Catalase Activity", "Freq"], 
                                                                Acclimation_Specific_Trait_Species_Count["Cortisol", "Freq"], 
                                                                Acclimation_Specific_Trait_Species_Count["Food Consumption", "Freq"],
                                                                Acclimation_Specific_Trait_Species_Count["Immune Defense", "Freq"],
                                                                Acclimation_Specific_Trait_Species_Count["Metabolic Rate", "Freq"], 
                                                                Acclimation_Specific_Trait_Species_Count["SOD Activity", "Freq"]), 
                                                  row.names = acclimation_specific_trait_rnames)

acclimation_specific_trait_study <- data.frame("Study" = c(Acclimation_Specific_Trait_Study_Count["Apparent Digestability Coefficient", "Freq"], 
                                                           Acclimation_Specific_Trait_Study_Count["Catalase Activity", "Freq"], 
                                                           Acclimation_Specific_Trait_Study_Count["Cortisol", "Freq"], 
                                                           Acclimation_Specific_Trait_Study_Count["Food Consumption", "Freq"],
                                                           Acclimation_Specific_Trait_Study_Count["Immune Defense", "Freq"],
                                                           Acclimation_Specific_Trait_Study_Count["Metabolic Rate", "Freq"], 
                                                           Acclimation_Specific_Trait_Study_Count["SOD Activity", "Freq"]), 
                                                  row.names = acclimation_specific_trait_rnames)

acclimation_specific_trait_table <- data.frame(estimate = Acclimation_Specific_Trait_Model_CVR_Estimates[,"estimate"], 
                                               lowerCL = Acclimation_Specific_Trait_Model_CVR_Estimates[,"ci.lb"], 
                                               upperCL = Acclimation_Specific_Trait_Model_CVR_Estimates[,"ci.ub"], 
                                               K = acclimation_specific_trait_k[,1], 
                                               group_no = acclimation_specific_trait_group_no[,1], 
                                               row.names = acclimation_specific_trait_rnames)
acclimation_specific_trait_table$name <- row.names(acclimation_specific_trait_table)

acclimation_specific_trait_raw_mean <- c(unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Apparent Digestability Coefficient") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Catalase Activity") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Cortisol") %>% 
                                                         select("InCVR"))), 
                                         unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Food Consumption") %>% 
                                                         select("InCVR"))),
                                         unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Immune Defense") %>% 
                                                         select("InCVR"))),
                                         unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Metabolic Rate") %>% 
                                                         select("InCVR"))),
                                         unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "SOD Activity") %>% 
                                                         select("InCVR"))))

acclimation_specific_trait_raw_name <- c(replicate(16, "Apparent Digestibility Coefficient"), 
                                         replicate(11, "Catalase Activity"), 
                                         replicate(13, "Cortisol Levels"),
                                         replicate(12, "Food Consumption"),
                                         replicate(14, "Immune Defense"),
                                         replicate(41, "Metabolic Rate"),
                                         replicate(11, "SOD Activity"))

acclimation_specific_trait_raw_df <- data.frame("Model" = acclimation_specific_trait_raw_name, 
                                                "Effect" = acclimation_specific_trait_raw_mean)

# Graph code - Combined

Acclimation_Specific_Trait_Order <- c("SOD Activity", "Metabolic Rate", "Immune Defense", 
                                      "Food Consumption", "Cortisol Levels", "Catalase Activity", 
                                      "Apparent Digestibility Coefficient")

density_acclimation_specific_trait_CVR <- acclimation_specific_trait_table %>% mutate(name = fct_relevel(name, Acclimation_Specific_Trait_Order)) %>%
                                          ggplot() +
                                          geom_density_ridges(data = acclimation_specific_trait_raw_df %>% mutate(Model = fct_relevel(Model, Acclimation_Specific_Trait_Order)), 
                                                              aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                              scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                          geom_linerange(aes(y = rev(seq(1, dim(acclimation_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                         size = 1) +
                                          geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                          size = 1, fatten = 2) +
                                          theme_bw() +
                                          guides(fill = "none", colour = "none") +
                                          labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                          theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                           vjust = c(-2.7, -0.8, -0.8, -0.8, -0.8, -0.8, -0.4))) +
                                          theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                          theme(axis.ticks = element_blank()) +
                                          theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                          theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                          scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                          scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                                          scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                                          coord_cartesian(xlim = c(-1, 1)) +
                                          annotate('text',  x = 1, y = (seq(1, dim(acclimation_specific_trait_table)[1], 1)+0.4),
                                          label= paste("italic(k)==", c(acclimation_specific_trait_table["SOD Activity", "K"],
                                                                        acclimation_specific_trait_table["Metabolic Rate", "K"],
                                                                        acclimation_specific_trait_table["Immune Defense", "K"], 
                                                                        acclimation_specific_trait_table["Food Consumption", "K"], 
                                                                        acclimation_specific_trait_table["Cortisol Levels", "K"], 
                                                                        acclimation_specific_trait_table["Catalase Activity", "K"],
                                                                        acclimation_specific_trait_table["Apparent Digestibility Coefficient", "K"]), "~","(", 
                                                                      c(acclimation_specific_trait_table["SOD Activity", "group_no"],
                                                                        acclimation_specific_trait_table["Metabolic Rate", "group_no"],
                                                                        acclimation_specific_trait_table["Immune Defense", "group_no"],
                                                                        acclimation_specific_trait_table["Food Consumption", "group_no"], 
                                                                        acclimation_specific_trait_table["Cortisol Levels", "group_no"], 
                                                                        acclimation_specific_trait_table["Catalase Activity", "group_no"],
                                                                        acclimation_specific_trait_table["Apparent Digestibility Coefficient", "group_no"]), 
                                                       ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                          geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["SOD Activity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Immune Defense", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                 paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Food Consumption", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Cortisol", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                 paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Catalase Activity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                 paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Apparent Digestability Coefficient", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                     x = -0.75, y = (seq(1, dim(acclimation_specific_trait_table)[1], 1)+0.4)), size = 3.5)

density_acclimation_specific_trait_CVR #(400x640)

# Preparing Graph - Part 1

acclimation_specific_trait_rnames_1 <- c("Apparent Digestibility Coefficient", "Catalase Activity", "Cortisol Levels", "Food Consumption")

acclimation_specific_trait_k_1 <- data.frame("k" = c(Acclimation_Specific_Trait_Exploration["Apparent Digestability Coefficient", "Freq"], 
                                                     Acclimation_Specific_Trait_Exploration["Catalase Activity", "Freq"], 
                                                     Acclimation_Specific_Trait_Exploration["Cortisol", "Freq"], 
                                                     Acclimation_Specific_Trait_Exploration["Food Consumption", "Freq"]), 
                                             row.names = acclimation_specific_trait_rnames_1)

acclimation_specific_trait_group_no_1 <- data.frame("Spp No." = c(Acclimation_Specific_Trait_Species_Count["Apparent Digestability Coefficient", "Freq"], 
                                                                  Acclimation_Specific_Trait_Species_Count["Catalase Activity", "Freq"], 
                                                                  Acclimation_Specific_Trait_Species_Count["Cortisol", "Freq"], 
                                                                  Acclimation_Specific_Trait_Species_Count["Food Consumption", "Freq"]), 
                                                    row.names = acclimation_specific_trait_rnames_1)

acclimation_specific_trait_study_1 <- data.frame("Study" = c(Acclimation_Specific_Trait_Study_Count["Apparent Digestability Coefficient", "Freq"], 
                                                             Acclimation_Specific_Trait_Study_Count["Catalase Activity", "Freq"], 
                                                             Acclimation_Specific_Trait_Study_Count["Cortisol", "Freq"], 
                                                             Acclimation_Specific_Trait_Study_Count["Food Consumption", "Freq"]), 
                                                 row.names = acclimation_specific_trait_rnames_1)

Acclimation_Specific_Trait_Model_CVR_Estimates_Reorder_1 <- Acclimation_Specific_Trait_Model_CVR_Estimates[c("Apparent Digestability Coefficient", "Catalase Activity", "Cortisol", "Food Consumption"), ]

acclimation_specific_trait_table_1 <- data.frame(estimate = Acclimation_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                                 lowerCL = Acclimation_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                                 upperCL = Acclimation_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                                 K = acclimation_specific_trait_k_1[,1], 
                                                 group_no = acclimation_specific_trait_group_no_1[,1], 
                                                 row.names = acclimation_specific_trait_rnames_1)
acclimation_specific_trait_table_1$name <- row.names(acclimation_specific_trait_table_1)

acclimation_specific_trait_raw_mean_1 <- c(unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Apparent Digestability Coefficient") %>% 
                                                         select("InCVR"))), 
                                           unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Catalase Activity") %>% 
                                                         select("InCVR"))), 
                                           unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Cortisol") %>% 
                                                         select("InCVR"))), 
                                           unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Food Consumption") %>% 
                                                         select("InCVR"))))

acclimation_specific_trait_raw_name_1 <- c(replicate(16, "Apparent Digestibility Coefficient"), 
                                           replicate(11, "Catalase Activity"), 
                                           replicate(13, "Cortisol Levels"),
                                           replicate(12, "Food Consumption"))

acclimation_specific_trait_raw_df_1 <- data.frame("Model" = acclimation_specific_trait_raw_name_1, 
                                                  "Effect" = acclimation_specific_trait_raw_mean_1)

# Graph code - Part 1

Acclimation_Specific_Trait_Order_1 <- c("Food Consumption", "Cortisol Levels", "Catalase Activity", 
                                        "Apparent Digestibility Coefficient")

density_acclimation_specific_trait_CVR_1 <- acclimation_specific_trait_table_1 %>% mutate(name = fct_relevel(name, Acclimation_Specific_Trait_Order_1)) %>%
                                            ggplot() +
                                            geom_density_ridges(data = acclimation_specific_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Acclimation_Specific_Trait_Order_1)), 
                                                                aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                            geom_linerange(aes(y = rev(seq(1, dim(acclimation_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                           size = 1) +
                                            geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                            size = 1, fatten = 2) +
                                            theme_bw() +
                                            guides(fill = "none", colour = "none") +
                                            labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                            theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                  vjust = c(-0.8, -0.8, -0.8, -0.4))) +
                                            theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                            theme(axis.ticks = element_blank()) +
                                            theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                            theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                            scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                            scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                                            scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                                            coord_cartesian(xlim = c(-1, 1)) +
                                            annotate('text',  x = 1, y = (seq(1, dim(acclimation_specific_trait_table_1)[1], 1)+0.4),
                                            label= paste("italic(k)==", c(acclimation_specific_trait_table_1["Food Consumption", "K"], 
                                                                          acclimation_specific_trait_table_1["Cortisol Levels", "K"], 
                                                                          acclimation_specific_trait_table_1["Catalase Activity", "K"],
                                                                          acclimation_specific_trait_table_1["Apparent Digestibility Coefficient", "K"]), "~","(", 
                                                                        c(acclimation_specific_trait_table_1["Food Consumption", "group_no"], 
                                                                          acclimation_specific_trait_table_1["Cortisol Levels", "group_no"], 
                                                                          acclimation_specific_trait_table_1["Catalase Activity", "group_no"],
                                                                          acclimation_specific_trait_table_1["Apparent Digestibility Coefficient", "group_no"]), 
                                                         ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                            geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Food Consumption", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Cortisol", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                   paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Catalase Activity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Apparent Digestability Coefficient", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                       x = -0.75, y = (seq(1, dim(acclimation_specific_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_acclimation_specific_trait_CVR_1 #(400x400)

# Preparing Graph - Part 2

acclimation_specific_trait_rnames_2 <- c("Immune Defense", "Metabolic Rate", "SOD Activity")

acclimation_specific_trait_k_2 <- data.frame("k" = c(Acclimation_Specific_Trait_Exploration["Immune Defense", "Freq"],
                                                     Acclimation_Specific_Trait_Exploration["Metabolic Rate", "Freq"], 
                                                     Acclimation_Specific_Trait_Exploration["SOD Activity", "Freq"]), 
                                             row.names = acclimation_specific_trait_rnames_2)

acclimation_specific_trait_group_no_2 <- data.frame("Spp No." = c(Acclimation_Specific_Trait_Species_Count["Immune Defense", "Freq"],
                                                                  Acclimation_Specific_Trait_Species_Count["Metabolic Rate", "Freq"], 
                                                                  Acclimation_Specific_Trait_Species_Count["SOD Activity", "Freq"]), 
                                                    row.names = acclimation_specific_trait_rnames_2)

acclimation_specific_trait_study_2 <- data.frame("Study" = c(Acclimation_Specific_Trait_Study_Count["Immune Defense", "Freq"],
                                                             Acclimation_Specific_Trait_Study_Count["Metabolic Rate", "Freq"], 
                                                             Acclimation_Specific_Trait_Study_Count["SOD Activity", "Freq"]), 
                                                 row.names = acclimation_specific_trait_rnames_2)

Acclimation_Specific_Trait_Model_CVR_Estimates_Reorder_2 <- Acclimation_Specific_Trait_Model_CVR_Estimates[c("Immune Defense", "Metabolic Rate", "SOD Activity"), ]

acclimation_specific_trait_table_2 <- data.frame(estimate = Acclimation_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                                 lowerCL = Acclimation_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                                 upperCL = Acclimation_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                                 K = acclimation_specific_trait_k_2[,1], 
                                                 group_no = acclimation_specific_trait_group_no_2[,1], 
                                                 row.names = acclimation_specific_trait_rnames_2)
acclimation_specific_trait_table_2$name <- row.names(acclimation_specific_trait_table_2)

acclimation_specific_trait_raw_mean_2 <- c(unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Immune Defense") %>% 
                                                         select("InCVR"))),
                                           unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "Metabolic Rate") %>% 
                                                         select("InCVR"))),
                                           unlist(unname(Acclimation_Specific_Trait_Data %>% filter(`Measurement` == "SOD Activity") %>% 
                                                         select("InCVR"))))

acclimation_specific_trait_raw_name_2 <- c(replicate(14, "Immune Defense"),
                                           replicate(41, "Metabolic Rate"),
                                           replicate(11, "SOD Activity"))

acclimation_specific_trait_raw_df_2 <- data.frame("Model" = acclimation_specific_trait_raw_name_2, 
                                                  "Effect" = acclimation_specific_trait_raw_mean_2)

# Graph code - Part 2

Acclimation_Specific_Trait_Order_2 <- c("SOD Activity", "Metabolic Rate", "Immune Defense")

density_acclimation_specific_trait_CVR_2 <- acclimation_specific_trait_table_2 %>% mutate(name = fct_relevel(name, Acclimation_Specific_Trait_Order_2)) %>%
                                            ggplot() +
                                            geom_density_ridges(data = acclimation_specific_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Acclimation_Specific_Trait_Order_2)), 
                                                                aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                            geom_linerange(aes(y = rev(seq(1, dim(acclimation_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                           size = 1) +
                                            geom_linerange(aes(y = rev(seq(1, dim(acclimation_specific_trait_table_2)[1], 1)), xmin = max(acclimation_specific_trait_raw_df_2$Effect)+0.25, xmax = 1.5, colour = name),
                                                           size = 1) +
                                            geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(acclimation_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                            size = 1, fatten = 2) +
                                            theme_bw() +
                                            guides(fill = "none", colour = "none") +
                                            labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                            theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                  vjust = c(-2.7, -0.8, -0.8))) +
                                            theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                            theme(axis.ticks = element_blank()) +
                                            theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                            theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                            scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                            scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                                            scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                                            coord_cartesian(xlim = c(-1, 1)) +
                                            annotate('text',  x = 1, y = (seq(1, dim(acclimation_specific_trait_table_2)[1], 1)+0.4),
                                            label= paste("italic(k)==", c(acclimation_specific_trait_table_2["SOD Activity", "K"],
                                                                          acclimation_specific_trait_table_2["Metabolic Rate", "K"],
                                                                          acclimation_specific_trait_table_2["Immune Defense", "K"]), "~","(", 
                                                                        c(acclimation_specific_trait_table_2["SOD Activity", "group_no"],
                                                                          acclimation_specific_trait_table_2["Metabolic Rate", "group_no"],
                                                                          acclimation_specific_trait_table_2["Immune Defense", "group_no"]), 
                                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                            geom_label(aes(label=c(paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["SOD Activity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Metabolic Rate", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Acclimation_Specific_Trait_Model_CVR_Estimates["Immune Defense", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                      x = -0.75, y = (seq(1, dim(acclimation_specific_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_acclimation_specific_trait_CVR_2 #(400x320)

##### Summary of Acclimation Plots #####

Acclimation_Layout <- rbind(c(1, 3), 
                            c(1, 3),
                            c(1, 3),
                            c(1, 3),
                            c(1, 3),
                            c(1, 3),
                            c(1, 4),
                            c(1, 4),
                            c(2, 4),
                            c(2, 5),
                            c(2, 5),
                            c(2, 5),
                            c(2, 5))

Acclimation_Combined_CVR <- grid.arrange(density_acclimation_specific_trait_CVR, density_acclimation_trait_CVR, density_acclimation_class_CVR, 
                                         density_acclimation_fluctuation_CVR, density_acclimation_stage_CVR, 
                                         layout_matrix = Acclimation_Layout)

Acclimation_Combined_CVR #(850 x 900 - does not include amplitude plot)

##### Developmental Subset Model - CVR #####
Developmental_Subset_Data <- Individual_Subset_Data %>% filter(Plasticity_Mechanism == "Developmental Plasticity")
Developmental_Species <- Developmental_Subset_Data %>% select("phylo") %>% unique()

Developmental_A_cor <- as.data.frame(A_cor)
Developmental_A_cor <- Developmental_A_cor[c(Developmental_Species$phylo), c(Developmental_Species$phylo)]
Developmental_A_cor <- as.matrix(Developmental_A_cor)

Developmental_VCV_InCVR <- make_VCV_matrix(Developmental_Subset_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  12ish minutes
  if(run){
    Developmental_Model_CVR <- metafor::rma.mv(InCVR ~ 1, V = Developmental_VCV_InCVR, test = "t", dfs = "contain",
                                               random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                             ~1|Shared_Animal_Number, ~1|Measurement), 
                                               R = list(phylo=Developmental_A_cor), data = Developmental_Subset_Data, method = "REML", sparse = TRUE, 
                                               control=list(rel.tol=1e-9))
    saveRDS(Developmental_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Developmental_Model_CVR.rds")
  } else {
            Developmental_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Developmental_Model_CVR.rds")})

Developmental_Model_CVR_rob <- robust(Developmental_Model_CVR, cluster = Developmental_Subset_Data$Study_ID, adjust = TRUE)

Developmental_Model_CVR_Estimates <- data.frame(estimate = Developmental_Model_CVR$b, 
                                                ci.lb = Developmental_Model_CVR$ci.lb, 
                                                ci.ub = Developmental_Model_CVR$ci.ub)
Developmental_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Model_CVR), 2))

#### Developmental Subset Model - Fluctuation Amplitude Meta-Regression - CVR ####
run <- FALSE
system.time( #  4ish minutes
  if(run){
    Developmental_Amplitude_Model_CVR <- metafor::rma.mv(InCVR, V = Developmental_VCV_InCVR, test = "t", dfs = "contain",
                                                         mods = ~ T2_Magnitude - 1,
                                                         random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                       ~1|Shared_Animal_Number, ~1|Measurement), 
                                                         R = list(phylo=Developmental_A_cor), data = Developmental_Subset_Data, method = "REML", sparse = TRUE, 
                                                         control=list(rel.tol=1e-9))
    saveRDS(Developmental_Amplitude_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Developmental_Amplitude_Model_CVR.rds")
  } else {
            Developmental_Amplitude_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Developmental_Amplitude_Model_CVR.rds")})

Developmental_Amplitude_Model_CVR_rob <- robust(Developmental_Amplitude_Model_CVR, cluster = Developmental_Subset_Data$Study_ID, adjust = TRUE)

Developmental_Amplitude_Model_CVR_Estimates <- data.frame(estimate = Developmental_Amplitude_Model_CVR$b, 
                                                          ci.lb = Developmental_Amplitude_Model_CVR$ci.lb, 
                                                          ci.ub = Developmental_Amplitude_Model_CVR$ci.ub)
Developmental_Amplitude_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Amplitude_Model_CVR), 2))

# Graph Preparing

Developmental_Plot_Data <- Developmental_Subset_Data
Developmental_Plot_Data <- Developmental_Plot_Data %>% mutate(n_category = ifelse(n_R1.1 <= 25, "25", 
                                                                           ifelse(n_R1.1 > 25 & n_R1.1 <= 50, "50", 
                                                                           ifelse(n_R1.1 > 50 & n_R1.1 <= 75, "75", "> 75"))))

# Graph Code

Developmental_Amplitude_Plot_CVR <- ggplot(Developmental_Plot_Data, aes(x = T2_Magnitude, y = InCVR)) + 
                                    geom_point(aes(x = T2_Magnitude, y = InCVR, 
                                               size = fct_relevel(n_category, c("25", "50", "75", "> 75"))), 
                                               shape = 21, fill = "#4292c6", alpha = 0.5) + 
                                    labs(x = "Fluctuation Amplitude (\u00B0C)", y = "Effect Size (lnCVR)", 
                                         size = "Sample Size", title = "Developmental Treatments") +
                                    theme_bw() +
                                    theme(plot.title = element_text(size = 12, colour ="black", face = "bold", hjust = 0.5, margin = margin(b = 10))) +
                                    theme(axis.text.y = element_text(size = 10, colour ="black", margin = margin(l = 5))) +
                                    theme(axis.text.x = element_text(size = 10, colour ="black", margin = margin(b = 10))) +
                                    theme(legend.position = "bottom", legend.direction = "horizontal") + 
                                    geom_hline(yintercept = Developmental_Model_CVR_Estimates$estimate, lty = 2) + 
                                    geom_smooth(method = "lm", linewidth = 1, se = F, colour = "#084594") +
                                    stat_poly_eq(formula = y ~ x, 
                                                 aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                                                 parse = TRUE) +
                                    coord_cartesian(xlim = c(0, 30), 
                                                    ylim = c(-2.5, 2.5))

Developmental_Amplitude_Plot_CVR #(400x400)

#### Developmental Subset Model - Type of Fluctuation Meta-Regression - CVR ####
Developmental_Fluctuation_Data <- Developmental_Subset_Data %>% filter(!is.na(Fluctuation_Category))
Developmental_Fluctuation_Exploration <- Developmental_Subset_Data %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Developmental_Fluctuation_Exploration) <- Developmental_Fluctuation_Exploration$Fluctuation_Category

Developmental_Fluctuation_Species_Count <- Developmental_Subset_Data %>% select("Scientific_Name", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                           filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Developmental_Fluctuation_Species_Count) <- Developmental_Fluctuation_Species_Count$Fluctuation_Category

Developmental_Fluctuation_Study_Count <- Developmental_Subset_Data %>% select("Study_ID", "Fluctuation_Category") %>% table() %>% data.frame() %>%
                                         filter(`Freq` != 0) %>% select("Fluctuation_Category") %>% table() %>% data.frame()
rownames(Developmental_Fluctuation_Study_Count) <- Developmental_Fluctuation_Study_Count$Fluctuation_Category

Developmental_Fluctuation_Species <- Developmental_Fluctuation_Data %>% select("phylo") %>% unique()

Developmental_Fluctuation_A_cor <- as.data.frame(A_cor)
Developmental_Fluctuation_A_cor <- Developmental_Fluctuation_A_cor[c(Developmental_Fluctuation_Species$phylo), c(Developmental_Fluctuation_Species$phylo)]
Developmental_Fluctuation_A_cor <- as.matrix(Developmental_Fluctuation_A_cor)

Developmental_Fluctuation_VCV_InCVR <- make_VCV_matrix(Developmental_Fluctuation_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  11ish minutes
  if(run){
    Developmental_Fluctuation_Model_CVR <- metafor::rma.mv(InCVR, V = Developmental_Fluctuation_VCV_InCVR, test = "t", dfs = "contain",
                                                           mods = ~ Fluctuation_Category - 1,
                                                           random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                         ~1|Shared_Animal_Number, ~1|Measurement), 
                                                           R = list(phylo=Developmental_Fluctuation_A_cor), data = Developmental_Fluctuation_Data, method = "REML", sparse = TRUE, 
                                                           control=list(rel.tol=1e-9))
    saveRDS(Developmental_Fluctuation_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Developmental_Fluctuation_Model_CVR.rds")
  } else {
            Developmental_Fluctuation_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Developmental_Fluctuation_Model_CVR.rds")})

Developmental_Fluctuation_Model_CVR_rob <- robust(Developmental_Fluctuation_Model_CVR, cluster = Developmental_Fluctuation_Data$Study_ID, adjust = TRUE)

Developmental_Fluctuation_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Developmental_Fluctuation_Model_CVR$b), 21, 100),
                                                            estimate = Developmental_Fluctuation_Model_CVR$b, 
                                                            ci.lb = Developmental_Fluctuation_Model_CVR$ci.lb, 
                                                            ci.ub = Developmental_Fluctuation_Model_CVR$ci.ub)
rownames(Developmental_Fluctuation_Model_CVR_Estimates) <- Developmental_Fluctuation_Model_CVR_Estimates$Category
Developmental_Fluctuation_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Fluctuation_Model_CVR), 2))

# Preparing Graph - Combined

developmental_fluctuation_rnames <- c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise", "Stochastic")

developmental_fluctuation_k <- data.frame("k" = c(Developmental_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                  Developmental_Fluctuation_Exploration["Alternating", "Freq"], 
                                                  Developmental_Fluctuation_Exploration["Stepwise", "Freq"], 
                                                  Developmental_Fluctuation_Exploration["Stochastic", "Freq"]), 
                                          row.names = developmental_fluctuation_rnames)

developmental_fluctuation_group_no <- data.frame("Spp No." = c(Developmental_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                               Developmental_Fluctuation_Species_Count["Alternating", "Freq"], 
                                                               Developmental_Fluctuation_Species_Count["Stepwise", "Freq"], 
                                                               Developmental_Fluctuation_Species_Count["Stochastic", "Freq"]), 
                                                 row.names = developmental_fluctuation_rnames)

developmental_fluctuation_study <- data.frame("Study" = c(Developmental_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                          Developmental_Fluctuation_Study_Count["Alternating", "Freq"], 
                                                          Developmental_Fluctuation_Study_Count["Stepwise", "Freq"], 
                                                          Developmental_Fluctuation_Study_Count["Stochastic", "Freq"]), 
                                                 row.names = developmental_fluctuation_rnames)

Developmental_Fluctuation_Model_CVR_Estimates_Reorder <- Developmental_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating", "Stepwise", "Stochastic"), ]

developmental_fluctuation_table <- data.frame(estimate = Developmental_Fluctuation_Model_CVR_Estimates_Reorder[,"estimate"], 
                                              lowerCL = Developmental_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.lb"], 
                                              upperCL = Developmental_Fluctuation_Model_CVR_Estimates_Reorder[,"ci.ub"], 
                                              K = developmental_fluctuation_k[,1], 
                                              group_no = developmental_fluctuation_group_no[,1], 
                                              row.names = developmental_fluctuation_rnames)
developmental_fluctuation_table$name <- row.names(developmental_fluctuation_table)

developmental_fluctuation_raw_mean <- c(unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                        select("InCVR"))), 
                                        unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                        select("InCVR"))), 
                                        unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                        select("InCVR"))), 
                                        unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stochastic") %>% 
                                                        select("InCVR"))))

developmental_fluctuation_raw_name <- c(replicate(330, "Sinusoidal (Sine Curve)"), 
                                        replicate(458, "Alternating"), 
                                        replicate(120, "Stepwise"), 
                                        replicate(17, "Stochastic"))

developmental_fluctuation_raw_df <- data.frame("Model" = developmental_fluctuation_raw_name, 
                                               "Effect" = developmental_fluctuation_raw_mean)

# Graph code - Combined

Developmental_Fluctuation_Order <- c("Stochastic", "Stepwise", 
                                     "Alternating", "Sinusoidal (Sine Curve)")

density_developmental_fluctuation_CVR <- developmental_fluctuation_table %>% mutate(name = fct_relevel(name, Developmental_Fluctuation_Order)) %>%
                                         ggplot() +
                                         geom_density_ridges(data = developmental_fluctuation_raw_df %>% mutate(Model = fct_relevel(Model, Developmental_Fluctuation_Order)), 
                                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                         geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                        size = 1) +
                                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_fluctuation_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                         size = 1, fatten = 2) +
                                         theme_bw() +
                                         guides(fill = "none", colour = "none") +
                                         labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                         theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                          vjust = c(-2.7, -2.7, -2.7, -0.8))) +
                                         theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                         theme(axis.ticks = element_blank()) +
                                         theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                         scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                         scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                         scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                         coord_cartesian(xlim = c(-1, 1)) +
                                         annotate('text',  x = 1, y = (seq(1, dim(developmental_fluctuation_table)[1], 1)+0.4),
                                         label= paste("italic(k)==", c(developmental_fluctuation_table["Stochastic", "K"], 
                                                                       developmental_fluctuation_table["Stepwise", "K"], 
                                                                       developmental_fluctuation_table["Alternating", "K"], 
                                                                       developmental_fluctuation_table["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                     c(developmental_fluctuation_table["Stochastic", "group_no"], 
                                                                       developmental_fluctuation_table["Stepwise", "group_no"], 
                                                                       developmental_fluctuation_table["Alternating", "group_no"], 
                                                                       developmental_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"]), 
                                                      ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                         geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Fluctuation_Model_CVR_Estimates["Stochastic", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                paste(format(round(mean(exp(Developmental_Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                paste(format(round(mean(exp(Developmental_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                paste(format(round(mean(exp(Developmental_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                    x = -0.75, y = (seq(1, dim(developmental_fluctuation_table)[1], 1)+0.4)), size = 3.5)

density_developmental_fluctuation_CVR #(400x400)

# Preparing Graph - Part 1

developmental_fluctuation_rnames_1 <- c("Sinusoidal (Sine Curve)", "Alternating")

developmental_fluctuation_k_1 <- data.frame("k" = c(Developmental_Fluctuation_Exploration["Sinusoidal (Sine Curve)", "Freq"], 
                                                    Developmental_Fluctuation_Exploration["Alternating", "Freq"]), 
                                            row.names = developmental_fluctuation_rnames_1)

developmental_fluctuation_group_no_1 <- data.frame("Spp No." = c(Developmental_Fluctuation_Species_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                                 Developmental_Fluctuation_Species_Count["Alternating", "Freq"]), 
                                                   row.names = developmental_fluctuation_rnames_1)

developmental_fluctuation_study_1 <- data.frame("Study" = c(Developmental_Fluctuation_Study_Count["Sinusoidal (Sine Curve)", "Freq"], 
                                                            Developmental_Fluctuation_Study_Count["Alternating", "Freq"]), 
                                                row.names = developmental_fluctuation_rnames_1)

Developmental_Fluctuation_Model_CVR_Estimates_Reorder_1 <- Developmental_Fluctuation_Model_CVR_Estimates[c("Sinusoidal (Sine Curve)", "Alternating"), ]

developmental_fluctuation_table_1 <- data.frame(estimate = Developmental_Fluctuation_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                                lowerCL = Developmental_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                                upperCL = Developmental_Fluctuation_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                                K = developmental_fluctuation_k_1[,1], 
                                                group_no = developmental_fluctuation_group_no_1[,1], 
                                                row.names = developmental_fluctuation_rnames_1)
developmental_fluctuation_table_1$name <- row.names(developmental_fluctuation_table_1)

developmental_fluctuation_raw_mean_1 <- c(unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Sinusoidal (Sine Curve)") %>% 
                                                        select("InCVR"))), 
                                          unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Alternating") %>% 
                                                        select("InCVR"))))

developmental_fluctuation_raw_name_1 <- c(replicate(330, "Sinusoidal (Sine Curve)"), 
                                          replicate(458, "Alternating"))

developmental_fluctuation_raw_df_1 <- data.frame("Model" = developmental_fluctuation_raw_name_1, 
                                                 "Effect" = developmental_fluctuation_raw_mean_1)

# Graph code - Part 1

Developmental_Fluctuation_Order_1 <- c("Alternating", "Sinusoidal (Sine Curve)")

density_developmental_fluctuation_CVR_1 <- developmental_fluctuation_table_1 %>% mutate(name = fct_relevel(name, Developmental_Fluctuation_Order_1)) %>%
                                           ggplot() +
                                           geom_density_ridges(data = developmental_fluctuation_raw_df_1 %>% mutate(Model = fct_relevel(Model, Developmental_Fluctuation_Order_1)), 
                                                               aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                   scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                           geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                          size = 1) +
                                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_fluctuation_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                           size = 1, fatten = 2) +
                                           theme_bw() +
                                           guides(fill = "none", colour = "none") +
                                           labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                           theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                 vjust = c(-2.7, -0.8))) +
                                           theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                           theme(axis.ticks = element_blank()) +
                                           theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                           theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                           scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                           scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                           scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                           coord_cartesian(xlim = c(-1, 1)) +
                                           annotate('text',  x = 1, y = (seq(1, dim(developmental_fluctuation_table_1)[1], 1)+0.4),
                                           label= paste("italic(k)==", c(developmental_fluctuation_table_1["Alternating", "K"], 
                                                                         developmental_fluctuation_table_1["Sinusoidal (Sine Curve)", "K"]), "~","(", 
                                                                       c(developmental_fluctuation_table_1["Alternating", "group_no"], 
                                                                         developmental_fluctuation_table_1["Sinusoidal (Sine Curve)", "group_no"]), 
                                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                           geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Fluctuation_Model_CVR_Estimates["Alternating", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                  paste(format(round(mean(exp(Developmental_Fluctuation_Model_CVR_Estimates["Sinusoidal (Sine Curve)", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                      x = -0.75, y = (seq(1, dim(developmental_fluctuation_table_1)[1], 1)+0.4)), size = 3.5)

density_developmental_fluctuation_CVR_1 #(400x240)

# Preparing Graph - Part 2

developmental_fluctuation_rnames_2 <- c("Stepwise", "Stochastic")

developmental_fluctuation_k_2 <- data.frame("k" = c(Developmental_Fluctuation_Exploration["Stepwise", "Freq"], 
                                                    Developmental_Fluctuation_Exploration["Stochastic", "Freq"]), 
                                            row.names = developmental_fluctuation_rnames_2)

developmental_fluctuation_group_no_2 <- data.frame("Spp No." = c(Developmental_Fluctuation_Species_Count["Stepwise", "Freq"], 
                                                                 Developmental_Fluctuation_Species_Count["Stochastic", "Freq"]), 
                                                   row.names = developmental_fluctuation_rnames_2)

developmental_fluctuation_study_2 <- data.frame("Study" = c(Developmental_Fluctuation_Study_Count["Stepwise", "Freq"], 
                                                            Developmental_Fluctuation_Study_Count["Stochastic", "Freq"]), 
                                                row.names = developmental_fluctuation_rnames_2)

Developmental_Fluctuation_Model_CVR_Estimates_Reorder_2 <- Developmental_Fluctuation_Model_CVR_Estimates[c("Stepwise", "Stochastic"), ]

developmental_fluctuation_table_2 <- data.frame(estimate = Developmental_Fluctuation_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                                lowerCL = Developmental_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                                upperCL = Developmental_Fluctuation_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                                K = developmental_fluctuation_k_2[,1], 
                                                group_no = developmental_fluctuation_group_no_2[,1], 
                                                row.names = developmental_fluctuation_rnames_2)
developmental_fluctuation_table_2$name <- row.names(developmental_fluctuation_table_2)

developmental_fluctuation_raw_mean_2 <- c(unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stepwise") %>% 
                                                        select("InCVR"))), 
                                          unlist(unname(Developmental_Fluctuation_Data %>% filter(`Fluctuation_Category` == "Stochastic") %>% 
                                                        select("InCVR"))))

developmental_fluctuation_raw_name_2 <- c(replicate(120, "Stepwise"), 
                                          replicate(17, "Stochastic"))

developmental_fluctuation_raw_df_2 <- data.frame("Model" = developmental_fluctuation_raw_name_2, 
                                                 "Effect" = developmental_fluctuation_raw_mean_2)

# Graph code - Part 2

Developmental_Fluctuation_Order_2 <- c("Stochastic", "Stepwise")

density_developmental_fluctuation_CVR_2 <- developmental_fluctuation_table_2 %>% mutate(name = fct_relevel(name, Developmental_Fluctuation_Order_2)) %>%
                                           ggplot() +
                                           geom_density_ridges(data = developmental_fluctuation_raw_df_2 %>% mutate(Model = fct_relevel(Model, Developmental_Fluctuation_Order_2)), 
                                                               aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                   scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                           geom_linerange(aes(y = rev(seq(1, dim(developmental_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                          size = 1) +
                                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_fluctuation_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                           size = 1, fatten = 2) +
                                           theme_bw() +
                                           guides(fill = "none", colour = "none") +
                                           labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                           theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                 vjust = c(-2.7, -2.7))) +
                                           theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                           theme(axis.ticks = element_blank()) +
                                           theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                           theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                           scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                           scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                           scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                           coord_cartesian(xlim = c(-1, 1)) +
                                           annotate('text',  x = 1, y = (seq(1, dim(developmental_fluctuation_table_2)[1], 1)+0.4),
                                           label= paste("italic(k)==", c(developmental_fluctuation_table_2["Stochastic", "K"], 
                                                                         developmental_fluctuation_table_2["Stepwise", "K"]), "~","(", 
                                                                       c(developmental_fluctuation_table_2["Stochastic", "group_no"], 
                                                                         developmental_fluctuation_table_2["Stepwise", "group_no"]), 
                                                        ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                           geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Fluctuation_Model_CVR_Estimates["Stochastic", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                  paste(format(round(mean(exp(Developmental_Fluctuation_Model_CVR_Estimates["Stepwise", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                      x = -0.75, y = (seq(1, dim(developmental_fluctuation_table_2)[1], 1)+0.4)), size = 3.5)

density_developmental_fluctuation_CVR_2 #(400x240)

##### Developmental Subset Model - Trait Meta-Regression - CVR #####
Developmental_Trait_Exploration <- Developmental_Subset_Data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Developmental_Trait_Exploration) <- Developmental_Trait_Exploration$Trait_Category

Developmental_Trait_Species_Count <- Developmental_Subset_Data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>% 
                                     filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Developmental_Trait_Species_Count) <- Developmental_Trait_Species_Count$Trait_Category

Developmental_Trait_Study_Count <- Developmental_Subset_Data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>% 
                                   filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Developmental_Trait_Study_Count) <- Developmental_Trait_Study_Count$Trait_Category

run <- FALSE
system.time( #  12ish minutes
  if(run){
    Developmental_Trait_Model_CVR <- metafor::rma.mv(InCVR, V = Developmental_VCV_InCVR, test = "t", dfs = "contain",
                                                     mods = ~ Trait_Category - 1,
                                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                   ~1|Shared_Animal_Number, ~1|Measurement), 
                                                     R = list(phylo=Developmental_A_cor), data = Developmental_Subset_Data, method = "REML", sparse = TRUE, 
                                                     control=list(rel.tol=1e-9))
    saveRDS(Developmental_Trait_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Developmental_Trait_Model_CVR.rds")
  } else {
            Developmental_Trait_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Developmental_Trait_Model_CVR.rds")})

Developmental_Trait_Model_CVR_rob <- robust(Developmental_Trait_Model_CVR, cluster = Developmental_Subset_Data$Study_ID, adjust = TRUE)

Developmental_Trait_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Developmental_Trait_Model_CVR$b), 15, 100),
                                                      estimate = Developmental_Trait_Model_CVR$b, 
                                                      ci.lb = Developmental_Trait_Model_CVR$ci.lb, 
                                                      ci.ub = Developmental_Trait_Model_CVR$ci.ub)
rownames(Developmental_Trait_Model_CVR_Estimates) <- Developmental_Trait_Model_CVR_Estimates$Category
Developmental_Trait_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Trait_Model_CVR), 2))

# Preparing Graph - Combined

developmental_trait_rnames <- c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits", 
                                "Morphology", "Physiological")

developmental_trait_k <- data.frame("k" = c(Developmental_Trait_Exploration["Behavioural", "Freq"], 
                                            Developmental_Trait_Exploration["Biochemical Assay", "Freq"], 
                                            Developmental_Trait_Exploration["Gene Expression", "Freq"], 
                                            Developmental_Trait_Exploration["Life-History Traits", "Freq"], 
                                            Developmental_Trait_Exploration["Morphology", "Freq"], 
                                            Developmental_Trait_Exploration["Physiological", "Freq"]), 
                                    row.names = developmental_trait_rnames)

developmental_trait_group_no <- data.frame("Spp No." = c(Developmental_Trait_Species_Count["Behavioural", "Freq"], 
                                                         Developmental_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                         Developmental_Trait_Species_Count["Gene Expression", "Freq"], 
                                                         Developmental_Trait_Species_Count["Life-History Traits", "Freq"],
                                                         Developmental_Trait_Species_Count["Morphology", "Freq"],
                                                         Developmental_Trait_Species_Count["Physiological", "Freq"]), 
                                           row.names = developmental_trait_rnames)

developmental_trait_study <- data.frame("Study" = c(Developmental_Trait_Study_Count["Behavioural", "Freq"], 
                                                    Developmental_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                    Developmental_Trait_Study_Count["Gene Expression", "Freq"], 
                                                    Developmental_Trait_Study_Count["Life-History Traits", "Freq"],
                                                    Developmental_Trait_Study_Count["Morphology", "Freq"],
                                                    Developmental_Trait_Study_Count["Physiological", "Freq"]), 
                                           row.names = developmental_trait_rnames)

developmental_trait_table <- data.frame(estimate = Developmental_Trait_Model_CVR_Estimates[,"estimate"], 
                                        lowerCL = Developmental_Trait_Model_CVR_Estimates[,"ci.lb"], 
                                        upperCL = Developmental_Trait_Model_CVR_Estimates[,"ci.ub"], 
                                        K = developmental_trait_k[,1], 
                                        group_no = developmental_trait_group_no[,1], 
                                        row.names = developmental_trait_rnames)
developmental_trait_table$name <- row.names(developmental_trait_table)

developmental_trait_raw_mean <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Behavioural") %>% 
                                                  select("InCVR"))), 
                                  unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                                  select("InCVR"))), 
                                  unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Gene Expression") %>% 
                                                  select("InCVR"))), 
                                  unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                                  select("InCVR"))), 
                                  unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                                  select("InCVR"))),
                                  unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                                  select("InCVR"))))

developmental_trait_raw_name <- c(replicate(14, "Behavioural"), 
                                  replicate(14, "Biochemical Assay"), 
                                  replicate(40, "Gene Expression"), 
                                  replicate(480, "Life-history Traits"), 
                                  replicate(376, "Morphology"),
                                  replicate(64, "Physiological"))

developmental_trait_raw_df <- data.frame("Model" = developmental_trait_raw_name, 
                                         "Effect" = developmental_trait_raw_mean)

# Graph code - Combined

Developmental_Trait_Order <- c("Physiological", "Morphology", "Life-history Traits",  
                               "Gene Expression", "Biochemical Assay", "Behavioural")

density_developmental_trait_CVR <- developmental_trait_table %>% mutate(name = fct_relevel(name, Developmental_Trait_Order)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = developmental_trait_raw_df %>% mutate(Model = fct_relevel(Model, Developmental_Trait_Order)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                    vjust = c(-2.7, -2.7, -0.8, -0.8, -0.8, -2.7))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                   scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(developmental_trait_table)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(developmental_trait_table["Physiological", "K"], 
                                                                 developmental_trait_table["Morphology", "K"], 
                                                                 developmental_trait_table["Life-history Traits", "K"],
                                                                 developmental_trait_table["Gene Expression", "K"],
                                                                 developmental_trait_table["Biochemical Assay", "K"],
                                                                 developmental_trait_table["Behavioural", "K"]), "~","(", 
                                                               c(developmental_trait_table["Physiological", "group_no"], 
                                                                 developmental_trait_table["Morphology", "group_no"], 
                                                                 developmental_trait_table["Life-history Traits", "group_no"],
                                                                 developmental_trait_table["Gene Expression", "group_no"],
                                                                 developmental_trait_table["Biochemical Assay", "group_no"],
                                                                 developmental_trait_table["Behavioural", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                          paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                          paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Gene Expression", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Behavioural", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(developmental_trait_table)[1], 1)+0.4)), size = 3.5)

density_developmental_trait_CVR #(400x560)

# Preparing Graph - Part 1

developmental_trait_rnames_1 <- c("Behavioural", "Biochemical Assay", "Gene Expression")

developmental_trait_k_1 <- data.frame("k" = c(Developmental_Trait_Exploration["Behavioural", "Freq"], 
                                              Developmental_Trait_Exploration["Biochemical Assay", "Freq"], 
                                              Developmental_Trait_Exploration["Gene Expression", "Freq"]), 
                                      row.names = developmental_trait_rnames_1)

developmental_trait_group_no_1 <- data.frame("Spp No." = c(Developmental_Trait_Species_Count["Behavioural", "Freq"], 
                                                           Developmental_Trait_Species_Count["Biochemical Assay", "Freq"], 
                                                           Developmental_Trait_Species_Count["Gene Expression", "Freq"]), 
                                             row.names = developmental_trait_rnames_1)

developmental_trait_study_1 <- data.frame("Study" = c(Developmental_Trait_Study_Count["Behavioural", "Freq"], 
                                                      Developmental_Trait_Study_Count["Biochemical Assay", "Freq"], 
                                                      Developmental_Trait_Study_Count["Gene Expression", "Freq"]), 
                                          row.names = developmental_trait_rnames_1)

Developmental_Trait_Model_CVR_Estimates_Reorder_1 <- Developmental_Trait_Model_CVR_Estimates[c("Behavioural", "Biochemical Assay", "Gene Expression"), ]

developmental_trait_table_1 <- data.frame(estimate = Developmental_Trait_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                          lowerCL = Developmental_Trait_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                          upperCL = Developmental_Trait_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                          K = developmental_trait_k_1[,1], 
                                          group_no = developmental_trait_group_no_1[,1], 
                                          row.names = developmental_trait_rnames_1)
developmental_trait_table_1$name <- row.names(developmental_trait_table_1)

developmental_trait_raw_mean_1 <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Behavioural") %>% 
                                                  select("InCVR"))), 
                                    unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Biochemical Assay") %>% 
                                                  select("InCVR"))), 
                                    unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Gene Expression") %>% 
                                                  select("InCVR"))))

developmental_trait_raw_name_1 <- c(replicate(14, "Behavioural"), 
                                    replicate(14, "Biochemical Assay"), 
                                    replicate(40, "Gene Expression"))

developmental_trait_raw_df_1 <- data.frame("Model" = developmental_trait_raw_name_1, 
                                           "Effect" = developmental_trait_raw_mean_1)

# Graph code - Part 1

Developmental_Trait_Order_1 <- c("Gene Expression", "Biochemical Assay", "Behavioural")

density_developmental_trait_CVR_1 <- developmental_trait_table_1 %>% mutate(name = fct_relevel(name, Developmental_Trait_Order_1)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = developmental_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Developmental_Trait_Order_1)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                     size = 1, fatten = 2) +
                                     theme_bw() +
                                     guides(fill = "none", colour = "none") +
                                     labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                           vjust = c(-0.8, -0.8, -2.7))) +
                                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                     theme(axis.ticks = element_blank()) +
                                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                     scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                     scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                     coord_cartesian(xlim = c(-1, 1)) +
                                     annotate('text',  x = 1, y = (seq(1, dim(developmental_trait_table_1)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(developmental_trait_table["Gene Expression", "K"],
                                                                   developmental_trait_table["Biochemical Assay", "K"],
                                                                   developmental_trait_table["Behavioural", "K"]), "~","(", 
                                                                 c(developmental_trait_table["Gene Expression", "group_no"],
                                                                   developmental_trait_table["Biochemical Assay", "group_no"],
                                                                   developmental_trait_table["Behavioural", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Gene Expression", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                            paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Biochemical Assay", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                            paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Behavioural", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.75, y = (seq(1, dim(developmental_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_developmental_trait_CVR_1 #(400x320)

# Preparing Graph - Part 2

developmental_trait_rnames_2 <- c("Life-history Traits", "Morphology", "Physiological")

developmental_trait_k_2 <- data.frame("k" = c(Developmental_Trait_Exploration["Life-History Traits", "Freq"], 
                                              Developmental_Trait_Exploration["Morphology", "Freq"], 
                                              Developmental_Trait_Exploration["Physiological", "Freq"]), 
                                      row.names = developmental_trait_rnames_2)

developmental_trait_group_no_2 <- data.frame("Spp No." = c(Developmental_Trait_Species_Count["Life-History Traits", "Freq"],
                                                           Developmental_Trait_Species_Count["Morphology", "Freq"],
                                                           Developmental_Trait_Species_Count["Physiological", "Freq"]), 
                                             row.names = developmental_trait_rnames_2)

developmental_trait_study_2 <- data.frame("Study" = c(Developmental_Trait_Study_Count["Life-History Traits", "Freq"],
                                                      Developmental_Trait_Study_Count["Morphology", "Freq"],
                                                      Developmental_Trait_Study_Count["Physiological", "Freq"]), 
                                          row.names = developmental_trait_rnames_2)

Developmental_Trait_Model_CVR_Estimates_Reorder_2 <- Developmental_Trait_Model_CVR_Estimates[c("Life-History Traits", "Morphology", "Physiological"), ]

developmental_trait_table_2 <- data.frame(estimate = Developmental_Trait_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                          lowerCL = Developmental_Trait_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                          upperCL = Developmental_Trait_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                          K = developmental_trait_k_2[,1], 
                                          group_no = developmental_trait_group_no_2[,1], 
                                          row.names = developmental_trait_rnames_2)
developmental_trait_table_2$name <- row.names(developmental_trait_table_2)

developmental_trait_raw_mean_2 <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Life-History Traits") %>% 
                                                  select("InCVR"))), 
                                    unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Morphology") %>% 
                                                  select("InCVR"))),
                                    unlist(unname(Developmental_Subset_Data %>% filter(`Trait_Category` == "Physiological") %>% 
                                                  select("InCVR"))))

developmental_trait_raw_name_2 <- c(replicate(480, "Life-history Traits"), 
                                    replicate(376, "Morphology"),
                                    replicate(64, "Physiological"))

developmental_trait_raw_df_2 <- data.frame("Model" = developmental_trait_raw_name_2, 
                                           "Effect" = developmental_trait_raw_mean_2)

# Graph code - Part 2

Developmental_Trait_Order_2 <- c("Physiological", "Morphology", "Life-history Traits")

density_developmental_trait_CVR_2 <- developmental_trait_table_2 %>% mutate(name = fct_relevel(name, Developmental_Trait_Order_2)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = developmental_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Developmental_Trait_Order_2)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(developmental_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                     size = 1, fatten = 2) +
                                     theme_bw() +
                                     guides(fill = "none", colour = "none") +
                                     labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                           vjust = c(-2.7, -2.7, -0.8))) +
                                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                     theme(axis.ticks = element_blank()) +
                                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                     scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                                     scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                                     coord_cartesian(xlim = c(-1, 1)) +
                                     annotate('text',  x = 1, y = (seq(1, dim(developmental_trait_table_2)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(developmental_trait_table_2["Physiological", "K"], 
                                                                   developmental_trait_table_2["Morphology", "K"], 
                                                                   developmental_trait_table_2["Life-history Traits", "K"]), "~","(", 
                                                                 c(developmental_trait_table_2["Physiological", "group_no"], 
                                                                   developmental_trait_table_2["Morphology", "group_no"], 
                                                                   developmental_trait_table_2["Life-history Traits", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Physiological", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                            paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Morphology", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                            paste(format(round(mean(exp(Developmental_Trait_Model_CVR_Estimates["Life-History Traits", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.75, y = (seq(1, dim(developmental_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_developmental_trait_CVR_2 #(400x320)

##### Developmental Subset Model - Exposure Time Meta-Regression - CVR #####
Developmental_Exposure_Exploration <- Developmental_Subset_Data %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Exploration) <- Developmental_Exposure_Exploration$Developmental_Exposure_Time_Category

Developmental_Exposure_Species_Count <- Developmental_Subset_Data %>% select("Scientific_Name", "Developmental_Exposure_Time_Category") %>% table() %>% data.frame() %>% 
                                        filter(`Freq` != 0) %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Species_Count) <- Developmental_Exposure_Species_Count$Developmental_Exposure_Time_Category

Developmental_Exposure_Study_Count <- Developmental_Subset_Data %>% select("Study_ID", "Developmental_Exposure_Time_Category") %>% table() %>% data.frame() %>% 
                                      filter(`Freq` != 0) %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Study_Count) <- Developmental_Exposure_Study_Count$Developmental_Exposure_Time_Category

run <- FALSE
system.time( #  13ish minutes
  if(run){
    Developmental_Exposure_Model_CVR <- metafor::rma.mv(InCVR, V = Developmental_VCV_InCVR, test = "t", dfs = "contain",
                                                        mods = ~ Developmental_Exposure_Time_Category - 1,
                                                        random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                      ~1|Shared_Animal_Number, ~1|Measurement), 
                                                        R = list(phylo=Developmental_A_cor), data = Developmental_Subset_Data, method = "REML", sparse = TRUE, 
                                                        control=list(rel.tol=1e-9))
    saveRDS(Developmental_Exposure_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Developmental_Exposure_Model_CVR.rds")
  } else {
            Developmental_Exposure_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Developmental_Exposure_Model_CVR.rds")})

Developmental_Exposure_Model_CVR_rob <- robust(Developmental_Exposure_Model_CVR, cluster = Developmental_Subset_Data$Study_ID, adjust = TRUE)

Developmental_Exposure_Model_CVR_Estimates <- data.frame(Category = substr(row.names(Developmental_Exposure_Model_CVR$b), 37, 100),
                                                         estimate = Developmental_Exposure_Model_CVR$b, 
                                                         ci.lb = Developmental_Exposure_Model_CVR$ci.lb, 
                                                         ci.ub = Developmental_Exposure_Model_CVR$ci.ub)
rownames(Developmental_Exposure_Model_CVR_Estimates) <- Developmental_Exposure_Model_CVR_Estimates$Category
Developmental_Exposure_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Exposure_Model_CVR), 2))

# Preparing Graph - Combined

developmental_exposure_rnames <- c("Embryo", "Juvenile", "Larva", "Pupa")

developmental_exposure_k <- data.frame("k" = c(Developmental_Exposure_Exploration["Embryo", "Freq"], 
                                               Developmental_Exposure_Exploration["Juvenile", "Freq"], 
                                               Developmental_Exposure_Exploration["Larvae", "Freq"], 
                                               Developmental_Exposure_Exploration["Pupae", "Freq"]), 
                                       row.names = developmental_exposure_rnames)

developmental_exposure_group_no <- data.frame("Spp No." = c(Developmental_Exposure_Species_Count["Embryo", "Freq"], 
                                                            Developmental_Exposure_Species_Count["Juvenile", "Freq"], 
                                                            Developmental_Exposure_Species_Count["Larvae", "Freq"], 
                                                            Developmental_Exposure_Species_Count["Pupae", "Freq"]), 
                                              row.names = developmental_exposure_rnames)

developmental_exposure_study <- data.frame("Study" = c(Developmental_Exposure_Study_Count["Embryo", "Freq"], 
                                                       Developmental_Exposure_Study_Count["Juvenile", "Freq"], 
                                                       Developmental_Exposure_Study_Count["Larvae", "Freq"], 
                                                       Developmental_Exposure_Study_Count["Pupae", "Freq"]), 
                                              row.names = developmental_exposure_rnames)

developmental_exposure_table <- data.frame(estimate = Developmental_Exposure_Model_CVR_Estimates[,"estimate"], 
                                           lowerCL = Developmental_Exposure_Model_CVR_Estimates[,"ci.lb"], 
                                           upperCL = Developmental_Exposure_Model_CVR_Estimates[,"ci.ub"], 
                                           K = developmental_exposure_k[,1], 
                                           group_no = developmental_exposure_group_no[,1], 
                                           row.names = developmental_exposure_rnames)
developmental_exposure_table$name <- row.names(developmental_exposure_table)

developmental_exposure_raw_mean <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Embryo") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Juvenile") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Larvae") %>% 
                                                     select("InCVR"))), 
                                     unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Pupae") %>% 
                                                     select("InCVR"))))

developmental_exposure_raw_name <- c(replicate(598, "Embryo"), 
                                     replicate(50, "Juvenile"), 
                                     replicate(311, "Larva"), 
                                     replicate(29, "Pupa"))

developmental_exposure_raw_df <- data.frame("Model" = developmental_exposure_raw_name, 
                                            "Effect" = developmental_exposure_raw_mean)

# Graph code - Combined

Developmental_Exposure_Order <- c("Pupa", "Larva", 
                                  "Juvenile", "Embryo")

density_developmental_exposure_CVR <- developmental_exposure_table %>% mutate(name = fct_relevel(name, Developmental_Exposure_Order)) %>%
                                      ggplot() +
                                      geom_density_ridges(data = developmental_exposure_raw_df %>% mutate(Model = fct_relevel(Model, Developmental_Exposure_Order)), 
                                                          aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                          scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                      geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                     size = 1) +
                                      geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_exposure_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                      size = 1, fatten = 2) +
                                      theme_bw() +
                                      guides(fill = "none", colour = "none") +
                                      labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                      theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                       vjust = c(-2.7, -2.7, -2.7, -2.7))) +
                                      theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                      theme(axis.ticks = element_blank()) +
                                      theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                      scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                      scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                      scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A")) +
                                      coord_cartesian(xlim = c(-1, 1)) +
                                      annotate('text',  x = 1, y = (seq(1, dim(developmental_exposure_table)[1], 1)+0.4),
                                      label= paste("italic(k)==", c(developmental_exposure_table["Pupa", "K"], 
                                                                    developmental_exposure_table["Larva", "K"], 
                                                                    developmental_exposure_table["Juvenile", "K"], 
                                                                    developmental_exposure_table["Embryo", "K"]), "~","(", 
                                                                  c(developmental_exposure_table["Pupa", "group_no"], 
                                                                    developmental_exposure_table["Larva", "group_no"], 
                                                                    developmental_exposure_table["Juvenile", "group_no"], 
                                                                    developmental_exposure_table["Embryo", "group_no"]), 
                                                   ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                      geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Exposure_Model_CVR_Estimates["Pupae", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                             paste(format(round(mean(exp(Developmental_Exposure_Model_CVR_Estimates["Larvae", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                             paste(format(round(mean(exp(Developmental_Exposure_Model_CVR_Estimates["Juvenile", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                             paste(format(round(mean(exp(Developmental_Exposure_Model_CVR_Estimates["Embryo", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                 x = -0.75, y = (seq(1, dim(developmental_exposure_table)[1], 1)+0.4)), size = 3.5)

density_developmental_exposure_CVR #(400x400)

# Preparing Graph - Part 1

developmental_exposure_rnames_1 <- c("Embryo", "Juvenile")

developmental_exposure_k_1 <- data.frame("k" = c(Developmental_Exposure_Exploration["Embryo", "Freq"], 
                                                 Developmental_Exposure_Exploration["Juvenile", "Freq"]), 
                                         row.names = developmental_exposure_rnames_1)

developmental_exposure_group_no_1 <- data.frame("Spp No." = c(Developmental_Exposure_Species_Count["Embryo", "Freq"], 
                                                              Developmental_Exposure_Species_Count["Juvenile", "Freq"]), 
                                                row.names = developmental_exposure_rnames_1)

developmental_exposure_study_1 <- data.frame("Study" = c(Developmental_Exposure_Study_Count["Embryo", "Freq"], 
                                                         Developmental_Exposure_Study_Count["Juvenile", "Freq"]), 
                                             row.names = developmental_exposure_rnames_1)

Developmental_Exposure_Model_CVR_Estimates_Reorder_1 <- Developmental_Exposure_Model_CVR_Estimates[c("Embryo", "Juvenile"), ]

developmental_exposure_table_1 <- data.frame(estimate = Developmental_Exposure_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                             lowerCL = Developmental_Exposure_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                             upperCL = Developmental_Exposure_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                             K = developmental_exposure_k_1[,1], 
                                             group_no = developmental_exposure_group_no_1[,1], 
                                             row.names = developmental_exposure_rnames_1)
developmental_exposure_table_1$name <- row.names(developmental_exposure_table_1)

developmental_exposure_raw_mean_1 <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Embryo") %>% 
                                                     select("InCVR"))), 
                                       unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Juvenile") %>% 
                                                     select("InCVR"))))

developmental_exposure_raw_name_1 <- c(replicate(598, "Embryo"), 
                                       replicate(50, "Juvenile"))

developmental_exposure_raw_df_1 <- data.frame("Model" = developmental_exposure_raw_name_1, 
                                              "Effect" = developmental_exposure_raw_mean_1)

# Graph code - Part 1

Developmental_Exposure_Order_1 <- c("Juvenile", "Embryo")

density_developmental_exposure_CVR_1 <- developmental_exposure_table_1 %>% mutate(name = fct_relevel(name, Developmental_Exposure_Order_1)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = developmental_exposure_raw_df_1 %>% mutate(Model = fct_relevel(Model, Developmental_Exposure_Order_1)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_exposure_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                        size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                              vjust = c(-2.7, -2.7))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                        scale_fill_manual(values = c("#3C5F8D", "#2B4E7A")) +
                                        coord_cartesian(xlim = c(-1, 1)) +
                                        annotate('text',  x = 1, y = (seq(1, dim(developmental_exposure_table_1)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(developmental_exposure_table_1["Juvenile", "K"], 
                                                                      developmental_exposure_table_1["Embryo", "K"]), "~","(", 
                                                                    c(developmental_exposure_table_1["Juvenile", "group_no"], 
                                                                      developmental_exposure_table_1["Embryo", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Exposure_Model_CVR_Estimates["Juvenile", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                               paste(format(round(mean(exp(Developmental_Exposure_Model_CVR_Estimates["Embryo", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.75, y = (seq(1, dim(developmental_exposure_table_1)[1], 1)+0.4)), size = 3.5)

density_developmental_exposure_CVR_1 #(400x240)

# Preparing Graph - Part 2

developmental_exposure_rnames_2 <- c("Larva", "Pupa")

developmental_exposure_k_2 <- data.frame("k" = c(Developmental_Exposure_Exploration["Larvae", "Freq"], 
                                                 Developmental_Exposure_Exploration["Pupae", "Freq"]), 
                                         row.names = developmental_exposure_rnames_2)

developmental_exposure_group_no_2 <- data.frame("Spp No." = c(Developmental_Exposure_Species_Count["Larvae", "Freq"], 
                                                              Developmental_Exposure_Species_Count["Pupae", "Freq"]), 
                                                row.names = developmental_exposure_rnames_2)

developmental_exposure_study_2 <- data.frame("Study" = c(Developmental_Exposure_Study_Count["Larvae", "Freq"], 
                                                         Developmental_Exposure_Study_Count["Pupae", "Freq"]), 
                                             row.names = developmental_exposure_rnames_2)

Developmental_Exposure_Model_CVR_Estimates_Reorder_2 <- Developmental_Exposure_Model_CVR_Estimates[c("Larva", "Pupa"), ]

developmental_exposure_table_2 <- data.frame(estimate = Developmental_Exposure_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                             lowerCL = Developmental_Exposure_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                             upperCL = Developmental_Exposure_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                             K = developmental_exposure_k_2[,1], 
                                             group_no = developmental_exposure_group_no_2[,1], 
                                             row.names = developmental_exposure_rnames_2)
developmental_exposure_table_2$name <- row.names(developmental_exposure_table_2)

developmental_exposure_raw_mean_2 <- c(unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Larvae") %>% 
                                                     select("InCVR"))), 
                                       unlist(unname(Developmental_Subset_Data %>% filter(`Developmental_Exposure_Time_Category` == "Pupae") %>% 
                                                     select("InCVR"))))

developmental_exposure_raw_name_2 <- c(replicate(311, "Larva"), 
                                       replicate(29, "Pupa"))

developmental_exposure_raw_df_2 <- data.frame("Model" = developmental_exposure_raw_name_2, 
                                              "Effect" = developmental_exposure_raw_mean_2)

# Graph code - Part 2

Developmental_Exposure_Order_2 <- c("Pupa", "Larva")

density_developmental_exposure_CVR_2 <- developmental_exposure_table_2 %>% mutate(name = fct_relevel(name, Developmental_Exposure_Order_2)) %>%
                                        ggplot() +
                                        geom_density_ridges(data = developmental_exposure_raw_df_2 %>% mutate(Model = fct_relevel(Model, Developmental_Exposure_Order_2)), 
                                                            aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                            scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                        geom_linerange(aes(y = rev(seq(1, dim(developmental_exposure_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                       size = 1) +
                                        geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_exposure_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                        size = 1, fatten = 2) +
                                        theme_bw() +
                                        guides(fill = "none", colour = "none") +
                                        labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                        theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                              vjust = c(-2.7, -2.7))) +
                                        theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                        theme(axis.ticks = element_blank()) +
                                        theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                        scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                        scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                        scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                        coord_cartesian(xlim = c(-1, 1)) +
                                        annotate('text',  x = 1, y = (seq(1, dim(developmental_exposure_table_2)[1], 1)+0.4),
                                        label= paste("italic(k)==", c(developmental_exposure_table_2["Pupa", "K"], 
                                                                      developmental_exposure_table_2["Larva", "K"]), "~","(", 
                                                                    c(developmental_exposure_table_2["Pupa", "group_no"], 
                                                                      developmental_exposure_table_2["Larva", "group_no"]), 
                                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                        geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Exposure_Model_CVR_Estimates["Pupae", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                               paste(format(round(mean(exp(Developmental_Exposure_Model_CVR_Estimates["Larvae", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                   x = -0.75, y = (seq(1, dim(developmental_exposure_table_2)[1], 1)+0.4)), size = 3.5)

density_developmental_exposure_CVR_2 #(400x240)

##### Developmental Subset Model - Class Meta-Regression - CVR #####
Developmental_Class_Exploration <- Developmental_Subset_Data %>% select("Class") %>% table() %>% data.frame()
rownames(Developmental_Class_Exploration) <- Developmental_Class_Exploration$Class

Developmental_Class_Data <- Developmental_Subset_Data %>% filter(Class != "Anthozoa" &
                                                                 Class != "Malacostraca")

Developmental_Class_Species_Count <- Developmental_Class_Data %>% select("Scientific_Name", "Class") %>% table() %>% data.frame() %>%
                                     filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Developmental_Class_Species_Count) <- Developmental_Class_Species_Count$Class

Developmental_Class_Study_Count <- Developmental_Class_Data %>% select("Study_ID", "Class") %>% table() %>% data.frame() %>%
                                   filter(`Freq` != 0) %>% select("Class") %>% table() %>% data.frame()
rownames(Developmental_Class_Study_Count) <- Developmental_Class_Study_Count$Class

Developmental_Class_Species <- Developmental_Class_Data %>% select("phylo") %>% unique()

Developmental_Class_A_cor <- as.data.frame(A_cor)
Developmental_Class_A_cor <- Developmental_Class_A_cor[c(Developmental_Class_Species$phylo), c(Developmental_Class_Species$phylo)]
Developmental_Class_A_cor <- as.matrix(Developmental_Class_A_cor)

Developmental_Class_VCV_InCVR <- make_VCV_matrix(Developmental_Class_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  6ish minutes
  if(run){
    Developmental_Class_Model_CVR <- metafor::rma.mv(InCVR, V = Developmental_Class_VCV_InCVR, test = "t", dfs = "contain",
                                                     mods = ~ Class - 1,
                                                     random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                   ~1|Shared_Animal_Number, ~1|Measurement), 
                                                     R = list(phylo=Developmental_Class_A_cor), data = Developmental_Class_Data, method = "REML", sparse = TRUE, 
                                                     control=list(rel.tol=1e-9))
    saveRDS(Developmental_Class_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Developmental_Class_Model_CVR.rds")
  } else {
            Developmental_Class_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Developmental_Class_Model_CVR.rds")})

Developmental_Class_Model_CVR_rob <- robust(Developmental_Class_Model_CVR, cluster = Developmental_Class_Data$Study_ID, adjust = TRUE)

Developmental_Class_Model_CVR_Estimates <- data.frame(Class = substr(row.names(Developmental_Class_Model_CVR$b), 6, 100),
                                                      estimate = Developmental_Class_Model_CVR$b, 
                                                      ci.lb = Developmental_Class_Model_CVR$ci.lb, 
                                                      ci.ub = Developmental_Class_Model_CVR$ci.ub)
rownames(Developmental_Class_Model_CVR_Estimates) <- Developmental_Class_Model_CVR_Estimates$Class
Developmental_Class_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Class_Model_CVR), 2))

# Preparing Graph - Combined

developmental_class_rnames <- c("Actinopteri", "Amphibia", "Arachnida", 
                                "Branchiopoda", "Insecta")

developmental_class_k <- data.frame("k" = c(Developmental_Class_Exploration["Actinopteri", "Freq"], 
                                            Developmental_Class_Exploration["Amphibia", "Freq"], 
                                            Developmental_Class_Exploration["Arachnida", "Freq"], 
                                            Developmental_Class_Exploration["Branchiopoda", "Freq"], 
                                            Developmental_Class_Exploration["Insecta", "Freq"]), 
                                    row.names = developmental_class_rnames)

developmental_class_group_no <- data.frame("Spp No." = c(Developmental_Class_Species_Count["Actinopteri", "Freq"], 
                                                         Developmental_Class_Species_Count["Amphibia", "Freq"],  
                                                         Developmental_Class_Species_Count["Arachnida", "Freq"],
                                                         Developmental_Class_Species_Count["Branchiopoda", "Freq"],
                                                         Developmental_Class_Species_Count["Insecta", "Freq"]), 
                                           row.names = developmental_class_rnames)

developmental_class_study <- data.frame("Study" = c(Developmental_Class_Study_Count["Actinopteri", "Freq"], 
                                                    Developmental_Class_Study_Count["Amphibia", "Freq"],  
                                                    Developmental_Class_Study_Count["Arachnida", "Freq"],
                                                    Developmental_Class_Study_Count["Branchiopoda", "Freq"],
                                                    Developmental_Class_Study_Count["Insecta", "Freq"]), 
                                           row.names = developmental_class_rnames)

developmental_class_table <- data.frame(estimate = Developmental_Class_Model_CVR_Estimates[,"estimate"], 
                                        lowerCL = Developmental_Class_Model_CVR_Estimates[,"ci.lb"], 
                                        upperCL = Developmental_Class_Model_CVR_Estimates[,"ci.ub"], 
                                        K = developmental_class_k[,1], 
                                        group_no = developmental_class_group_no[,1], 
                                        row.names = developmental_class_rnames)
developmental_class_table$name <- row.names(developmental_class_table)

developmental_class_raw_mean <- c(unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Actinopteri") %>% 
                                                  select("InCVR"))), 
                                  unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Amphibia") %>% 
                                                  select("InCVR"))), 
                                  unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                                  select("InCVR"))), 
                                  unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Branchiopoda") %>% 
                                                  select("InCVR"))),
                                  unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                                  select("InCVR"))))

developmental_class_raw_name <- c(replicate(45, "Actinopteri"), 
                                  replicate(61, "Amphibia"), 
                                  replicate(102, "Arachnida"), 
                                  replicate(21, "Branchiopoda"),
                                  replicate(511, "Insecta"))

developmental_class_raw_df <- data.frame("Model" = developmental_class_raw_name, 
                                         "Effect" = developmental_class_raw_mean)

# Graph code - Combined

Developmental_Class_Order <- c("Insecta", "Branchiopoda","Arachnida", "Amphibia", "Actinopteri")

density_developmental_class_CVR <- developmental_class_table %>% mutate(name = fct_relevel(name, Developmental_Class_Order)) %>%
                                   ggplot() +
                                   geom_density_ridges(data = developmental_class_raw_df %>% mutate(Model = fct_relevel(Model, Developmental_Class_Order)), 
                                                       aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                       scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                   geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                  size = 1) +
                                   geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_class_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                   size = 1, fatten = 2) +
                                   theme_bw() +
                                   guides(fill = "none", colour = "none") +
                                   labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                   theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                    vjust = c(-2.7, -2.7, -2.7, -2.7, -2.7))) +
                                   theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                   theme(axis.ticks = element_blank()) +
                                   theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                   scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                   scale_colour_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                   scale_fill_manual(values = c("#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                   coord_cartesian(xlim = c(-1, 1)) +
                                   annotate('text',  x = 1, y = (seq(1, dim(developmental_class_table)[1], 1)+0.4),
                                   label= paste("italic(k)==", c(developmental_class_table["Insecta", "K"], 
                                                                 developmental_class_table["Branchiopoda", "K"],
                                                                 developmental_class_table["Arachnida", "K"], 
                                                                 developmental_class_table["Amphibia", "K"],
                                                                 developmental_class_table["Actinopteri", "K"]), "~","(", 
                                                               c(developmental_class_table["Insecta", "group_no"],
                                                                 developmental_class_table["Branchiopoda", "group_no"],
                                                                 developmental_class_table["Arachnida", "group_no"], 
                                                                 developmental_class_table["Amphibia", "group_no"],
                                                                 developmental_class_table["Actinopteri", "group_no"]), 
                                                ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                   geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Class_Model_CVR_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Developmental_Class_Model_CVR_Estimates["Branchiopoda", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Developmental_Class_Model_CVR_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                          paste(format(round(mean(exp(Developmental_Class_Model_CVR_Estimates["Amphibia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                          paste(format(round(mean(exp(Developmental_Class_Model_CVR_Estimates["Actinopteri", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                              x = -0.75, y = (seq(1, dim(developmental_class_table)[1], 1)+0.4)), size = 3.5)

density_developmental_class_CVR #(400x480)

# Preparing Graph - Part 1

developmental_class_rnames_1 <- c("Actinopteri", "Amphibia", "Arachnida")

developmental_class_k_1 <- data.frame("k" = c(Developmental_Class_Exploration["Actinopteri", "Freq"], 
                                              Developmental_Class_Exploration["Amphibia", "Freq"], 
                                              Developmental_Class_Exploration["Arachnida", "Freq"]), 
                                      row.names = developmental_class_rnames_1)

developmental_class_group_no_1 <- data.frame("Spp No." = c(Developmental_Class_Species_Count["Actinopteri", "Freq"], 
                                                           Developmental_Class_Species_Count["Amphibia", "Freq"],  
                                                           Developmental_Class_Species_Count["Arachnida", "Freq"]), 
                                             row.names = developmental_class_rnames_1)

developmental_class_study_1 <- data.frame("Study" = c(Developmental_Class_Study_Count["Actinopteri", "Freq"], 
                                                      Developmental_Class_Study_Count["Amphibia", "Freq"],  
                                                      Developmental_Class_Study_Count["Arachnida", "Freq"]), 
                                          row.names = developmental_class_rnames_1)

Developmental_Class_Model_CVR_Estimates_Reorder_1 <- Developmental_Class_Model_CVR_Estimates[c("Actinopteri", "Amphibia", "Arachnida"), ]

developmental_class_table_1 <- data.frame(estimate = Developmental_Class_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                          lowerCL = Developmental_Class_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                          upperCL = Developmental_Class_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                          K = developmental_class_k_1[,1], 
                                          group_no = developmental_class_group_no_1[,1], 
                                          row.names = developmental_class_rnames_1)
developmental_class_table_1$name <- row.names(developmental_class_table_1)

developmental_class_raw_mean_1 <- c(unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Actinopteri") %>% 
                                                  select("InCVR"))), 
                                    unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Amphibia") %>% 
                                                  select("InCVR"))), 
                                    unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Arachnida") %>% 
                                                  select("InCVR"))))

developmental_class_raw_name_1 <- c(replicate(45, "Actinopteri"), 
                                    replicate(61, "Amphibia"), 
                                    replicate(102, "Arachnida"))

developmental_class_raw_df_1 <- data.frame("Model" = developmental_class_raw_name_1, 
                                           "Effect" = developmental_class_raw_mean_1)

# Graph code - Part 1

Developmental_Class_Order_1 <- c("Arachnida", "Amphibia", "Actinopteri")

density_developmental_class_CVR_1 <- developmental_class_table_1 %>% mutate(name = fct_relevel(name, Developmental_Class_Order_1)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = developmental_class_raw_df_1 %>% mutate(Model = fct_relevel(Model, Developmental_Class_Order_1)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_class_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                     size = 1, fatten = 2) +
                                     theme_bw() +
                                     guides(fill = "none", colour = "none") +
                                     labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                           vjust = c(-2.7, -2.7, -2.7))) +
                                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                     theme(axis.ticks = element_blank()) +
                                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                     scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                     scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B")) +
                                     coord_cartesian(xlim = c(-1, 1)) +
                                     annotate('text',  x = 1, y = (seq(1, dim(developmental_class_table_1)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(developmental_class_table_1["Arachnida", "K"], 
                                                                   developmental_class_table_1["Amphibia", "K"],
                                                                   developmental_class_table_1["Actinopteri", "K"]), "~","(", 
                                                                 c(developmental_class_table_1["Arachnida", "group_no"], 
                                                                   developmental_class_table_1["Amphibia", "group_no"],
                                                                   developmental_class_table_1["Actinopteri", "group_no"]), 
                                                  ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Class_Model_CVR_Estimates["Arachnida", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                            paste(format(round(mean(exp(Developmental_Class_Model_CVR_Estimates["Amphibia", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                            paste(format(round(mean(exp(Developmental_Class_Model_CVR_Estimates["Actinopteri", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.75, y = (seq(1, dim(developmental_class_table_1)[1], 1)+0.4)), size = 3.5)

density_developmental_class_CVR_1 #(400x320)

# Preparing Graph - Part 2

developmental_class_rnames_2 <- c("Branchiopoda", "Insecta")

developmental_class_k_2 <- data.frame("k" = c(Developmental_Class_Exploration["Branchiopoda", "Freq"], 
                                              Developmental_Class_Exploration["Insecta", "Freq"]), 
                                      row.names = developmental_class_rnames_2)

developmental_class_group_no_2 <- data.frame("Spp No." = c(Developmental_Class_Species_Count["Branchiopoda", "Freq"],
                                                           Developmental_Class_Species_Count["Insecta", "Freq"]), 
                                             row.names = developmental_class_rnames_2)

developmental_class_study_2 <- data.frame("Study" = c(Developmental_Class_Study_Count["Branchiopoda", "Freq"],
                                                      Developmental_Class_Study_Count["Insecta", "Freq"]), 
                                          row.names = developmental_class_rnames_2)

Developmental_Class_Model_CVR_Estimates_Reorder_2 <- Developmental_Class_Model_CVR_Estimates[c("Branchiopoda", "Insecta"), ]

developmental_class_table_2 <- data.frame(estimate = Developmental_Class_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                          lowerCL = Developmental_Class_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                          upperCL = Developmental_Class_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                          K = developmental_class_k_2[,1], 
                                          group_no = developmental_class_group_no_2[,1], 
                                          row.names = developmental_class_rnames_2)
developmental_class_table_2$name <- row.names(developmental_class_table_2)

developmental_class_raw_mean_2 <- c(unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Branchiopoda") %>% 
                                                  select("InCVR"))),
                                    unlist(unname(Developmental_Class_Data %>% filter(`Class` == "Insecta") %>% 
                                                  select("InCVR"))))

developmental_class_raw_name_2 <- c(replicate(21, "Branchiopoda"),
                                    replicate(511, "Insecta"))

developmental_class_raw_df_2 <- data.frame("Model" = developmental_class_raw_name_2, 
                                           "Effect" = developmental_class_raw_mean_2)

# Graph code - Part 2

Developmental_Class_Order_2 <- c("Insecta", "Branchiopoda")

density_developmental_class_CVR_2 <- developmental_class_table_2 %>% mutate(name = fct_relevel(name, Developmental_Class_Order_2)) %>%
                                     ggplot() +
                                     geom_density_ridges(data = developmental_class_raw_df_2 %>% mutate(Model = fct_relevel(Model, Developmental_Class_Order_2)), 
                                                         aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                             scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                     geom_linerange(aes(y = rev(seq(1, dim(developmental_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                    size = 1) +
                                     geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_class_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                     size = 1, fatten = 2) +
                                     theme_bw() +
                                     guides(fill = "none", colour = "none") +
                                     labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                     theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                           vjust = c(-2.7, -2.7))) +
                                     theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                     theme(axis.ticks = element_blank()) +
                                     theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                     scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                     scale_colour_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                     scale_fill_manual(values = c("#5D7AA1", "#4A6E9C")) +
                                     coord_cartesian(xlim = c(-1, 1)) +
                                     annotate('text',  x = 1, y = (seq(1, dim(developmental_class_table_2)[1], 1)+0.4),
                                     label= paste("italic(k)==", c(developmental_class_table_2["Insecta", "K"], 
                                                                   developmental_class_table_2["Branchiopoda", "K"]), "~","(", 
                                                                 c(developmental_class_table_2["Insecta", "group_no"],
                                                                   developmental_class_table_2["Branchiopoda", "group_no"]), 
                                                 ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                     geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Class_Model_CVR_Estimates["Insecta", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                            paste(format(round(mean(exp(Developmental_Class_Model_CVR_Estimates["Branchiopoda", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                x = -0.75, y = (seq(1, dim(developmental_class_table_2)[1], 1)+0.4)), size = 3.5)

density_developmental_class_CVR_2 #(400x240)

##### Developmental Subset Model - Specific Trait Meta-Regression - CVR #####
Developmental_Specific_Trait_Exploration <- Developmental_Subset_Data %>% select("Measurement") %>% table() %>% data.frame()
Developmental_Specific_Trait_Exploration <- Developmental_Specific_Trait_Exploration %>% filter(Freq > 10)
rownames(Developmental_Specific_Trait_Exploration) <- Developmental_Specific_Trait_Exploration$Measurement

Developmental_Specific_Trait_Data <- Developmental_Subset_Data %>% filter(Measurement == "Development Time"| 
                                                                          Measurement == "Fecundity"|
                                                                          Measurement == "Food Consumption"|
                                                                          Measurement == "Head Width"|
                                                                          Measurement == "Length"|
                                                                          Measurement == "Locomotor Performance"|
                                                                          Measurement == "Longevity"|
                                                                          Measurement == "Mass"|
                                                                          Measurement == "Tail Length")

Developmental_Specific_Trait_Species_Count <- Developmental_Specific_Trait_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>% 
                                              filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Developmental_Specific_Trait_Species_Count) <- Developmental_Specific_Trait_Species_Count$Measurement

Developmental_Specific_Trait_Study_Count <- Developmental_Specific_Trait_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>% 
                                            filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Developmental_Specific_Trait_Study_Count) <- Developmental_Specific_Trait_Study_Count$Measurement

Developmental_Specific_Trait_Species <- Developmental_Specific_Trait_Data %>% select("phylo") %>% unique()

Developmental_Specific_Trait_A_cor <- as.data.frame(A_cor)
Developmental_Specific_Trait_A_cor <- Developmental_Specific_Trait_A_cor[c(Developmental_Specific_Trait_Species$phylo), c(Developmental_Specific_Trait_Species$phylo)]
Developmental_Specific_Trait_A_cor <- as.matrix(Developmental_Specific_Trait_A_cor)

Developmental_Specific_Trait_VCV_InCVR <- make_VCV_matrix(Developmental_Specific_Trait_Data, V = "v_InCVR", cluster = "Shared_Control_Number")

run <- FALSE
system.time( #  5ish minutes
  if(run){
    Developmental_Specific_Trait_Model_CVR <- metafor::rma.mv(InCVR, V = Developmental_Specific_Trait_VCV_InCVR, test = "t", dfs = "contain",
                                                              mods = ~ Measurement - 1,
                                                              random = list(~1|phylo, ~1|Study_ID, ~1|obs, ~1|Scientific_Name, 
                                                                            ~1|Shared_Animal_Number), 
                                                              R = list(phylo=Developmental_Specific_Trait_A_cor), data = Developmental_Specific_Trait_Data, method = "REML", sparse = TRUE, 
                                                              control=list(rel.tol=1e-9))
    saveRDS(Developmental_Specific_Trait_Model_CVR, "./3.Data_Analysis/2.Outputs/Models/Developmental_Specific_Trait_Model_CVR.rds")
  } else {
            Developmental_Specific_Trait_Model_CVR <- readRDS("./3.Data_Analysis/2.Outputs/Models/Developmental_Specific_Trait_Model_CVR.rds")})

Developmental_Specific_Trait_Model_CVR_rob <- robust(Developmental_Specific_Trait_Model_CVR, cluster = Developmental_Specific_Trait_Data$Study_ID, adjust = TRUE)

Developmental_Specific_Trait_Model_CVR_Estimates <- data.frame(Trait = substr(row.names(Developmental_Specific_Trait_Model_CVR$b), 12, 100),
                                                               estimate = Developmental_Specific_Trait_Model_CVR$b, 
                                                               ci.lb = Developmental_Specific_Trait_Model_CVR$ci.lb, 
                                                               ci.ub = Developmental_Specific_Trait_Model_CVR$ci.ub)
rownames(Developmental_Specific_Trait_Model_CVR_Estimates) <- Developmental_Specific_Trait_Model_CVR_Estimates$Trait
Developmental_Specific_Trait_Model_CVR_i2 <- data.frame(round(orchaRd::i2_ml(Developmental_Specific_Trait_Model_CVR), 2))

# Preparing Graph - Combined

developmental_specific_trait_rnames <- c("Development Time", "Fecundity", "Food Consumption", "Head Width", 
                                         "Length", "Locomotor Performance", "Longevity", "Mass", "Tail Length")

developmental_specific_trait_k <- data.frame("k" = c(Developmental_Specific_Trait_Exploration["Development Time", "Freq"], 
                                                     Developmental_Specific_Trait_Exploration["Fecundity", "Freq"], 
                                                     Developmental_Specific_Trait_Exploration["Food Consumption", "Freq"],
                                                     Developmental_Specific_Trait_Exploration["Head Width", "Freq"],
                                                     Developmental_Specific_Trait_Exploration["Length", "Freq"], 
                                                     Developmental_Specific_Trait_Exploration["Locomotor Performance", "Freq"], 
                                                     Developmental_Specific_Trait_Exploration["Longevity", "Freq"], 
                                                     Developmental_Specific_Trait_Exploration["Mass", "Freq"], 
                                                     Developmental_Specific_Trait_Exploration["Tail Length", "Freq"]), 
                                             row.names = developmental_specific_trait_rnames)

developmental_specific_trait_group_no <- data.frame("Spp No." = c(Developmental_Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                                  Developmental_Specific_Trait_Species_Count["Fecundity", "Freq"], 
                                                                  Developmental_Specific_Trait_Species_Count["Food Consumption", "Freq"],
                                                                  Developmental_Specific_Trait_Species_Count["Head Width", "Freq"],
                                                                  Developmental_Specific_Trait_Species_Count["Length", "Freq"], 
                                                                  Developmental_Specific_Trait_Species_Count["Locomotor Performance", "Freq"], 
                                                                  Developmental_Specific_Trait_Species_Count["Longevity", "Freq"], 
                                                                  Developmental_Specific_Trait_Species_Count["Mass", "Freq"],  
                                                                  Developmental_Specific_Trait_Species_Count["Tail Length", "Freq"]), 
                                                    row.names = developmental_specific_trait_rnames)

developmental_specific_trait_study <- data.frame("Study" = c(Developmental_Specific_Trait_Study_Count["Development Time", "Freq"], 
                                                             Developmental_Specific_Trait_Study_Count["Fecundity", "Freq"], 
                                                             Developmental_Specific_Trait_Study_Count["Food Consumption", "Freq"],
                                                             Developmental_Specific_Trait_Study_Count["Head Width", "Freq"],
                                                             Developmental_Specific_Trait_Study_Count["Length", "Freq"], 
                                                             Developmental_Specific_Trait_Study_Count["Locomotor Performance", "Freq"], 
                                                             Developmental_Specific_Trait_Study_Count["Longevity", "Freq"], 
                                                             Developmental_Specific_Trait_Study_Count["Mass", "Freq"],  
                                                             Developmental_Specific_Trait_Study_Count["Tail Length", "Freq"]), 
                                                    row.names = developmental_specific_trait_rnames)

developmental_specific_trait_table <- data.frame(estimate = Developmental_Specific_Trait_Model_CVR_Estimates[,"estimate"], 
                                                 lowerCL = Developmental_Specific_Trait_Model_CVR_Estimates[,"ci.lb"], 
                                                 upperCL = Developmental_Specific_Trait_Model_CVR_Estimates[,"ci.ub"], 
                                                 K = developmental_specific_trait_k[,1], 
                                                 group_no = developmental_specific_trait_group_no[,1], 
                                                 row.names = developmental_specific_trait_rnames)
developmental_specific_trait_table$name <- row.names(developmental_specific_trait_table)

developmental_specific_trait_raw_mean <- c(unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                                           select("InCVR"))), 
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Fecundity") %>% 
                                                           select("InCVR"))),
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Food Consumption") %>% 
                                                           select("InCVR"))), 
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Head Width") %>% 
                                                           select("InCVR"))),
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                                           select("InCVR"))), 
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Locomotor Performance") %>% 
                                                           select("InCVR"))), 
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Longevity") %>% 
                                                           select("InCVR"))), 
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                                           select("InCVR"))),  
                                           unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Tail Length") %>% 
                                                           select("InCVR"))))

developmental_specific_trait_raw_name <- c(replicate(287, "Development Time"), 
                                           replicate(67, "Fecundity"),
                                           replicate(12, "Food Consumption"),
                                           replicate(14, "Head Width"),
                                           replicate(89, "Length"), 
                                           replicate(27, "Locomotor Performance"), 
                                           replicate(91, "Longevity"), 
                                           replicate(116, "Mass"),  
                                           replicate(26, "Tail Length"))

developmental_specific_trait_raw_df <- data.frame("Model" = developmental_specific_trait_raw_name, 
                                                  "Effect" = developmental_specific_trait_raw_mean)

# Graph code - Combined

Developmental_Specific_Trait_Order <- c("Tail Length", "Mass", "Longevity", "Locomotor Performance", "Length", 
                                        "Head Width", "Food Consumption", "Fecundity", "Development Time")

density_developmental_specific_trait_CVR <- developmental_specific_trait_table %>% mutate(name = fct_relevel(name, Developmental_Specific_Trait_Order)) %>%
                                            ggplot() +
                                            geom_density_ridges(data = developmental_specific_trait_raw_df %>% mutate(Model = fct_relevel(Model, Developmental_Specific_Trait_Order)), 
                                                                aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                            geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                           size = 1) +
                                            geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_specific_trait_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                            size = 1, fatten = 2) +
                                            theme_bw() +
                                            guides(fill = "none", colour = "none") +
                                            labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                            theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                                             vjust = c(-2.7, -2.7, -2.7, -0.8, -2.7, 
                                                                                       -2.7, -0.8, -2.7, -0.8))) +
                                            theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                            theme(axis.ticks = element_blank()) +
                                            theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                            theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                            scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                            scale_colour_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692", "#3D5E89", 
                                                                           "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                            scale_fill_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692", "#3D5E89", 
                                                                         "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                            coord_cartesian(xlim = c(-3, 3)) +
                                            annotate('text',  x = 3, y = (seq(1, dim(developmental_specific_trait_table)[1], 1)+0.4),
                                            label= paste("italic(k)==", c(developmental_specific_trait_table["Tail Length", "K"],
                                                                          developmental_specific_trait_table["Mass", "K"],
                                                                          developmental_specific_trait_table["Longevity", "K"], 
                                                                          developmental_specific_trait_table["Locomotor Performance", "K"], 
                                                                          developmental_specific_trait_table["Length", "K"],
                                                                          developmental_specific_trait_table["Head Width", "K"],
                                                                          developmental_specific_trait_table["Food Consumption", "K"], 
                                                                          developmental_specific_trait_table["Fecundity", "K"],
                                                                          developmental_specific_trait_table["Development Time", "K"]), "~","(", 
                                                                        c(developmental_specific_trait_table["Tail Length", "group_no"], 
                                                                          developmental_specific_trait_table["Mass", "group_no"],
                                                                          developmental_specific_trait_table["Longevity", "group_no"], 
                                                                          developmental_specific_trait_table["Locomotor Performance", "group_no"], 
                                                                          developmental_specific_trait_table["Length", "group_no"], 
                                                                          developmental_specific_trait_table["Head Width", "group_no"],
                                                                          developmental_specific_trait_table["Food Consumption", "group_no"], 
                                                                          developmental_specific_trait_table["Fecundity", "group_no"],
                                                                          developmental_specific_trait_table["Development Time", "group_no"]), 
                                                         ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                            geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Tail Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Longevity", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                   paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Locomotor Performance", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                   paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Head Width", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Food Consumption", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                   paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Fecundity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                   paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                       x = -2.5, y = (seq(1, dim(developmental_specific_trait_table)[1], 1)+0.4)), size = 3.5)

density_developmental_specific_trait_CVR #(400x800)

# Preparing Graph - Part 1

developmental_specific_trait_rnames_1 <- c("Development Time", "Fecundity", "Food Consumption", "Head Width", "Length")

developmental_specific_trait_k_1 <- data.frame("k" = c(Developmental_Specific_Trait_Exploration["Development Time", "Freq"], 
                                                       Developmental_Specific_Trait_Exploration["Fecundity", "Freq"], 
                                                       Developmental_Specific_Trait_Exploration["Food Consumption", "Freq"],
                                                       Developmental_Specific_Trait_Exploration["Head Width", "Freq"],
                                                       Developmental_Specific_Trait_Exploration["Length", "Freq"]), 
                                               row.names = developmental_specific_trait_rnames_1)

developmental_specific_trait_group_no_1 <- data.frame("Spp No." = c(Developmental_Specific_Trait_Species_Count["Development Time", "Freq"], 
                                                                    Developmental_Specific_Trait_Species_Count["Fecundity", "Freq"], 
                                                                    Developmental_Specific_Trait_Species_Count["Food Consumption", "Freq"],
                                                                    Developmental_Specific_Trait_Species_Count["Head Width", "Freq"],
                                                                    Developmental_Specific_Trait_Species_Count["Length", "Freq"]), 
                                                      row.names = developmental_specific_trait_rnames_1)

developmental_specific_trait_study_1 <- data.frame("Study" = c(Developmental_Specific_Trait_Study_Count["Development Time", "Freq"], 
                                                               Developmental_Specific_Trait_Study_Count["Fecundity", "Freq"], 
                                                               Developmental_Specific_Trait_Study_Count["Food Consumption", "Freq"],
                                                               Developmental_Specific_Trait_Study_Count["Head Width", "Freq"],
                                                               Developmental_Specific_Trait_Study_Count["Length", "Freq"]), 
                                                   row.names = developmental_specific_trait_rnames_1)

Developmental_Specific_Trait_Model_CVR_Estimates_Reorder_1 <- Developmental_Specific_Trait_Model_CVR_Estimates[c("Development Time", "Fecundity", "Food Consumption", "Head Width", "Length"), ]

developmental_specific_trait_table_1 <- data.frame(estimate = Developmental_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"estimate"], 
                                                   lowerCL = Developmental_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"ci.lb"], 
                                                   upperCL = Developmental_Specific_Trait_Model_CVR_Estimates_Reorder_1[,"ci.ub"], 
                                                   K = developmental_specific_trait_k_1[,1], 
                                                   group_no = developmental_specific_trait_group_no_1[,1], 
                                                   row.names = developmental_specific_trait_rnames_1)
developmental_specific_trait_table_1$name <- row.names(developmental_specific_trait_table_1)

developmental_specific_trait_raw_mean_1 <- c(unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Development Time") %>% 
                                                           select("InCVR"))), 
                                             unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Fecundity") %>% 
                                                           select("InCVR"))),
                                             unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Food Consumption") %>% 
                                                           select("InCVR"))), 
                                             unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Head Width") %>% 
                                                           select("InCVR"))),
                                             unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Length") %>% 
                                                           select("InCVR"))))

developmental_specific_trait_raw_name_1 <- c(replicate(287, "Development Time"), 
                                             replicate(67, "Fecundity"),
                                             replicate(12, "Food Consumption"),
                                             replicate(14, "Head Width"),
                                             replicate(89, "Length"))

developmental_specific_trait_raw_df_1 <- data.frame("Model" = developmental_specific_trait_raw_name_1, 
                                                    "Effect" = developmental_specific_trait_raw_mean_1)

# Graph code - Part 1

Developmental_Specific_Trait_Order_1 <- c("Length", "Head Width", "Food Consumption", "Fecundity", "Development Time")

density_developmental_specific_trait_CVR_1 <- developmental_specific_trait_table_1 %>% mutate(name = fct_relevel(name, Developmental_Specific_Trait_Order_1)) %>%
                                              ggplot() +
                                              geom_density_ridges(data = developmental_specific_trait_raw_df_1 %>% mutate(Model = fct_relevel(Model, Developmental_Specific_Trait_Order_1)), 
                                                                  aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                  scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                              geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                             size = 1) +
                                              geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_specific_trait_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                              size = 1, fatten = 2) +
                                              theme_bw() +
                                              guides(fill = "none", colour = "none") +
                                              labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                              theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                    vjust = c(-2.7, -2.7, -0.8, -2.7, -0.8))) +
                                              theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                              theme(axis.ticks = element_blank()) +
                                              theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                              theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                              scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                              scale_colour_manual(values = c("#3D5E89", "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                              scale_fill_manual(values = c("#3D5E89", "#234373", "#1C375F", "#0D2A51", "#0F2643")) +
                                              coord_cartesian(xlim = c(-3, 3)) +
                                              annotate('text',  x = 3, y = (seq(1, dim(developmental_specific_trait_table_1)[1], 1)+0.4),
                                              label= paste("italic(k)==", c(developmental_specific_trait_table_1["Length", "K"],
                                                                            developmental_specific_trait_table_1["Head Width", "K"],
                                                                            developmental_specific_trait_table_1["Food Consumption", "K"], 
                                                                            developmental_specific_trait_table_1["Fecundity", "K"],
                                                                            developmental_specific_trait_table_1["Development Time", "K"]), "~","(", 
                                                                          c(developmental_specific_trait_table_1["Length", "group_no"], 
                                                                            developmental_specific_trait_table_1["Head Width", "group_no"],
                                                                            developmental_specific_trait_table_1["Food Consumption", "group_no"], 
                                                                            developmental_specific_trait_table_1["Fecundity", "group_no"],
                                                                            developmental_specific_trait_table_1["Development Time", "group_no"]), 
                                                           ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                              geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                     paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Head Width", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                     paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Food Consumption", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                     paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Fecundity", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                     paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Development Time", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                         x = -2.5, y = (seq(1, dim(developmental_specific_trait_table_1)[1], 1)+0.4)), size = 3.5)

density_developmental_specific_trait_CVR_1 #(400x480)

# Preparing Graph - Part 2

developmental_specific_trait_rnames_2 <- c("Locomotor Performance", "Longevity", "Mass", "Tail Length")

developmental_specific_trait_k_2 <- data.frame("k" = c(Developmental_Specific_Trait_Exploration["Locomotor Performance", "Freq"], 
                                                       Developmental_Specific_Trait_Exploration["Longevity", "Freq"], 
                                                       Developmental_Specific_Trait_Exploration["Mass", "Freq"], 
                                                       Developmental_Specific_Trait_Exploration["Tail Length", "Freq"]), 
                                               row.names = developmental_specific_trait_rnames_2)

developmental_specific_trait_group_no_2 <- data.frame("Spp No." = c(Developmental_Specific_Trait_Species_Count["Locomotor Performance", "Freq"], 
                                                                    Developmental_Specific_Trait_Species_Count["Longevity", "Freq"], 
                                                                    Developmental_Specific_Trait_Species_Count["Mass", "Freq"],  
                                                                    Developmental_Specific_Trait_Species_Count["Tail Length", "Freq"]), 
                                                      row.names = developmental_specific_trait_rnames_2)

developmental_specific_trait_study_2 <- data.frame("Study" = c(Developmental_Specific_Trait_Study_Count["Locomotor Performance", "Freq"], 
                                                               Developmental_Specific_Trait_Study_Count["Longevity", "Freq"], 
                                                               Developmental_Specific_Trait_Study_Count["Mass", "Freq"],  
                                                               Developmental_Specific_Trait_Study_Count["Tail Length", "Freq"]), 
                                                   row.names = developmental_specific_trait_rnames_2)

Developmental_Specific_Trait_Model_CVR_Estimates_Reorder_2 <- Developmental_Specific_Trait_Model_CVR_Estimates[c("Locomotor Performance", "Longevity", "Mass", "Tail Length"), ]

developmental_specific_trait_table_2 <- data.frame(estimate = Developmental_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"estimate"], 
                                                   lowerCL = Developmental_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"ci.lb"], 
                                                   upperCL = Developmental_Specific_Trait_Model_CVR_Estimates_Reorder_2[,"ci.ub"], 
                                                   K = developmental_specific_trait_k_2[,1], 
                                                   group_no = developmental_specific_trait_group_no_2[,1], 
                                                   row.names = developmental_specific_trait_rnames_2)
developmental_specific_trait_table_2$name <- row.names(developmental_specific_trait_table_2)

developmental_specific_trait_raw_mean_2 <- c(unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Locomotor Performance") %>% 
                                                           select("InCVR"))), 
                                             unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Longevity") %>% 
                                                           select("InCVR"))), 
                                             unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Mass") %>% 
                                                           select("InCVR"))),  
                                             unlist(unname(Developmental_Specific_Trait_Data %>% filter(`Measurement` == "Tail Length") %>% 
                                                           select("InCVR"))))

developmental_specific_trait_raw_name_2 <- c(replicate(27, "Locomotor Performance"), 
                                             replicate(91, "Longevity"), 
                                             replicate(116, "Mass"),  
                                             replicate(26, "Tail Length"))

developmental_specific_trait_raw_df_2 <- data.frame("Model" = developmental_specific_trait_raw_name_2, 
                                                    "Effect" = developmental_specific_trait_raw_mean_2)

# Graph code - Part 2

Developmental_Specific_Trait_Order_2 <- c("Tail Length", "Mass", "Longevity", "Locomotor Performance")

density_developmental_specific_trait_CVR_2 <- developmental_specific_trait_table_2 %>% mutate(name = fct_relevel(name, Developmental_Specific_Trait_Order_2)) %>%
                                              ggplot() +
                                              geom_density_ridges(data = developmental_specific_trait_raw_df_2 %>% mutate(Model = fct_relevel(Model, Developmental_Specific_Trait_Order_2)), 
                                                                  aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                                  scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                                              geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                                             size = 1) +
                                              geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table_2)[1], 1)), xmin = min(developmental_specific_trait_raw_df_2$Effect)-0.1, xmax = -3.5, colour = name),
                                                             size = 1) +
                                              geom_linerange(aes(y = rev(seq(1, dim(developmental_specific_trait_table_2)[1], 1)), xmin = max(developmental_specific_trait_raw_df_2$Effect)+0.15, xmax = 3.5, colour = name),
                                                             size = 1) +
                                              geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(developmental_specific_trait_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                                              size = 1, fatten = 2) +
                                              theme_bw() +
                                              guides(fill = "none", colour = "none") +
                                              labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                                              theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                    vjust = c(-2.7, -2.7, -2.7, -0.8))) +
                                              theme(axis.text.x = element_text(margin = margin(b = 5))) +
                                              theme(axis.ticks = element_blank()) +
                                              theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                              theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                                              scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                                              scale_colour_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692")) +
                                              scale_fill_manual(values = c("#6582A9", "#607C9F", "#4A6E9C", "#446692")) +
                                              coord_cartesian(xlim = c(-3, 3)) +
                                              annotate('text',  x = 3, y = (seq(1, dim(developmental_specific_trait_table_2)[1], 1)+0.4),
                                              label= paste("italic(k)==", c(developmental_specific_trait_table_2["Tail Length", "K"],
                                                                            developmental_specific_trait_table_2["Mass", "K"],
                                                                            developmental_specific_trait_table_2["Longevity", "K"], 
                                                                            developmental_specific_trait_table_2["Locomotor Performance", "K"]), "~","(", 
                                                                          c(developmental_specific_trait_table_2["Tail Length", "group_no"], 
                                                                            developmental_specific_trait_table_2["Mass", "group_no"],
                                                                            developmental_specific_trait_table_2["Longevity", "group_no"], 
                                                                            developmental_specific_trait_table_2["Locomotor Performance", "group_no"]), 
                                                           ")"), parse = TRUE, hjust = "right", size = 3.5) +
                                              geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Tail Length", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                     paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Mass", "estimate"])-1)*100, 2), nsmall = 2), "%"),
                                                                     paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Longevity", "estimate"])-1)*100, 2), nsmall = 2), "%"), 
                                                                     paste(format(round(mean(exp(Developmental_Specific_Trait_Model_CVR_Estimates["Locomotor Performance", "estimate"])-1)*100, 2), nsmall = 2), "%")), 
                                                         x = -2.5, y = (seq(1, dim(developmental_specific_trait_table_2)[1], 1)+0.4)), size = 3.5)

density_developmental_specific_trait_CVR_2 #(400x400)

##### Summary of Developmental Plots #####

Developmental_Layout <- rbind(c(1, 3), 
                              c(1, 3),
                              c(1, 3),
                              c(1, 3),
                              c(1, 3),
                              c(1, 3),
                              c(1, 4),
                              c(1, 4),
                              c(1, 4),
                              c(2, 4),
                              c(2, 5),
                              c(2, 5),
                              c(2, 5),
                              c(2, 5))

Developmental_Combined_CVR <- grid.arrange(density_developmental_specific_trait_CVR, density_developmental_class_CVR, density_developmental_trait_CVR, 
                                           density_developmental_fluctuation_CVR, density_developmental_exposure_CVR, 
                                           layout_matrix = Developmental_Layout)

Developmental_Combined_CVR #(850 x 1500 - does not include amplitude plot)

##### Meta-analytic Models (Intercept Only) - Subset Graph #####

# Preparing Data - Combined

intercept_rnames <- c("Overall", "Population-level Traits", "Individual-level Traits", "Aquatic", 
                      "Terrestrial", "Acclimation", "Developmental")

intercept_means_df <- data.frame("Mean" = c(Overall_Model_CVR_Estimates$estimate, 
                                            Population_Model_CVR_Estimates$estimate, 
                                            Individual_Model_CVR_Estimates$estimate, 
                                            Aquatic_Model_CVR_Estimates$estimate, 
                                            Terrestrial_Model_CVR_Estimates$estimate, 
                                            Acclimation_Model_CVR_Estimates$estimate, 
                                            Developmental_Model_CVR_Estimates$estimate), 
                                 row.names = intercept_rnames)

intercept_low_df <- data.frame("Low CI" = c(Overall_Model_CVR_Estimates$ci.lb, 
                                            Population_Model_CVR_Estimates$ci.lb, 
                                            Individual_Model_CVR_Estimates$ci.lb, 
                                            Aquatic_Model_CVR_Estimates$ci.lb, 
                                            Terrestrial_Model_CVR_Estimates$ci.lb, 
                                            Acclimation_Model_CVR_Estimates$ci.lb, 
                                            Developmental_Model_CVR_Estimates$ci.lb), 
                               row.names = intercept_rnames)

intercept_high_df <- data.frame("High CI" = c(Overall_Model_CVR_Estimates$ci.ub, 
                                              Population_Model_CVR_Estimates$ci.ub, 
                                              Individual_Model_CVR_Estimates$ci.ub, 
                                              Aquatic_Model_CVR_Estimates$ci.ub, 
                                              Terrestrial_Model_CVR_Estimates$ci.ub, 
                                              Acclimation_Model_CVR_Estimates$ci.ub, 
                                              Developmental_Model_CVR_Estimates$ci.ub), 
                                row.names = intercept_rnames)

intercept_k <- data.frame("k" = c(length(data$Effect_Size_ID), 
                                  length(Population_Subset_Data$Effect_Size_ID), 
                                  length(Individual_Subset_Data$Effect_Size_ID), 
                                  length(Aquatic_Subset_Data$Effect_Size_ID), 
                                  length(Terrestrial_Subset_Data$Effect_Size_ID), 
                                  length(Acclimation_Subset_Data$Effect_Size_ID),
                                  length(Developmental_Subset_Data$Effect_Size_ID)), 
                          row.names = intercept_rnames)

intercept_group_no <- data.frame("Spp No." = c(length(unique(data$Scientific_Name)), 
                                               length(unique(Population_Subset_Data$Scientific_Name)), 
                                               length(unique(Individual_Subset_Data$Scientific_Name)), 
                                               length(unique(Aquatic_Subset_Data$Scientific_Name)), 
                                               length(unique(Terrestrial_Subset_Data$Scientific_Name)), 
                                               length(unique(Acclimation_Subset_Data$Scientific_Name)),
                                               length(unique(Developmental_Subset_Data$Scientific_Name))),
                                 row.names = intercept_rnames)

intercept_study <- data.frame("Study" = c(length(unique(data$Study_ID)), 
                                          length(unique(Population_Subset_Data$Study_ID)), 
                                          length(unique(Individual_Subset_Data$Study_ID)), 
                                          length(unique(Aquatic_Subset_Data$Study_ID)), 
                                          length(unique(Terrestrial_Subset_Data$Study_ID)), 
                                          length(unique(Acclimation_Subset_Data$Study_ID)),
                                          length(unique(Developmental_Subset_Data$Study_ID))),
                              row.names = intercept_rnames)

intercept_table <- data.frame(estimate = intercept_means_df[,1], 
                              lowerCL = intercept_low_df[,1], 
                              upperCL = intercept_high_df[,1], 
                              K = intercept_k[,1], 
                              group_no = intercept_group_no[,1], 
                              row.names = intercept_rnames)
intercept_table$name <- row.names(intercept_table)

intercept_raw_mean <- c(data$InCVR, 
                        Population_Subset_Data$InCVR, 
                        Individual_Subset_Data$InCVR, 
                        Aquatic_Subset_Data$InCVR, 
                        Terrestrial_Subset_Data$InCVR, 
                        Acclimation_Subset_Data$InCVR, 
                        Developmental_Subset_Data$InCVR)

intercept_raw_name <- c(replicate(1492, "Overall"), 
                        replicate(168, "Population-level Traits"), 
                        replicate(1324, "Individual-level Traits"), 
                        replicate(395, "Aquatic"), 
                        replicate(929, "Terrestrial"), 
                        replicate(336, "Acclimation"), 
                        replicate(988, "Developmental"))

intercept_raw_df <- data.frame("Model" = intercept_raw_name, 
                               "Effect" = intercept_raw_mean)
# Graph code - Combined

Intercept_Order <- c("Developmental", "Acclimation", "Terrestrial", "Aquatic", 
                     "Individual-level Traits", "Population-level Traits", "Overall")

density_intercept_CVR <- intercept_table %>% mutate(name = fct_relevel(name, Intercept_Order)) %>%
                         ggplot() +
                         geom_density_ridges(data = intercept_raw_df %>% mutate(Model = fct_relevel(Model, Intercept_Order)), 
                                             aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                             scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                         geom_linerange(aes(y = rev(seq(1, dim(intercept_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                        size = 1) +
                         geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(intercept_table)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                         size = 1, fatten = 2) +
                         theme_bw() +
                         guides(fill = "none", colour = "none") +
                         labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                         theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                               vjust = c(-2.7, -2.7, -2.7, -2.7, -0.8, -0.8, -2.7))) +
                         theme(axis.text.x = element_text(margin = margin(b = 5))) +
                         theme(axis.ticks = element_blank()) +
                         theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                         theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                         scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                         scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                         scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C", "#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                         coord_cartesian(xlim = c(-1, 1)) +
                         annotate('text',  x = 1, y = (seq(1, dim(intercept_table)[1], 1)+0.4),
                         label= paste("italic(k)==", c(intercept_table["Developmental", "K"], 
                                                       intercept_table["Acclimation", "K"], 
                                                       intercept_table["Terrestrial", "K"], 
                                                       intercept_table["Aquatic", "K"], 
                                                       intercept_table["Individual-level Traits", "K"], 
                                                       intercept_table["Population-level Traits", "K"], 
                                                       intercept_table["Overall", "K"]), "~","(", 
                                                     c(intercept_table["Developmental", "group_no"], 
                                                       intercept_table["Acclimation", "group_no"], 
                                                       intercept_table["Terrestrial", "group_no"], 
                                                       intercept_table["Aquatic", "group_no"], 
                                                       intercept_table["Individual-level Traits", "group_no"], 
                                                       intercept_table["Population-level Traits", "group_no"], 
                                                       intercept_table["Overall", "group_no"]), 
                                     ")"), parse = TRUE, hjust = "right", size = 3.5) +
                         geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"), 
                                                paste(format(round(mean(exp(Acclimation_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"),
                                                paste(format(round(mean(exp(Terrestrial_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"), 
                                                paste(format(round(mean(exp(Aquatic_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"),
                                                paste(format(round(mean(exp(Individual_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"), 
                                                paste(format(round(mean(exp(Population_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"),
                                                paste(format(round(mean(exp(Overall_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%")), 
                                    x = -0.75, y = (seq(1, dim(intercept_table)[1], 1)+0.4)), size = 3.5)

density_intercept_CVR #(400x400)

# Preparing Data - Part 1

intercept_rnames_1 <- c("Overall", "Population-level Traits", "Individual-level Traits", "Aquatic")

intercept_means_df_1 <- data.frame("Mean" = c(Overall_Model_CVR_Estimates$estimate, 
                                              Population_Model_CVR_Estimates$estimate, 
                                              Individual_Model_CVR_Estimates$estimate, 
                                              Aquatic_Model_CVR_Estimates$estimate), 
                                 row.names = intercept_rnames_1)

intercept_low_df_1 <- data.frame("Low CI" = c(Overall_Model_CVR_Estimates$ci.lb, 
                                              Population_Model_CVR_Estimates$ci.lb, 
                                              Individual_Model_CVR_Estimates$ci.lb, 
                                              Aquatic_Model_CVR_Estimates$ci.lb), 
                               row.names = intercept_rnames_1)

intercept_high_df_1 <- data.frame("High CI" = c(Overall_Model_CVR_Estimates$ci.ub, 
                                                Population_Model_CVR_Estimates$ci.ub, 
                                                Individual_Model_CVR_Estimates$ci.ub, 
                                                Aquatic_Model_CVR_Estimates$ci.ub), 
                                row.names = intercept_rnames_1)

intercept_k_1 <- data.frame("k" = c(length(data$Effect_Size_ID), 
                                    length(Population_Subset_Data$Effect_Size_ID), 
                                    length(Individual_Subset_Data$Effect_Size_ID), 
                                    length(Aquatic_Subset_Data$Effect_Size_ID)), 
                            row.names = intercept_rnames_1)

intercept_group_no_1 <- data.frame("Spp No." = c(length(unique(data$Scientific_Name)), 
                                                 length(unique(Population_Subset_Data$Scientific_Name)), 
                                                 length(unique(Individual_Subset_Data$Scientific_Name)), 
                                                 length(unique(Aquatic_Subset_Data$Scientific_Name))),
                                   row.names = intercept_rnames_1)

intercept_study_1 <- data.frame("Study" = c(length(unique(data$Study_ID)), 
                                            length(unique(Population_Subset_Data$Study_ID)), 
                                            length(unique(Individual_Subset_Data$Study_ID)), 
                                            length(unique(Aquatic_Subset_Data$Study_ID))),
                                   row.names = intercept_rnames_1)

intercept_table_1 <- data.frame(estimate = intercept_means_df_1[,1], 
                                lowerCL = intercept_low_df_1[,1], 
                                upperCL = intercept_high_df_1[,1], 
                                K = intercept_k_1[,1], 
                                group_no = intercept_group_no_1[,1], 
                                row.names = intercept_rnames_1)
intercept_table_1$name <- row.names(intercept_table_1)

intercept_raw_mean_1 <- c(data$InCVR, 
                          Population_Subset_Data$InCVR, 
                          Individual_Subset_Data$InCVR, 
                          Aquatic_Subset_Data$InCVR)

intercept_raw_name_1 <- c(replicate(1492, "Overall"), 
                          replicate(168, "Population-level Traits"), 
                          replicate(1324, "Individual-level Traits"), 
                          replicate(395, "Aquatic"))

intercept_raw_df_1 <- data.frame("Model" = intercept_raw_name_1, 
                                 "Effect" = intercept_raw_mean_1)
# Graph code - Part 1

Intercept_Order_1 <- c("Aquatic", "Individual-level Traits", 
                       "Population-level Traits", "Overall")

density_intercept_CVR_1 <- intercept_table_1 %>% mutate(name = fct_relevel(name, Intercept_Order_1)) %>%
                           ggplot() +
                           geom_density_ridges(data = intercept_raw_df_1 %>% mutate(Model = fct_relevel(Model, Intercept_Order_1)), 
                                               aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                               scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                           geom_linerange(aes(y = rev(seq(1, dim(intercept_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                          size = 1) +
                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(intercept_table_1)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                           size = 1, fatten = 2) +
                           theme_bw() +
                           guides(fill = "none", colour = "none") +
                           labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                           theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                                            vjust = c(-2.7, -0.8, -0.8, -2.7))) +
                           theme(axis.text.x = element_text(margin = margin(b = 5))) +
                           theme(axis.ticks = element_blank()) +
                           theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                           theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                           scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                           scale_colour_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                           scale_fill_manual(values = c("#3C5F8D", "#2B4E7A", "#1B3D6B", "#0D2A51")) +
                           coord_cartesian(xlim = c(-1, 1)) +
                           annotate('text',  x = 1, y = (seq(1, dim(intercept_table_1)[1], 1)+0.4),
                           label= paste("italic(k)==", c(intercept_table_1["Aquatic", "K"], 
                                                         intercept_table_1["Individual-level Traits", "K"], 
                                                         intercept_table_1["Population-level Traits", "K"], 
                                                         intercept_table_1["Overall", "K"]), "~","(", 
                                                       c(intercept_table_1["Aquatic", "group_no"], 
                                                         intercept_table_1["Individual-level Traits", "group_no"], 
                                                         intercept_table_1["Population-level Traits", "group_no"], 
                                                         intercept_table_1["Overall", "group_no"]), 
                                       ")"), parse = TRUE, hjust = "right", size = 3.5) +
                           geom_label(aes(label=c(paste(format(round(mean(exp(Aquatic_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"),
                                                  paste(format(round(mean(exp(Individual_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"), 
                                                  paste(format(round(mean(exp(Population_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"),
                                                  paste(format(round(mean(exp(Overall_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%")), 
                                      x = -0.75, y = (seq(1, dim(intercept_table_1)[1], 1)+0.4)), size = 3.5)

density_intercept_CVR_1 #(400x400) 

# Preparing Data - Part 2

intercept_rnames_2 <- c("Terrestrial", "Acclimation", "Developmental")

intercept_means_df_2 <- data.frame("Mean" = c(Terrestrial_Model_CVR_Estimates$estimate, 
                                              Acclimation_Model_CVR_Estimates$estimate, 
                                              Developmental_Model_CVR_Estimates$estimate), 
                                   row.names = intercept_rnames_2)

intercept_low_df_2 <- data.frame("Low CI" = c(Terrestrial_Model_CVR_Estimates$ci.lb, 
                                              Acclimation_Model_CVR_Estimates$ci.lb, 
                                              Developmental_Model_CVR_Estimates$ci.lb), 
                                 row.names = intercept_rnames_2)

intercept_high_df_2 <- data.frame("High CI" = c(Terrestrial_Model_CVR_Estimates$ci.ub, 
                                                Acclimation_Model_CVR_Estimates$ci.ub, 
                                                Developmental_Model_CVR_Estimates$ci.ub), 
                                  row.names = intercept_rnames_2)

intercept_k_2 <- data.frame("k" = c(length(Terrestrial_Subset_Data$Effect_Size_ID), 
                                    length(Acclimation_Subset_Data$Effect_Size_ID),
                                    length(Developmental_Subset_Data$Effect_Size_ID)), 
                            row.names = intercept_rnames_2)

intercept_group_no_2 <- data.frame("Spp No." = c(length(unique(Terrestrial_Subset_Data$Scientific_Name)), 
                                                length(unique(Acclimation_Subset_Data$Scientific_Name)),
                                                length(unique(Developmental_Subset_Data$Scientific_Name))),
                                   row.names = intercept_rnames_2)

intercept_study_2 <- data.frame("Study" = c(length(unique(Terrestrial_Subset_Data$Study_ID)), 
                                            length(unique(Acclimation_Subset_Data$Study_ID)),
                                            length(unique(Developmental_Subset_Data$Study_ID))),
                                row.names = intercept_rnames_2)

intercept_table_2 <- data.frame(estimate = intercept_means_df_2[,1], 
                                lowerCL = intercept_low_df_2[,1], 
                                upperCL = intercept_high_df_2[,1], 
                                K = intercept_k_2[,1], 
                                group_no = intercept_group_no_2[,1], 
                                row.names = intercept_rnames_2)
intercept_table_2$name <- row.names(intercept_table_2)

intercept_raw_mean_2 <- c(Terrestrial_Subset_Data$InCVR, 
                          Acclimation_Subset_Data$InCVR, 
                          Developmental_Subset_Data$InCVR)

intercept_raw_name_2 <- c(replicate(929, "Terrestrial"), 
                          replicate(336, "Acclimation"), 
                          replicate(988, "Developmental"))

intercept_raw_df_2 <- data.frame("Model" = intercept_raw_name_2, 
                                 "Effect" = intercept_raw_mean_2)
# Graph code - Part 2
Intercept_Order_2 <- c("Developmental", "Acclimation", "Terrestrial")

density_intercept_CVR_2 <- intercept_table_2 %>% mutate(name = fct_relevel(name, Intercept_Order_2)) %>%
                           ggplot() +
                           geom_density_ridges(data = intercept_raw_df_2 %>% mutate(Model = fct_relevel(Model, Intercept_Order_2)), 
                                               aes(x = Effect, y = Model, colour = Model, fill = Model), 
                                                   scale = 0.8, alpha = 0.3, size = 1, inherit.aes = FALSE) +
                           geom_linerange(aes(y = rev(seq(1, dim(intercept_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, colour = name),
                                              size = 1) +
                           geom_pointrange(aes(x = estimate, y = rev(seq(1, dim(intercept_table_2)[1], 1)-0.1), xmin = lowerCL, xmax = upperCL, fill = name, colour = name), 
                                               size = 1, fatten = 2) +
                           theme_bw() +
                           guides(fill = "none", colour = "none") +
                           labs(x = TeX("Effect Size (lnCVR)"), y = "") +
                           theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, 
                                 vjust = c(-2.7, -2.7, -2.7))) +
                           theme(axis.text.x = element_text(margin = margin(b = 5))) +
                           theme(axis.ticks = element_blank()) +
                           theme(panel.grid.major.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                           theme(panel.grid.minor.x = element_line(colour = rgb(235, 235, 235, 150, maxColorValue = 500))) +
                           scale_y_discrete(expand = expansion(add = c(0.2, 1)), labels = function(x) str_wrap(x, width = 13)) +
                           scale_colour_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                           scale_fill_manual(values = c("#6582A9", "#5D7AA1", "#4A6E9C")) +
                           coord_cartesian(xlim = c(-1, 1)) +
                           annotate('text',  x = 1, y = (seq(1, dim(intercept_table_2)[1], 1)+0.4),
                           label= paste("italic(k)==", c(intercept_table_2["Developmental", "K"], 
                                                         intercept_table_2["Acclimation", "K"], 
                                                         intercept_table_2["Terrestrial", "K"]), "~","(", 
                                                       c(intercept_table_2["Developmental", "group_no"], 
                                                         intercept_table_2["Acclimation", "group_no"], 
                                                         intercept_table_2["Terrestrial", "group_no"]), 
                                         ")"), parse = TRUE, hjust = "right", size = 3.5) +
                           geom_label(aes(label=c(paste(format(round(mean(exp(Developmental_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"), 
                                                  paste(format(round(mean(exp(Acclimation_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%"),
                                                  paste(format(round(mean(exp(Terrestrial_Model_CVR_Estimates$estimate)-1)*100, 2), nsmall = 2), "%")), 
                                          x = -0.75, y = (seq(1, dim(intercept_table_2)[1], 1)+0.4)), size = 3.5)

density_intercept_CVR_2 #(400x320)

##### Supplementary Material Tables #####

# Consistency Changes - Studies, Species and Effect Sizes Counts

Pre_Data <- read.csv("./3.Data_Analysis/2.Outputs/Data/Pre_Data.csv")

Individual_Pre_Data <- Pre_Data %>% filter(Trait_Category != "Population")
Developmental_Pre_Data <- Individual_Pre_Data %>% filter(Plasticity_Mechanism == "Developmental Plasticity")

Developmental_Exposure_Time_Studies <- Developmental_Pre_Data %>% select("Study_ID", "Developmental_Exposure_Time") %>% table() %>% data.frame() %>% 
                                       filter(`Freq` != 0) %>% select("Developmental_Exposure_Time") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Studies) <- Developmental_Exposure_Time_Studies$Developmental_Exposure_Time
colnames(Developmental_Exposure_Time_Studies) <- c("Developmental_Exposure_Time", "Study")

Developmental_Exposure_Time_Species <- Developmental_Pre_Data %>% select("Scientific_Name", "Developmental_Exposure_Time") %>% table() %>% data.frame() %>% 
                                       filter(`Freq` != 0) %>% select("Developmental_Exposure_Time") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Species) <- Developmental_Exposure_Time_Species$Developmental_Exposure_Time
colnames(Developmental_Exposure_Time_Species) <- c("Developmental_Exposure_Time", "Species")

Developmental_Exposure_Time_Effects <- Developmental_Pre_Data %>% select("Developmental_Exposure_Time") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Effects) <- Developmental_Exposure_Time_Effects$Developmental_Exposure_Time
colnames(Developmental_Exposure_Time_Effects) <- c("Developmental_Exposure_Time", "Effect Sizes")

Developmental_Exposure_Time_Final_Counts <- Developmental_Exposure_Time_Studies %>% 
                                            left_join(Developmental_Exposure_Time_Species, by = "Developmental_Exposure_Time") %>% 
                                            left_join(Developmental_Exposure_Time_Effects, by = "Developmental_Exposure_Time")

write.csv(Developmental_Exposure_Time_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Developmental_Exposure_Time_Final_Counts.csv", row.names = FALSE)

Acclimation_Pre_Data <- Individual_Pre_Data %>% filter(Plasticity_Mechanism == "Acclimation")

Acclimation_LH_Stage_Studies <- Acclimation_Pre_Data %>% select("Study_ID", "Acclimation_Life.History_Stage") %>% table() %>% data.frame() %>% 
                                filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Studies) <- Acclimation_LH_Stage_Studies$Acclimation_Life.History_Stage
colnames(Acclimation_LH_Stage_Studies) <- c("Acclimation_Life.History_Stage", "Study")

Acclimation_LH_Stage_Species <- Acclimation_Pre_Data %>% select("Scientific_Name", "Acclimation_Life.History_Stage") %>% table() %>% data.frame() %>% 
                                filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Species) <- Acclimation_LH_Stage_Species$Acclimation_Life.History_Stage
colnames(Acclimation_LH_Stage_Species) <- c("Acclimation_Life.History_Stage", "Species")

Acclimation_LH_Stage_Effects <- Acclimation_Pre_Data %>% select("Acclimation_Life.History_Stage") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Effects) <- Acclimation_LH_Stage_Effects$Acclimation_Life.History_Stage
colnames(Acclimation_LH_Stage_Effects) <- c("Acclimation_Life.History_Stage", "Effect Sizes")

Acclimation_LH_Stage_Final_Counts <- Acclimation_LH_Stage_Studies %>% 
                                     left_join(Acclimation_LH_Stage_Species, by = "Acclimation_Life.History_Stage") %>% 
                                     left_join(Acclimation_LH_Stage_Effects, by = "Acclimation_Life.History_Stage")

write.csv(Acclimation_LH_Stage_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Acclimation_LH_Stage_Final_Counts.csv", row.names = FALSE)

Measurement_Studies <- Pre_Data %>% select("Study_ID", "Measurement") %>% table() %>% data.frame() %>% 
                       filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Measurement_Studies) <- Measurement_Studies$Measurement
colnames(Measurement_Studies) <- c("Measurement", "Study")

Measurement_Species <- Pre_Data %>% select("Scientific_Name", "Measurement") %>% table() %>% data.frame() %>% 
                       filter(`Freq` != 0) %>% select("Measurement") %>% table() %>% data.frame()
rownames(Measurement_Species) <- Measurement_Species$Measurement
colnames(Measurement_Species) <- c("Measurement", "Species")

Measurement_Effects <- Pre_Data %>% select("Measurement") %>% table() %>% data.frame()
rownames(Measurement_Effects) <- Measurement_Effects$Measurement
colnames(Measurement_Effects) <- c("Measurement", "Effect Sizes")

Measurement_Final_Counts <- Measurement_Studies %>% 
                            left_join(Measurement_Species, by = "Measurement") %>% 
                            left_join(Measurement_Effects, by = "Measurement")

write.csv(Measurement_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Measurement_Final_Counts.csv", row.names = FALSE)

# Category - Studies, Species and Effect Sizes

Developmental_Exposure_Time_Category_Studies <- Developmental_Subset_Data %>% select("Study_ID", "Developmental_Exposure_Time_Category") %>% table() %>% data.frame() %>% 
                                                filter(`Freq` != 0) %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Category_Studies) <- Developmental_Exposure_Time_Category_Studies$Developmental_Exposure_Time_Category
colnames(Developmental_Exposure_Time_Category_Studies) <- c("Developmental_Exposure_Time_Category", "Study")

Developmental_Exposure_Time_Category_Species <- Developmental_Subset_Data %>% select("Scientific_Name", "Developmental_Exposure_Time_Category") %>% table() %>% data.frame() %>% 
                                                filter(`Freq` != 0) %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Category_Species) <- Developmental_Exposure_Time_Category_Species$Developmental_Exposure_Time_Category
colnames(Developmental_Exposure_Time_Category_Species) <- c("Developmental_Exposure_Time_Category", "Species")

Developmental_Exposure_Time_Category_Effects <- Developmental_Subset_Data %>% select("Developmental_Exposure_Time_Category") %>% table() %>% data.frame()
rownames(Developmental_Exposure_Time_Category_Effects) <- Developmental_Exposure_Time_Category_Effects$Developmental_Exposure_Time_Category
colnames(Developmental_Exposure_Time_Category_Effects) <- c("Developmental_Exposure_Time_Category", "Effect Sizes")

Developmental_Exposure_Time_Category_Final_Counts <- Developmental_Exposure_Time_Category_Studies %>% 
                                                     left_join(Developmental_Exposure_Time_Category_Species, by = "Developmental_Exposure_Time_Category") %>% 
                                                     left_join(Developmental_Exposure_Time_Category_Effects, by = "Developmental_Exposure_Time_Category")

write.csv(Developmental_Exposure_Time_Category_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Developmental_Exposure_Time_Category_Final_Counts.csv", row.names = FALSE)

Acclimation_LH_Stage_Category_Studies <- Acclimation_Subset_Data %>% select("Study_ID", "Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame() %>% 
                                         filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Category_Studies) <- Acclimation_LH_Stage_Category_Studies$Acclimation_Life.History_Stage_Category
colnames(Acclimation_LH_Stage_Category_Studies) <- c("Acclimation_Life.History_Stage_Category", "Study")

Acclimation_LH_Stage_Category_Species <- Acclimation_Subset_Data %>% select("Scientific_Name", "Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame() %>% 
                                         filter(`Freq` != 0) %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Category_Species) <- Acclimation_LH_Stage_Category_Species$Acclimation_Life.History_Stage_Category
colnames(Acclimation_LH_Stage_Category_Species) <- c("Acclimation_Life.History_Stage_Category", "Species")

Acclimation_LH_Stage_Category_Effects <- Acclimation_Subset_Data %>% select("Acclimation_Life.History_Stage_Category") %>% table() %>% data.frame()
rownames(Acclimation_LH_Stage_Category_Effects) <- Acclimation_LH_Stage_Category_Effects$Acclimation_Life.History_Stage_Category
colnames(Acclimation_LH_Stage_Category_Effects) <- c("Acclimation_Life.History_Stage_Category", "Effect Sizes")

Acclimation_LH_Stage_Category_Final_Counts <- Acclimation_LH_Stage_Category_Studies %>% 
                                              left_join(Acclimation_LH_Stage_Category_Species, by = "Acclimation_Life.History_Stage_Category") %>% 
                                              left_join(Acclimation_LH_Stage_Category_Effects, by = "Acclimation_Life.History_Stage_Category")

write.csv(Acclimation_LH_Stage_Category_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Acclimation_LH_Stage_Category_Final_Counts.csv", row.names = FALSE)

Measurement_Category_Studies <- data %>% select("Study_ID", "Trait_Category") %>% table() %>% data.frame() %>% 
                                filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Measurement_Category_Studies) <- Measurement_Category_Studies$Trait_Category
colnames(Measurement_Category_Studies) <- c("Trait_Category", "Study")

Measurement_Category_Species <- data %>% select("Scientific_Name", "Trait_Category") %>% table() %>% data.frame() %>% 
                                filter(`Freq` != 0) %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Measurement_Category_Species) <- Measurement_Category_Species$Trait_Category
colnames(Measurement_Category_Species) <- c("Trait_Category", "Species")

Measurement_Category_Effects <- data %>% select("Trait_Category") %>% table() %>% data.frame()
rownames(Measurement_Category_Effects) <- Measurement_Category_Effects$Trait_Category
colnames(Measurement_Category_Effects) <- c("Trait_Category", "Effect Sizes")

Measurement_Category_Final_Counts <- Measurement_Category_Studies %>% 
                                     left_join(Measurement_Category_Species, by = "Trait_Category") %>% 
                                     left_join(Measurement_Category_Effects, by = "Trait_Category")

write.csv(Measurement_Category_Final_Counts, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Measurement_Category_Final_Counts.csv", row.names = FALSE)

# Phylogenetic Tree with labels

labelled_tree <- tree

Scientific_Name_Effects <- data %>% select("Scientific_Name") %>% table() %>% data.frame()
rownames(Scientific_Name_Effects) <- Scientific_Name_Effects$Scientific_Name
colnames(Scientific_Name_Effects) <- c("Scientific_Name", "Effect_Sizes")
Scientific_Name_Effects <- Scientific_Name_Effects[c(labelled_tree$tip.label), ]

Scientific_Name_Studies <- data %>% select("Study_ID", "Scientific_Name") %>% table() %>% data.frame() %>% 
                           filter(`Freq` != 0) %>% select("Scientific_Name") %>% table() %>% data.frame()
rownames(Scientific_Name_Studies) <- Scientific_Name_Studies$Scientific_Name
colnames(Scientific_Name_Studies) <- c("Scientific_Name", "Study")
Scientific_Name_Studies <- Scientific_Name_Studies[c(labelled_tree$tip.label), ]

labelled_tree$tip.label <- paste(labelled_tree$tip.label, " ", Scientific_Name_Effects$Effect_Sizes, "(", Scientific_Name_Studies$Study, ")")
node.depth(labelled_tree, method = 2)
plot(labelled_tree, node.color = "#183357")

# Raw Data Table

Raw_Overall_CVR <- data.frame("Overall" = c("MLMA"),
                              "Studies" = c(intercept_study["Overall", "Study"]), 
                              "Species" = c(intercept_table["Overall", "group_no"]), 
                              "Effect Sizes" = c(intercept_table["Overall", "K"]),
                              "Estimate" = c(Overall_Model_CVR$b[1]),
                              "CI Low" = c(Overall_Model_CVR$ci.lb), 
                              "CI High" = c(Overall_Model_CVR$ci.ub), 
                              "df" = c(Overall_Model_CVR$ddf[[1]]), 
                              "p-value" = c(Overall_Model_CVR$pval))

Raw_Amplitude_CVR <- data.frame("Overall" = c("Fluctuation Amplitude"),
                                "Studies" = c(length(unique(data$Study_ID))), 
                                "Species" = c(length(unique(data$Scientific_Name))), 
                                "Effect Sizes" = c(length(data$Effect_Size_ID)),
                                "Estimate" = c(Amplitude_Model_CVR$b[1]),
                                "CI Low" = c(Amplitude_Model_CVR$ci.lb), 
                                "CI High" = c(Amplitude_Model_CVR$ci.ub), 
                                "df" = c(Amplitude_Model_CVR$ddf[[1]]), 
                                "p-value" = c(Amplitude_Model_CVR$pval))

Raw_Fluctuation_Type_CVR <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating", 
                                                              "Stepwise", "Stochastic"),
                                       "Studies" = c(fluctuation_study["Sinusoidal (Sine Curve)", "Study"], fluctuation_study["Alternating", "Study"], 
                                                     fluctuation_study["Stepwise", "Study"], fluctuation_study["Stochastic", "Study"]), 
                                       "Species" = c(fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], fluctuation_table["Alternating", "group_no"], 
                                                     fluctuation_table["Stepwise", "group_no"], fluctuation_table["Stochastic", "group_no"]), 
                                       "Effect Sizes" = c(fluctuation_table["Sinusoidal (Sine Curve)", "K"], fluctuation_table["Alternating", "K"], 
                                                          fluctuation_table["Stepwise", "K"], fluctuation_table["Stochastic", "K"]),
                                       "Estimate" = c(Fluctuation_Model_CVR$b[[2]], Fluctuation_Model_CVR$b[[1]],
                                                      Fluctuation_Model_CVR$b[[3]], Fluctuation_Model_CVR$b[[4]]),
                                       "CI Low" = c(Fluctuation_Model_CVR$ci.lb[2], Fluctuation_Model_CVR$ci.lb[1], 
                                                    Fluctuation_Model_CVR$ci.lb[3], Fluctuation_Model_CVR$ci.lb[4]), 
                                       "CI High" = c(Fluctuation_Model_CVR$ci.ub[2], Fluctuation_Model_CVR$ci.ub[1], 
                                                     Fluctuation_Model_CVR$ci.ub[3], Fluctuation_Model_CVR$ci.ub[4]), 
                                       "df" = c(Fluctuation_Model_CVR$ddf[[2]], Fluctuation_Model_CVR$ddf[[1]], 
                                                Fluctuation_Model_CVR$ddf[[3]], Fluctuation_Model_CVR$ddf[[4]]), 
                                       "p-value" = c(Fluctuation_Model_CVR$pval[2], Fluctuation_Model_CVR$pval[1], 
                                                     Fluctuation_Model_CVR$pval[3], Fluctuation_Model_CVR$pval[4]))

Raw_Trait_CVR <- data.frame("Phenotypic Trait Categories" = c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits",  
                                                              "Morphology", "Physiological", "Population"),
                            "Studies" = c(trait_study["Behavioural", "Study"], trait_study["Biochemical Assay", "Study"], trait_study["Gene Expression", "Study"], trait_study["Life-history Traits", "Study"],
                                          trait_study["Morphology", "Study"], trait_study["Physiological", "Study"], trait_study["Population", "Study"]), 
                            "Species" = c(trait_table["Behavioural", "group_no"], trait_table["Biochemical Assay", "group_no"], trait_table["Gene Expression", "group_no"], trait_table["Life-history Traits", "group_no"],
                                          trait_table["Morphology", "group_no"], trait_table["Physiological", "group_no"], trait_table["Population", "group_no"]), 
                            "Effect Sizes" = c(trait_table["Behavioural", "K"], trait_table["Biochemical Assay", "K"], trait_table["Gene Expression", "K"], trait_table["Life-history Traits", "K"],
                                               trait_table["Morphology", "K"], trait_table["Physiological", "K"], trait_table["Population", "K"]),
                            "Estimate" = c(Trait_Model_CVR$b[[1]], Trait_Model_CVR$b[[2]], Trait_Model_CVR$b[[3]], Trait_Model_CVR$b[[4]],
                                           Trait_Model_CVR$b[[5]], Trait_Model_CVR$b[[6]], Trait_Model_CVR$b[[7]]),
                            "CI Low" = c(Trait_Model_CVR$ci.lb[1], Trait_Model_CVR$ci.lb[2], Trait_Model_CVR$ci.lb[3], Trait_Model_CVR$ci.lb[4], 
                                         Trait_Model_CVR$ci.lb[5], Trait_Model_CVR$ci.lb[6], Trait_Model_CVR$ci.lb[7]), 
                            "CI High" = c(Trait_Model_CVR$ci.ub[1], Trait_Model_CVR$ci.ub[2], Trait_Model_CVR$ci.ub[3], Trait_Model_CVR$ci.ub[4], 
                                          Trait_Model_CVR$ci.ub[5], Trait_Model_CVR$ci.ub[6], Trait_Model_CVR$ci.ub[7]), 
                            "df" = c(Trait_Model_CVR$ddf[[1]], Trait_Model_CVR$ddf[[2]], Trait_Model_CVR$ddf[[3]], Trait_Model_CVR$ddf[[4]], 
                                     Trait_Model_CVR$ddf[[5]], Trait_Model_CVR$ddf[[6]], Trait_Model_CVR$ddf[[7]]), 
                            "p-value" = c(Trait_Model_CVR$pval[1], Trait_Model_CVR$pval[2], Trait_Model_CVR$pval[3], Trait_Model_CVR$pval[4], 
                                          Trait_Model_CVR$pval[5], Trait_Model_CVR$pval[6], Trait_Model_CVR$pval[7]))

Raw_Class_CVR <- data.frame("Taxonomic Class" = c("Actinopteri", "Amphibia", "Anthozoa", "Arachnida", "Bivalvia", 
                                                  "Branchiopoda", "Gastropoda", "Holothuroidea", "Insecta", "Malacostraca"),
                            "Studies" = c(class_study["Actinopteri", "Study"], class_study["Amphibia", "Study"], class_study["Anthozoa", "Study"], 
                                          class_study["Arachnida", "Study"], class_study["Bivalvia", "Study"], class_study["Branchiopoda", "Study"], 
                                          class_study["Gastropoda", "Study"], class_study["Holothuroidea", "Study"], class_study["Insecta", "Study"], 
                                          class_study["Malacostraca", "Study"]), 
                            "Species" = c(class_table["Actinopteri", "group_no"], class_table["Amphibia", "group_no"], class_table["Anthozoa", "group_no"], 
                                          class_table["Arachnida", "group_no"], class_table["Bivalvia", "group_no"], class_table["Branchiopoda", "group_no"], 
                                          class_table["Gastropoda", "group_no"], class_table["Holothuroidea", "group_no"], class_table["Insecta", "group_no"], 
                                          class_table["Malacostraca", "group_no"]), 
                            "Effect Sizes" = c(class_table["Actinopteri", "K"], class_table["Amphibia", "K"], class_table["Anthozoa", "K"], 
                                               class_table["Arachnida", "K"], class_table["Bivalvia", "K"], class_table["Branchiopoda", "K"], 
                                               class_table["Gastropoda", "K"], class_table["Holothuroidea", "K"], class_table["Insecta", "K"], 
                                               class_table["Malacostraca", "K"]),
                            "Estimate" = c(Class_Model_CVR$b[[1]], Class_Model_CVR$b[[2]], Class_Model_CVR$b[[3]], 
                                           Class_Model_CVR$b[[4]], Class_Model_CVR$b[[5]], Class_Model_CVR$b[[6]], 
                                           Class_Model_CVR$b[[7]], Class_Model_CVR$b[[8]], Class_Model_CVR$b[[9]], 
                                           Class_Model_CVR$b[[10]]),
                            "CI Low" = c(Class_Model_CVR$ci.lb[1], Class_Model_CVR$ci.lb[2], Class_Model_CVR$ci.lb[3], 
                                         Class_Model_CVR$ci.lb[4], Class_Model_CVR$ci.lb[5], Class_Model_CVR$ci.lb[6], 
                                         Class_Model_CVR$ci.lb[7], Class_Model_CVR$ci.lb[8], Class_Model_CVR$ci.lb[9], 
                                         Class_Model_CVR$ci.lb[10]), 
                            "CI High" = c(Class_Model_CVR$ci.ub[1], Class_Model_CVR$ci.ub[2], Class_Model_CVR$ci.ub[3], 
                                          Class_Model_CVR$ci.ub[4], Class_Model_CVR$ci.ub[5], Class_Model_CVR$ci.ub[6], 
                                          Class_Model_CVR$ci.ub[7], Class_Model_CVR$ci.ub[8], Class_Model_CVR$ci.ub[9], 
                                          Class_Model_CVR$ci.ub[10]), 
                            "df" = c(Class_Model_CVR$ddf[[1]], Class_Model_CVR$ddf[[2]], Class_Model_CVR$ddf[[3]], 
                                     Class_Model_CVR$ddf[[4]], Class_Model_CVR$ddf[[5]], Class_Model_CVR$ddf[[6]], 
                                     Class_Model_CVR$ddf[[7]], Class_Model_CVR$ddf[[8]], Class_Model_CVR$ddf[[9]], 
                                     Class_Model_CVR$ddf[[10]]), 
                            "p-value" = c(Class_Model_CVR$pval[1], Class_Model_CVR$pval[2], Class_Model_CVR$pval[3], 
                                          Class_Model_CVR$pval[4], Class_Model_CVR$pval[5], Class_Model_CVR$pval[6], 
                                          Class_Model_CVR$pval[7], Class_Model_CVR$pval[8], Class_Model_CVR$pval[9], 
                                          Class_Model_CVR$pval[10]))

Raw_Specific_Trait_CVR <- data.frame("Specific Phenotypic Traits" = c("Apparent Digestibility Coefficient", "Catalase Activity", "Cortisol Levels", "Development Time", "Fecundity", 
                                                                      "Food Consumption", "Head Width", "hsp70", "Immune Defense", "Length", 
                                                                      "Locomotor Performance", "Longevity", "Mass", "Metabolic Rate", "Mortality", 
                                                                      "PO Activity", "Reproductive Rate", "Sex", "SOD Activity", "Survival", 
                                                                      "Tail Length", "Triglyceride"),
                                     "Studies" = c(specific_trait_study["Apparent Digestability Coefficient", "Study"], specific_trait_study["Catalase Activity", "Study"], specific_trait_study["Cortisol", "Study"], 
                                                   specific_trait_study["Development Time", "Study"], specific_trait_study["Fecundity", "Study"], specific_trait_study["Food Consumption", "Study"], 
                                                   specific_trait_study["Head Width", "Study"], specific_trait_study["hsp70", "Study"], specific_trait_study["Immune Defense", "Study"], 
                                                   specific_trait_study["Length", "Study"], specific_trait_study["Locomotor Performance", "Study"], 
                                                   specific_trait_study_2["Longevity", "Study"], specific_trait_study_2["Mass", "Study"], specific_trait_study_2["Metabolic Rate", "Study"], 
                                                   specific_trait_study_2["Mortality", "Study"], specific_trait_study_2["PO Activity", "Study"], specific_trait_study_2["Reproductive Rate", "Study"], 
                                                   specific_trait_study_2["Sex", "Study"], specific_trait_study_2["SOD Activity", "Study"], specific_trait_study_2["Survival", "Study"], 
                                                   specific_trait_study_2["Tail Length", "Study"], specific_trait_study_2["Triglyceride", "Study"]), 
                                     "Species" = c(specific_trait_table["Apparent Digestability Coefficient", "group_no"], specific_trait_table["Catalase Activity", "group_no"], specific_trait_table["Cortisol", "group_no"], 
                                                   specific_trait_table["Development Time", "group_no"], specific_trait_table["Fecundity", "group_no"], specific_trait_table["Food Consumption", "group_no"], 
                                                   specific_trait_table["Head Width", "group_no"], specific_trait_table["hsp70", "group_no"], specific_trait_table["Immune Defense", "group_no"], 
                                                   specific_trait_table["Length", "group_no"], specific_trait_table["Locomotor Performance", "group_no"], 
                                                   specific_trait_table_2["Longevity", "group_no"], specific_trait_table_2["Mass", "group_no"], specific_trait_table_2["Metabolic Rate", "group_no"], 
                                                   specific_trait_table_2["Mortality", "group_no"], specific_trait_table_2["PO Activity", "group_no"], specific_trait_table_2["Reproductive Rate", "group_no"], 
                                                   specific_trait_table_2["Sex", "group_no"], specific_trait_table_2["SOD Activity", "group_no"], specific_trait_table_2["Survival", "group_no"], 
                                                   specific_trait_table_2["Tail Length", "group_no"], specific_trait_table_2["Triglyceride", "group_no"]), 
                                     "Effect Sizes" = c(specific_trait_table["Apparent Digestability Coefficient", "K"], specific_trait_table["Catalase Activity", "K"], specific_trait_table["Cortisol", "K"], 
                                                        specific_trait_table["Development Time", "K"], specific_trait_table["Fecundity", "K"], specific_trait_table["Food Consumption", "K"], 
                                                        specific_trait_table["Head Width", "K"], specific_trait_table["hsp70", "K"], specific_trait_table["Immune Defense", "K"], 
                                                        specific_trait_table["Length", "K"], specific_trait_table["Locomotor Performance", "K"], 
                                                        specific_trait_table_2["Longevity", "K"], specific_trait_table_2["Mass", "K"], specific_trait_table_2["Metabolic Rate", "K"], 
                                                        specific_trait_table_2["Mortality", "K"], specific_trait_table_2["PO Activity", "K"], specific_trait_table_2["Reproductive Rate", "K"], 
                                                        specific_trait_table_2["Sex", "K"], specific_trait_table_2["SOD Activity", "K"], specific_trait_table_2["Survival", "K"], 
                                                        specific_trait_table_2["Tail Length", "K"], specific_trait_table_2["Triglyceride", "K"]),
                                     "Estimate" = c(Specific_Trait_Model_CVR$b[[1]], Specific_Trait_Model_CVR$b[[2]], Specific_Trait_Model_CVR$b[[3]], 
                                                    Specific_Trait_Model_CVR$b[[4]], Specific_Trait_Model_CVR$b[[5]], Specific_Trait_Model_CVR$b[[6]], 
                                                    Specific_Trait_Model_CVR$b[[8]], Specific_Trait_Model_CVR$b[[9]], Specific_Trait_Model_CVR$b[[10]], 
                                                    Specific_Trait_Model_CVR$b[[11]], Specific_Trait_Model_CVR$b[[12]], Specific_Trait_Model_CVR$b[[13]], 
                                                    Specific_Trait_Model_CVR$b[[14]], Specific_Trait_Model_CVR$b[[15]], Specific_Trait_Model_CVR$b[[16]], 
                                                    Specific_Trait_Model_CVR$b[[17]], Specific_Trait_Model_CVR$b[[18]], Specific_Trait_Model_CVR$b[[7]], 
                                                    Specific_Trait_Model_CVR$b[[19]], Specific_Trait_Model_CVR$b[[20]], Specific_Trait_Model_CVR$b[[21]], 
                                                    Specific_Trait_Model_CVR$b[[22]]),
                                     "CI Low" = c(Specific_Trait_Model_CVR$ci.lb[1], Specific_Trait_Model_CVR$ci.lb[2], Specific_Trait_Model_CVR$ci.lb[3], 
                                                  Specific_Trait_Model_CVR$ci.lb[4], Specific_Trait_Model_CVR$ci.lb[5], Specific_Trait_Model_CVR$ci.lb[6], 
                                                  Specific_Trait_Model_CVR$ci.lb[8], Specific_Trait_Model_CVR$ci.lb[9], Specific_Trait_Model_CVR$ci.lb[10], 
                                                  Specific_Trait_Model_CVR$ci.lb[11], Specific_Trait_Model_CVR$ci.lb[12], Specific_Trait_Model_CVR$ci.lb[13], 
                                                  Specific_Trait_Model_CVR$ci.lb[14], Specific_Trait_Model_CVR$ci.lb[15], Specific_Trait_Model_CVR$ci.lb[16], 
                                                  Specific_Trait_Model_CVR$ci.lb[17], Specific_Trait_Model_CVR$ci.lb[18], Specific_Trait_Model_CVR$ci.lb[7], 
                                                  Specific_Trait_Model_CVR$ci.lb[19], Specific_Trait_Model_CVR$ci.lb[20], Specific_Trait_Model_CVR$ci.lb[21], 
                                                  Specific_Trait_Model_CVR$ci.lb[22]), 
                                     "CI High" = c(Specific_Trait_Model_CVR$ci.ub[1], Specific_Trait_Model_CVR$ci.ub[2], Specific_Trait_Model_CVR$ci.ub[3], 
                                                   Specific_Trait_Model_CVR$ci.ub[4], Specific_Trait_Model_CVR$ci.ub[5], Specific_Trait_Model_CVR$ci.ub[6], 
                                                   Specific_Trait_Model_CVR$ci.ub[8], Specific_Trait_Model_CVR$ci.ub[9], Specific_Trait_Model_CVR$ci.ub[10], 
                                                   Specific_Trait_Model_CVR$ci.ub[11], Specific_Trait_Model_CVR$ci.ub[12], Specific_Trait_Model_CVR$ci.ub[13], 
                                                   Specific_Trait_Model_CVR$ci.ub[14], Specific_Trait_Model_CVR$ci.ub[15], Specific_Trait_Model_CVR$ci.ub[16], 
                                                   Specific_Trait_Model_CVR$ci.ub[17], Specific_Trait_Model_CVR$ci.ub[18], Specific_Trait_Model_CVR$ci.ub[7], 
                                                   Specific_Trait_Model_CVR$ci.ub[19], Specific_Trait_Model_CVR$ci.ub[20], Specific_Trait_Model_CVR$ci.ub[21], 
                                                   Specific_Trait_Model_CVR$ci.ub[22]), 
                                     "df" = c(Specific_Trait_Model_CVR$ddf[[1]], Specific_Trait_Model_CVR$ddf[[2]], Specific_Trait_Model_CVR$ddf[[3]], 
                                              Specific_Trait_Model_CVR$ddf[[4]], Specific_Trait_Model_CVR$ddf[[5]], Specific_Trait_Model_CVR$ddf[[6]], 
                                              Specific_Trait_Model_CVR$ddf[[8]], Specific_Trait_Model_CVR$ddf[[9]], Specific_Trait_Model_CVR$ddf[[10]], 
                                              Specific_Trait_Model_CVR$ddf[[11]], Specific_Trait_Model_CVR$ddf[[12]], Specific_Trait_Model_CVR$ddf[[13]], 
                                              Specific_Trait_Model_CVR$ddf[[14]], Specific_Trait_Model_CVR$ddf[[15]], Specific_Trait_Model_CVR$ddf[[16]], 
                                              Specific_Trait_Model_CVR$ddf[[17]], Specific_Trait_Model_CVR$ddf[[18]], Specific_Trait_Model_CVR$ddf[[7]], 
                                              Specific_Trait_Model_CVR$ddf[[19]], Specific_Trait_Model_CVR$ddf[[20]], Specific_Trait_Model_CVR$ddf[[21]], 
                                              Specific_Trait_Model_CVR$ddf[[22]]), 
                                     "p-value" = c(Specific_Trait_Model_CVR$pval[1], Specific_Trait_Model_CVR$pval[2], Specific_Trait_Model_CVR$pval[3], 
                                                   Specific_Trait_Model_CVR$pval[4], Specific_Trait_Model_CVR$pval[5], Specific_Trait_Model_CVR$pval[6], 
                                                   Specific_Trait_Model_CVR$pval[8], Specific_Trait_Model_CVR$pval[9], Specific_Trait_Model_CVR$pval[10], 
                                                   Specific_Trait_Model_CVR$pval[11], Specific_Trait_Model_CVR$pval[12], Specific_Trait_Model_CVR$pval[13], 
                                                   Specific_Trait_Model_CVR$pval[14], Specific_Trait_Model_CVR$pval[15], Specific_Trait_Model_CVR$pval[16], 
                                                   Specific_Trait_Model_CVR$pval[17], Specific_Trait_Model_CVR$pval[18], Specific_Trait_Model_CVR$pval[7], 
                                                   Specific_Trait_Model_CVR$pval[19], Specific_Trait_Model_CVR$pval[20], Specific_Trait_Model_CVR$pval[21], 
                                                   Specific_Trait_Model_CVR$pval[22]))

Raw_Population_CVR <- data.frame("Population-level Traits" = c("MLMA"),
                                 "Studies" = c(intercept_study["Population-level Traits", "Study"]), 
                                 "Species" = c(intercept_table["Population-level Traits", "group_no"]), 
                                 "Effect Sizes" = c(intercept_table["Population-level Traits", "K"]),
                                 "Estimate" = c(Population_Model_CVR$b[1]),
                                 "CI Low" = c(Population_Model_CVR$ci.lb), 
                                 "CI High" = c(Population_Model_CVR$ci.ub), 
                                 "df" = c(Population_Model_CVR$ddf[[1]]), 
                                 "p-value" = c(Population_Model_CVR$pval))

Raw_Population_Amplitude_CVR <- data.frame("Population-level Traits" = c("Fluctuation Amplitude"),
                                           "Studies" = c(length(unique(Population_Subset_Data$Study_ID))), 
                                           "Species" = c(length(unique(Population_Subset_Data$Scientific_Name))), 
                                           "Effect Sizes" = c(length(Population_Subset_Data$Effect_Size_ID)),
                                           "Estimate" = c(Population_Amplitude_Model_CVR$b[1]),
                                           "CI Low" = c(Population_Amplitude_Model_CVR$ci.lb), 
                                           "CI High" = c(Population_Amplitude_Model_CVR$ci.ub), 
                                           "df" = c(Population_Amplitude_Model_CVR$ddf[[1]]), 
                                           "p-value" = c(Population_Amplitude_Model_CVR$pval))

Raw_Population_Fluctuation_Type_CVR <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating", 
                                                                         "Stepwise"),
                                                  "Studies" = c(population_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], population_fluctuation_study["Alternating", "Study"], 
                                                                population_fluctuation_study["Stepwise", "Study"]), 
                                                  "Species" = c(population_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], population_fluctuation_table["Alternating", "group_no"], 
                                                                population_fluctuation_table["Stepwise", "group_no"]), 
                                                  "Effect Sizes" = c(population_fluctuation_table["Sinusoidal (Sine Curve)", "K"], population_fluctuation_table["Alternating", "K"], 
                                                                     population_fluctuation_table["Stepwise", "K"]),
                                                  "Estimate" = c(Population_Fluctuation_Model_CVR$b[[2]], Population_Fluctuation_Model_CVR$b[[1]],
                                                                 Population_Fluctuation_Model_CVR$b[[3]]),
                                                  "CI Low" = c(Population_Fluctuation_Model_CVR$ci.lb[2], Population_Fluctuation_Model_CVR$ci.lb[1], 
                                                               Population_Fluctuation_Model_CVR$ci.lb[3]), 
                                                  "CI High" = c(Population_Fluctuation_Model_CVR$ci.ub[2], Population_Fluctuation_Model_CVR$ci.ub[1], 
                                                                Population_Fluctuation_Model_CVR$ci.ub[3]), 
                                                  "df" = c(Population_Fluctuation_Model_CVR$ddf[[2]], Population_Fluctuation_Model_CVR$ddf[[1]], 
                                                           Population_Fluctuation_Model_CVR$ddf[[3]]), 
                                                  "p-value" = c(Population_Fluctuation_Model_CVR$pval[2], Population_Fluctuation_Model_CVR$pval[1], 
                                                                Population_Fluctuation_Model_CVR$pval[3]))

Raw_Population_Class_CVR <- data.frame("Taxonomic Class" = c("Actinopteri", "Arachnida", "Insecta"),
                                       "Studies" = c(population_class_study["Actinopteri", "Study"], population_class_study["Arachnida", "Study"], population_class_study["Insecta", "Study"]), 
                                       "Species" = c(population_class_table["Actinopteri", "group_no"], population_class_table["Arachnida", "group_no"], population_class_table["Insecta", "group_no"]), 
                                       "Effect Sizes" = c(population_class_table["Actinopteri", "K"], population_class_table["Arachnida", "K"], population_class_table["Insecta", "K"]),
                                       "Estimate" = c(Population_Class_Model_CVR$b[[1]], Population_Class_Model_CVR$b[[2]], Population_Class_Model_CVR$b[[3]]),
                                       "CI Low" = c(Population_Class_Model_CVR$ci.lb[1], Population_Class_Model_CVR$ci.lb[2], Population_Class_Model_CVR$ci.lb[3]), 
                                       "CI High" = c(Population_Class_Model_CVR$ci.ub[1], Population_Class_Model_CVR$ci.ub[2], Population_Class_Model_CVR$ci.ub[3]), 
                                       "df" = c(Population_Class_Model_CVR$ddf[[1]], Population_Class_Model_CVR$ddf[[2]], Population_Class_Model_CVR$ddf[[3]]), 
                                       "p-value" = c(Population_Class_Model_CVR$pval[1], Population_Class_Model_CVR$pval[2], Population_Class_Model_CVR$pval[3]))

Raw_Individual_CVR <- data.frame("Individual-level Traits" = c("MLMA"),
                                 "Studies" = c(intercept_study["Individual-level Traits", "Study"]), 
                                 "Species" = c(intercept_table["Individual-level Traits", "group_no"]), 
                                 "Effect Sizes" = c(intercept_table["Individual-level Traits", "K"]),
                                 "Estimate" = c(Individual_Model_CVR$b[1]),
                                 "CI Low" = c(Individual_Model_CVR$ci.lb), 
                                 "CI High" = c(Individual_Model_CVR$ci.ub), 
                                 "df" = c(Individual_Model_CVR$ddf[[1]]), 
                                 "p-value" = c(Individual_Model_CVR$pval))

Raw_Individual_Amplitude_CVR <- data.frame("Individual-level Traits" = c("Fluctuation Amplitude"),
                                           "Studies" = c(length(unique(Individual_Subset_Data$Study_ID))), 
                                           "Species" = c(length(unique(Individual_Subset_Data$Scientific_Name))), 
                                           "Effect Sizes" = c(length(Individual_Subset_Data$Effect_Size_ID)),
                                           "Estimate" = c(Individual_Amplitude_Model_CVR$b[1]),
                                           "CI Low" = c(Individual_Amplitude_Model_CVR$ci.lb), 
                                           "CI High" = c(Individual_Amplitude_Model_CVR$ci.ub), 
                                           "df" = c(Individual_Amplitude_Model_CVR$ddf[[1]]), 
                                           "p-value" = c(Individual_Amplitude_Model_CVR$pval))

Raw_Individual_Fluctuation_Type_CVR <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating", 
                                                                         "Stepwise", "Stochastic"),
                                                  "Studies" = c(individual_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], individual_fluctuation_study["Alternating", "Study"], 
                                                                individual_fluctuation_study["Stepwise", "Study"], individual_fluctuation_study["Stochastic", "Study"]), 
                                                  "Species" = c(individual_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], individual_fluctuation_table["Alternating", "group_no"], 
                                                                individual_fluctuation_table["Stepwise", "group_no"], individual_fluctuation_table["Stochastic", "group_no"]), 
                                                  "Effect Sizes" = c(individual_fluctuation_table["Sinusoidal (Sine Curve)", "K"], individual_fluctuation_table["Alternating", "K"], 
                                                                     individual_fluctuation_table["Stepwise", "K"], individual_fluctuation_table["Stochastic", "K"]),
                                                  "Estimate" = c(Individual_Fluctuation_Model_CVR$b[[2]], Individual_Fluctuation_Model_CVR$b[[1]],
                                                                 Individual_Fluctuation_Model_CVR$b[[3]], Individual_Fluctuation_Model_CVR$b[[4]]),
                                                  "CI Low" = c(Individual_Fluctuation_Model_CVR$ci.lb[2], Individual_Fluctuation_Model_CVR$ci.lb[1], 
                                                               Individual_Fluctuation_Model_CVR$ci.lb[3], Individual_Fluctuation_Model_CVR$ci.lb[4]), 
                                                  "CI High" = c(Individual_Fluctuation_Model_CVR$ci.ub[2], Individual_Fluctuation_Model_CVR$ci.ub[1], 
                                                                Individual_Fluctuation_Model_CVR$ci.ub[3], Individual_Fluctuation_Model_CVR$ci.ub[4]), 
                                                  "df" = c(Individual_Fluctuation_Model_CVR$ddf[[2]], Individual_Fluctuation_Model_CVR$ddf[[1]], 
                                                           Individual_Fluctuation_Model_CVR$ddf[[3]], Individual_Fluctuation_Model_CVR$ddf[[4]]), 
                                                  "p-value" = c(Individual_Fluctuation_Model_CVR$pval[2], Individual_Fluctuation_Model_CVR$pval[1], 
                                                                Individual_Fluctuation_Model_CVR$pval[3], Individual_Fluctuation_Model_CVR$pval[4]))

Raw_Individual_Class_CVR <- data.frame("Taxonomic Class" = c("Actinopteri", "Amphibia", "Arachnida", "Bivalvia", 
                                                             "Branchiopoda", "Gastropoda", "Holothuroidea", "Insecta", "Malacostraca"),
                                       "Studies" = c(individual_class_study["Actinopteri", "Study"], individual_class_study["Amphibia", "Study"], individual_class_study["Arachnida", "Study"], 
                                                     individual_class_study["Bivalvia", "Study"], individual_class_study["Branchiopoda", "Study"], individual_class_study["Gastropoda", "Study"], 
                                                 individual_class_study["Holothuroidea", "Study"], individual_class_study["Insecta", "Study"], individual_class_study["Malacostraca", "Study"]), 
                                       "Species" = c(individual_class_table["Actinopteri", "group_no"], individual_class_table["Amphibia", "group_no"], individual_class_table["Arachnida", "group_no"], 
                                                     individual_class_table["Bivalvia", "group_no"], individual_class_table["Branchiopoda", "group_no"], individual_class_table["Gastropoda", "group_no"], 
                                                     individual_class_table["Holothuroidea", "group_no"], individual_class_table["Insecta", "group_no"], individual_class_table["Malacostraca", "group_no"]), 
                                       "Effect Sizes" = c(individual_class_table["Actinopteri", "K"], individual_class_table["Amphibia", "K"], individual_class_table["Arachnida", "K"], 
                                                          individual_class_table["Bivalvia", "K"], individual_class_table["Branchiopoda", "K"], individual_class_table["Gastropoda", "K"], 
                                                          individual_class_table["Holothuroidea", "K"], individual_class_table["Insecta", "K"], individual_class_table["Malacostraca", "K"]),
                                       "Estimate" = c(Individual_Class_Model_CVR$b[[1]], Individual_Class_Model_CVR$b[[2]], Individual_Class_Model_CVR$b[[3]], 
                                                      Individual_Class_Model_CVR$b[[4]], Individual_Class_Model_CVR$b[[5]], Individual_Class_Model_CVR$b[[6]], 
                                                      Individual_Class_Model_CVR$b[[7]], Individual_Class_Model_CVR$b[[8]], Individual_Class_Model_CVR$b[[9]]),
                                       "CI Low" = c(Individual_Class_Model_CVR$ci.lb[1], Individual_Class_Model_CVR$ci.lb[2], Individual_Class_Model_CVR$ci.lb[3], 
                                                    Individual_Class_Model_CVR$ci.lb[4], Individual_Class_Model_CVR$ci.lb[5], Individual_Class_Model_CVR$ci.lb[6], 
                                                    Individual_Class_Model_CVR$ci.lb[7], Individual_Class_Model_CVR$ci.lb[8], Individual_Class_Model_CVR$ci.lb[9]), 
                                       "CI High" = c(Individual_Class_Model_CVR$ci.ub[1], Individual_Class_Model_CVR$ci.ub[2], Individual_Class_Model_CVR$ci.ub[3], 
                                                     Individual_Class_Model_CVR$ci.ub[4], Individual_Class_Model_CVR$ci.ub[5], Individual_Class_Model_CVR$ci.ub[6], 
                                                     Individual_Class_Model_CVR$ci.ub[7], Individual_Class_Model_CVR$ci.ub[8], Individual_Class_Model_CVR$ci.ub[9]), 
                                       "df" = c(Individual_Class_Model_CVR$ddf[[1]], Individual_Class_Model_CVR$ddf[[2]], Individual_Class_Model_CVR$ddf[[3]], 
                                                Individual_Class_Model_CVR$ddf[[4]], Individual_Class_Model_CVR$ddf[[5]], Individual_Class_Model_CVR$ddf[[6]], 
                                                Individual_Class_Model_CVR$ddf[[7]], Individual_Class_Model_CVR$ddf[[8]], Individual_Class_Model_CVR$ddf[[9]]), 
                                       "p-value" = c(Individual_Class_Model_CVR$pval[1], Individual_Class_Model_CVR$pval[2], Individual_Class_Model_CVR$pval[3], 
                                                     Individual_Class_Model_CVR$pval[4], Individual_Class_Model_CVR$pval[5], Individual_Class_Model_CVR$pval[6], 
                                                     Individual_Class_Model_CVR$pval[7], Individual_Class_Model_CVR$pval[8], Individual_Class_Model_CVR$pval[9]))

Raw_Aquatic_CVR <- data.frame("Aquatic Organisms" = c("MLMA"),
                              "Studies" = c(intercept_study["Aquatic", "Study"]), 
                              "Species" = c(intercept_table["Aquatic", "group_no"]), 
                              "Effect Sizes" = c(intercept_table["Aquatic", "K"]),
                              "Estimate" = c(Aquatic_Model_CVR$b[1]),
                              "CI Low" = c(Aquatic_Model_CVR$ci.lb), 
                              "CI High" = c(Aquatic_Model_CVR$ci.ub), 
                              "df" = c(Aquatic_Model_CVR$ddf[[1]]), 
                              "p-value" = c(Aquatic_Model_CVR$pval))

Raw_Aquatic_Amplitude_CVR <- data.frame("Aquatic Organisms" = c("Fluctuation Amplitude"),
                                        "Studies" = c(length(unique(Aquatic_Subset_Data$Study_ID))), 
                                        "Species" = c(length(unique(Aquatic_Subset_Data$Scientific_Name))), 
                                        "Effect Sizes" = c(length(Aquatic_Subset_Data$Effect_Size_ID)),
                                        "Estimate" = c(Aquatic_Amplitude_Model_CVR$b[1]),
                                        "CI Low" = c(Aquatic_Amplitude_Model_CVR$ci.lb), 
                                        "CI High" = c(Aquatic_Amplitude_Model_CVR$ci.ub), 
                                        "df" = c(Aquatic_Amplitude_Model_CVR$ddf[[1]]), 
                                        "p-value" = c(Aquatic_Amplitude_Model_CVR$pval))

Raw_Aquatic_Fluctuation_Type_CVR <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating"),
                                           "Studies" = c(aquatic_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], aquatic_fluctuation_study["Alternating", "Study"]), 
                                           "Species" = c(aquatic_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], aquatic_fluctuation_table["Alternating", "group_no"]), 
                                           "Effect Sizes" = c(aquatic_fluctuation_table["Sinusoidal (Sine Curve)", "K"], aquatic_fluctuation_table["Alternating", "K"]),
                                           "Estimate" = c(Aquatic_Fluctuation_Model_CVR$b[[2]], Aquatic_Fluctuation_Model_CVR$b[[1]]),
                                           "CI Low" = c(Aquatic_Fluctuation_Model_CVR$ci.lb[2], Aquatic_Fluctuation_Model_CVR$ci.lb[1]), 
                                           "CI High" = c(Aquatic_Fluctuation_Model_CVR$ci.ub[2], Aquatic_Fluctuation_Model_CVR$ci.ub[1]), 
                                           "df" = c(Aquatic_Fluctuation_Model_CVR$ddf[[2]], Aquatic_Fluctuation_Model_CVR$ddf[[1]]), 
                                           "p-value" = c(Aquatic_Fluctuation_Model_CVR$pval[2], Aquatic_Fluctuation_Model_CVR$pval[1]))

Raw_Aquatic_Trait_CVR <- data.frame("Phenotypic Trait Categories" = c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits",  
                                                                      "Morphology", "Physiological"),
                                    "Studies" = c(aquatic_trait_study["Behavioural", "Study"], aquatic_trait_study["Biochemical Assay", "Study"], aquatic_trait_study["Gene Expression", "Study"], aquatic_trait_study["Life-history Traits", "Study"],
                                                  aquatic_trait_study["Morphology", "Study"], aquatic_trait_study["Physiological", "Study"]), 
                                    "Species" = c(aquatic_trait_table["Behavioural", "group_no"], aquatic_trait_table["Biochemical Assay", "group_no"], aquatic_trait_table["Gene Expression", "group_no"], aquatic_trait_table["Life-history Traits", "group_no"],
                                                  aquatic_trait_table["Morphology", "group_no"], aquatic_trait_table["Physiological", "group_no"]), 
                                    "Effect Sizes" = c(aquatic_trait_table["Behavioural", "K"], aquatic_trait_table["Biochemical Assay", "K"], aquatic_trait_table["Gene Expression", "K"], aquatic_trait_table["Life-history Traits", "K"],
                                                       aquatic_trait_table["Morphology", "K"], aquatic_trait_table["Physiological", "K"]),
                                    "Estimate" = c(Aquatic_Trait_Model_CVR$b[[1]], Aquatic_Trait_Model_CVR$b[[2]], Aquatic_Trait_Model_CVR$b[[3]], Aquatic_Trait_Model_CVR$b[[4]],
                                                   Aquatic_Trait_Model_CVR$b[[5]], Aquatic_Trait_Model_CVR$b[[6]]),
                                    "CI Low" = c(Aquatic_Trait_Model_CVR$ci.lb[1], Aquatic_Trait_Model_CVR$ci.lb[2], Aquatic_Trait_Model_CVR$ci.lb[3], Aquatic_Trait_Model_CVR$ci.lb[4], 
                                                 Aquatic_Trait_Model_CVR$ci.lb[5], Aquatic_Trait_Model_CVR$ci.lb[6]), 
                                    "CI High" = c(Aquatic_Trait_Model_CVR$ci.ub[1], Aquatic_Trait_Model_CVR$ci.ub[2], Aquatic_Trait_Model_CVR$ci.ub[3], Aquatic_Trait_Model_CVR$ci.ub[4], 
                                                  Aquatic_Trait_Model_CVR$ci.ub[5], Aquatic_Trait_Model_CVR$ci.ub[6]), 
                                    "df" = c(Aquatic_Trait_Model_CVR$ddf[[1]], Aquatic_Trait_Model_CVR$ddf[[2]], Aquatic_Trait_Model_CVR$ddf[[3]], Aquatic_Trait_Model_CVR$ddf[[4]], 
                                             Aquatic_Trait_Model_CVR$ddf[[5]], Aquatic_Trait_Model_CVR$ddf[[6]]), 
                                    "p-value" = c(Aquatic_Trait_Model_CVR$pval[1], Aquatic_Trait_Model_CVR$pval[2], Aquatic_Trait_Model_CVR$pval[3], Aquatic_Trait_Model_CVR$pval[4], 
                                                  Aquatic_Trait_Model_CVR$pval[5], Aquatic_Trait_Model_CVR$pval[6]))

Raw_Aquatic_Plasticity_CVR <- data.frame("Exposure Type" = c("Acclimation", "Developmental"),
                                         "Studies" = c(aquatic_plasticity_study["Acclimation", "Study"], aquatic_plasticity_study["Development", "Study"]), 
                                         "Species" = c(aquatic_plasticity_table["Acclimation", "group_no"], aquatic_plasticity_table["Development", "group_no"]), 
                                         "Effect Sizes" = c(aquatic_plasticity_table["Acclimation", "K"], aquatic_plasticity_table["Development", "K"]),
                                         "Estimate" = c(Aquatic_Plasticity_Model_CVR$b[[1]], Aquatic_Plasticity_Model_CVR$b[[2]]),
                                         "CI Low" = c(Aquatic_Plasticity_Model_CVR$ci.lb[1], Aquatic_Plasticity_Model_CVR$ci.lb[2]), 
                                         "CI High" = c(Aquatic_Plasticity_Model_CVR$ci.ub[1], Aquatic_Plasticity_Model_CVR$ci.ub[2]), 
                                         "df" = c(Aquatic_Plasticity_Model_CVR$ddf[[1]], Aquatic_Plasticity_Model_CVR$ddf[[2]]), 
                                         "p-value" = c(Aquatic_Plasticity_Model_CVR$pval[1], Aquatic_Plasticity_Model_CVR$pval[2]))

Raw_Aquatic_Specific_Trait_CVR <- data.frame("Specific Phenotypic Traits" = c("Apparent Digestibility Coefficient", "Cortisol Levels", "Development Time", 
                                                                              "Food Consumption", "Immune Defense", "Length", 
                                                                              "Locomotor Performance", "Mass", "Metabolic Rate"),
                                             "Studies" = c(aquatic_specific_trait_study["Apparent Digestability Coefficient", "Study"], aquatic_specific_trait_study["Cortisol", "Study"], 
                                                           aquatic_specific_trait_study["Development Time", "Study"], aquatic_specific_trait_study["Food Consumption", "Study"], 
                                                           aquatic_specific_trait_study["Immune Defense", "Study"], aquatic_specific_trait_study["Length", "Study"], 
                                                           aquatic_specific_trait_study["Locomotor Performance", "Study"], aquatic_specific_trait_study["Mass", "Study"], 
                                                           aquatic_specific_trait_study["Metabolic Rate", "Study"]), 
                                             "Species" = c(aquatic_specific_trait_table["Apparent Digestability Coefficient", "group_no"], aquatic_specific_trait_table["Cortisol", "group_no"], 
                                                           aquatic_specific_trait_table["Development Time", "group_no"], aquatic_specific_trait_table["Food Consumption", "group_no"], 
                                                           aquatic_specific_trait_table["Immune Defense", "group_no"], aquatic_specific_trait_table["Length", "group_no"], 
                                                           aquatic_specific_trait_table["Locomotor Performance", "group_no"], aquatic_specific_trait_table["Mass", "group_no"], 
                                                           aquatic_specific_trait_table["Metabolic Rate", "group_no"]), 
                                             "Effect Sizes" = c(aquatic_specific_trait_table["Apparent Digestability Coefficient", "K"], aquatic_specific_trait_table["Cortisol", "K"], 
                                                                aquatic_specific_trait_table["Development Time", "K"], aquatic_specific_trait_table["Food Consumption", "K"], 
                                                                aquatic_specific_trait_table["Immune Defense", "K"], aquatic_specific_trait_table["Length", "K"], 
                                                                aquatic_specific_trait_table["Locomotor Performance", "K"], aquatic_specific_trait_table["Mass", "K"], 
                                                                aquatic_specific_trait_table["Metabolic Rate", "K"]),
                                             "Estimate" = c(Aquatic_Specific_Trait_Model_CVR$b[[1]], Aquatic_Specific_Trait_Model_CVR$b[[2]], Aquatic_Specific_Trait_Model_CVR$b[[3]], 
                                                            Aquatic_Specific_Trait_Model_CVR$b[[4]], Aquatic_Specific_Trait_Model_CVR$b[[5]], Aquatic_Specific_Trait_Model_CVR$b[[6]], 
                                                            Aquatic_Specific_Trait_Model_CVR$b[[7]], Aquatic_Specific_Trait_Model_CVR$b[[8]], Aquatic_Specific_Trait_Model_CVR$b[[9]]),
                                             "CI Low" = c(Aquatic_Specific_Trait_Model_CVR$ci.lb[1], Aquatic_Specific_Trait_Model_CVR$ci.lb[2], Aquatic_Specific_Trait_Model_CVR$ci.lb[3], 
                                                          Aquatic_Specific_Trait_Model_CVR$ci.lb[4], Aquatic_Specific_Trait_Model_CVR$ci.lb[5], Aquatic_Specific_Trait_Model_CVR$ci.lb[6], 
                                                          Aquatic_Specific_Trait_Model_CVR$ci.lb[7], Aquatic_Specific_Trait_Model_CVR$ci.lb[8], Aquatic_Specific_Trait_Model_CVR$ci.lb[9]), 
                                             "CI High" = c(Aquatic_Specific_Trait_Model_CVR$ci.ub[1], Aquatic_Specific_Trait_Model_CVR$ci.ub[2], Aquatic_Specific_Trait_Model_CVR$ci.ub[3], 
                                                           Aquatic_Specific_Trait_Model_CVR$ci.ub[4], Aquatic_Specific_Trait_Model_CVR$ci.ub[5], Aquatic_Specific_Trait_Model_CVR$ci.ub[6], 
                                                           Aquatic_Specific_Trait_Model_CVR$ci.ub[7], Aquatic_Specific_Trait_Model_CVR$ci.ub[8], Aquatic_Specific_Trait_Model_CVR$ci.ub[9]), 
                                             "df" = c(Aquatic_Specific_Trait_Model_CVR$ddf[[1]], Aquatic_Specific_Trait_Model_CVR$ddf[[2]], Aquatic_Specific_Trait_Model_CVR$ddf[[3]], 
                                                      Aquatic_Specific_Trait_Model_CVR$ddf[[4]], Aquatic_Specific_Trait_Model_CVR$ddf[[5]], Aquatic_Specific_Trait_Model_CVR$ddf[[6]], 
                                                      Aquatic_Specific_Trait_Model_CVR$ddf[[7]], Aquatic_Specific_Trait_Model_CVR$ddf[[8]], Aquatic_Specific_Trait_Model_CVR$ddf[[9]]), 
                                             "p-value" = c(Aquatic_Specific_Trait_Model_CVR$pval[1], Aquatic_Specific_Trait_Model_CVR$pval[2], Aquatic_Specific_Trait_Model_CVR$pval[3], 
                                                           Aquatic_Specific_Trait_Model_CVR$pval[4], Aquatic_Specific_Trait_Model_CVR$pval[5], Aquatic_Specific_Trait_Model_CVR$pval[6], 
                                                           Aquatic_Specific_Trait_Model_CVR$pval[7], Aquatic_Specific_Trait_Model_CVR$pval[8], Aquatic_Specific_Trait_Model_CVR$pval[9]))

Raw_Terrestrial_CVR <- data.frame("Terrestrial Organisms" = c("MLMA"),
                                  "Studies" = c(intercept_study["Terrestrial", "Study"]), 
                                  "Species" = c(intercept_table["Terrestrial", "group_no"]), 
                                  "Effect Sizes" = c(intercept_table["Terrestrial", "K"]),
                                  "Estimate" = c(Terrestrial_Model_CVR$b[1]),
                                  "CI Low" = c(Terrestrial_Model_CVR$ci.lb), 
                                  "CI High" = c(Terrestrial_Model_CVR$ci.ub), 
                                  "df" = c(Terrestrial_Model_CVR$ddf[[1]]), 
                                  "p-value" = c(Terrestrial_Model_CVR$pval))

Raw_Terrestrial_Amplitude_CVR <- data.frame("Terrestrial Organisms" = c("Fluctuation Amplitude"),
                                            "Studies" = c(length(unique(Terrestrial_Subset_Data$Study_ID))), 
                                            "Species" = c(length(unique(Terrestrial_Subset_Data$Scientific_Name))), 
                                            "Effect Sizes" = c(length(Terrestrial_Subset_Data$Effect_Size_ID)),
                                            "Estimate" = c(Terrestrial_Amplitude_Model_CVR$b[1]),
                                            "CI Low" = c(Terrestrial_Amplitude_Model_CVR$ci.lb), 
                                            "CI High" = c(Terrestrial_Amplitude_Model_CVR$ci.ub), 
                                            "df" = c(Terrestrial_Amplitude_Model_CVR$ddf[[1]]), 
                                            "p-value" = c(Terrestrial_Amplitude_Model_CVR$pval))

Raw_Terrestrial_Fluctuation_Type_CVR <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating", 
                                                                          "Stepwise", "Stochastic"),
                                                   "Studies" = c(terrestrial_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], terrestrial_fluctuation_study["Alternating", "Study"], 
                                                                 terrestrial_fluctuation_study["Stepwise", "Study"], terrestrial_fluctuation_study["Stochastic", "Study"]), 
                                                   "Species" = c(terrestrial_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], terrestrial_fluctuation_table["Alternating", "group_no"], 
                                                                 terrestrial_fluctuation_table["Stepwise", "group_no"], terrestrial_fluctuation_table["Stochastic", "group_no"]), 
                                                   "Effect Sizes" = c(terrestrial_fluctuation_table["Sinusoidal (Sine Curve)", "K"], terrestrial_fluctuation_table["Alternating", "K"], 
                                                                      terrestrial_fluctuation_table["Stepwise", "K"], terrestrial_fluctuation_table["Stochastic", "K"]),
                                                   "Estimate" = c(Terrestrial_Fluctuation_Model_CVR$b[[2]], Terrestrial_Fluctuation_Model_CVR$b[[1]],
                                                                  Terrestrial_Fluctuation_Model_CVR$b[[3]], Terrestrial_Fluctuation_Model_CVR$b[[4]]),
                                                   "CI Low" = c(Terrestrial_Fluctuation_Model_CVR$ci.lb[2], Terrestrial_Fluctuation_Model_CVR$ci.lb[1], 
                                                                Terrestrial_Fluctuation_Model_CVR$ci.lb[3], Terrestrial_Fluctuation_Model_CVR$ci.lb[4]), 
                                                   "CI High" = c(Terrestrial_Fluctuation_Model_CVR$ci.ub[2], Terrestrial_Fluctuation_Model_CVR$ci.ub[1], 
                                                                 Terrestrial_Fluctuation_Model_CVR$ci.ub[3], Terrestrial_Fluctuation_Model_CVR$ci.ub[4]), 
                                                   "df" = c(Terrestrial_Fluctuation_Model_CVR$ddf[[2]], Terrestrial_Fluctuation_Model_CVR$ddf[[1]], 
                                                            Terrestrial_Fluctuation_Model_CVR$ddf[[3]], Terrestrial_Fluctuation_Model_CVR$ddf[[4]]), 
                                                   "p-value" = c(Terrestrial_Fluctuation_Model_CVR$pval[2], Terrestrial_Fluctuation_Model_CVR$pval[1], 
                                                                 Terrestrial_Fluctuation_Model_CVR$pval[3], Terrestrial_Fluctuation_Model_CVR$pval[4]))

Raw_Terrestrial_Trait_CVR <- data.frame("Phenotypic Trait Categories" = c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits",  
                                                                          "Morphology", "Physiological"),
                                        "Studies" = c(terrestrial_trait_study["Behavioural", "Study"], terrestrial_trait_study["Biochemical Assay", "Study"], terrestrial_trait_study["Gene Expression", "Study"], terrestrial_trait_study["Life-history Traits", "Study"],
                                                      terrestrial_trait_study["Morphology", "Study"], terrestrial_trait_study["Physiological", "Study"]), 
                                        "Species" = c(terrestrial_trait_table["Behavioural", "group_no"], terrestrial_trait_table["Biochemical Assay", "group_no"], terrestrial_trait_table["Gene Expression", "group_no"], terrestrial_trait_table["Life-history Traits", "group_no"],
                                                      terrestrial_trait_table["Morphology", "group_no"], terrestrial_trait_table["Physiological", "group_no"]), 
                                        "Effect Sizes" = c(terrestrial_trait_table["Behavioural", "K"], terrestrial_trait_table["Biochemical Assay", "K"], terrestrial_trait_table["Gene Expression", "K"], terrestrial_trait_table["Life-history Traits", "K"],
                                                           terrestrial_trait_table["Morphology", "K"], terrestrial_trait_table["Physiological", "K"]),
                                        "Estimate" = c(Terrestrial_Trait_Model_CVR$b[[1]], Terrestrial_Trait_Model_CVR$b[[2]], Terrestrial_Trait_Model_CVR$b[[3]], Terrestrial_Trait_Model_CVR$b[[4]],
                                                       Terrestrial_Trait_Model_CVR$b[[5]], Terrestrial_Trait_Model_CVR$b[[6]]),
                                        "CI Low" = c(Terrestrial_Trait_Model_CVR$ci.lb[1], Terrestrial_Trait_Model_CVR$ci.lb[2], Terrestrial_Trait_Model_CVR$ci.lb[3], Terrestrial_Trait_Model_CVR$ci.lb[4], 
                                                     Terrestrial_Trait_Model_CVR$ci.lb[5], Terrestrial_Trait_Model_CVR$ci.lb[6]), 
                                        "CI High" = c(Terrestrial_Trait_Model_CVR$ci.ub[1], Terrestrial_Trait_Model_CVR$ci.ub[2], Terrestrial_Trait_Model_CVR$ci.ub[3], Terrestrial_Trait_Model_CVR$ci.ub[4], 
                                                      Terrestrial_Trait_Model_CVR$ci.ub[5], Terrestrial_Trait_Model_CVR$ci.ub[6]), 
                                        "df" = c(Terrestrial_Trait_Model_CVR$ddf[[1]], Terrestrial_Trait_Model_CVR$ddf[[2]], Terrestrial_Trait_Model_CVR$ddf[[3]], Terrestrial_Trait_Model_CVR$ddf[[4]], 
                                                 Terrestrial_Trait_Model_CVR$ddf[[5]], Terrestrial_Trait_Model_CVR$ddf[[6]]), 
                                        "p-value" = c(Terrestrial_Trait_Model_CVR$pval[1], Terrestrial_Trait_Model_CVR$pval[2], Terrestrial_Trait_Model_CVR$pval[3], Terrestrial_Trait_Model_CVR$pval[4], 
                                                      Terrestrial_Trait_Model_CVR$pval[5], Terrestrial_Trait_Model_CVR$pval[6]))

Raw_Terrestrial_Plasticity_CVR <- data.frame("Exposure Type" = c("Acclimation", "Developmental"),
                                             "Studies" = c(terrestrial_plasticity_study["Acclimation", "Study"], terrestrial_plasticity_study["Development", "Study"]), 
                                             "Species" = c(terrestrial_plasticity_table["Acclimation", "group_no"], terrestrial_plasticity_table["Development", "group_no"]), 
                                             "Effect Sizes" = c(terrestrial_plasticity_table["Acclimation", "K"], terrestrial_plasticity_table["Development", "K"]),
                                             "Estimate" = c(Terrestrial_Plasticity_Model_CVR$b[[1]], Terrestrial_Plasticity_Model_CVR$b[[2]]),
                                             "CI Low" = c(Terrestrial_Plasticity_Model_CVR$ci.lb[1], Terrestrial_Plasticity_Model_CVR$ci.lb[2]), 
                                             "CI High" = c(Terrestrial_Plasticity_Model_CVR$ci.ub[1], Terrestrial_Plasticity_Model_CVR$ci.ub[2]), 
                                             "df" = c(Terrestrial_Plasticity_Model_CVR$ddf[[1]], Terrestrial_Plasticity_Model_CVR$ddf[[2]]), 
                                             "p-value" = c(Terrestrial_Plasticity_Model_CVR$pval[1], Terrestrial_Plasticity_Model_CVR$pval[2]))

Raw_Terrestrial_Specific_Trait_CVR <- data.frame("Specific Phenotypic Traits" = c("Development Time", "Fecundity", "Food Consumption", "Head Width", 
                                                                                  "Length", "Locomotor Performance", "Longevity", "Mass", 
                                                                                  "Metabolic Rate", "PO Activity", "Reproductive Rate", "Tail Length"),
                                                 "Studies" = c(terrestrial_specific_trait_study["Development Time", "Study"], terrestrial_specific_trait_study["Fecundity", "Study"], terrestrial_specific_trait_study["Food Consumption", "Study"], 
                                                               terrestrial_specific_trait_study["Head Width", "Study"], terrestrial_specific_trait_study["Length", "Study"], terrestrial_specific_trait_study["Locomotor Performance", "Study"], 
                                                               terrestrial_specific_trait_study["Longevity", "Study"], terrestrial_specific_trait_study["Mass", "Study"], terrestrial_specific_trait_study["Metabolic Rate", "Study"], 
                                                               terrestrial_specific_trait_study["PO Activity", "Study"], terrestrial_specific_trait_study["Reproductive Rate", "Study"], terrestrial_specific_trait_study["Tail Length", "Study"]), 
                                                 "Species" = c(terrestrial_specific_trait_table["Development Time", "group_no"], terrestrial_specific_trait_table["Fecundity", "group_no"], terrestrial_specific_trait_table["Food Consumption", "group_no"], 
                                                               terrestrial_specific_trait_table["Head Width", "group_no"], terrestrial_specific_trait_table["Length", "group_no"], terrestrial_specific_trait_table["Locomotor Performance", "group_no"], 
                                                               terrestrial_specific_trait_table["Longevity", "group_no"], terrestrial_specific_trait_table["Mass", "group_no"], terrestrial_specific_trait_table["Metabolic Rate", "group_no"], 
                                                               terrestrial_specific_trait_table["PO Activity", "group_no"], terrestrial_specific_trait_table["Reproductive Rate", "group_no"], terrestrial_specific_trait_table["Tail Length", "group_no"]), 
                                                 "Effect Sizes" = c(terrestrial_specific_trait_table["Development Time", "K"], terrestrial_specific_trait_table["Fecundity", "K"], terrestrial_specific_trait_table["Food Consumption", "K"], 
                                                                    terrestrial_specific_trait_table["Head Width", "K"], terrestrial_specific_trait_table["Length", "K"], terrestrial_specific_trait_table["Locomotor Performance", "K"], 
                                                                    terrestrial_specific_trait_table["Longevity", "K"], terrestrial_specific_trait_table["Mass", "K"], terrestrial_specific_trait_table["Metabolic Rate", "K"], 
                                                                    terrestrial_specific_trait_table["PO Activity", "K"], terrestrial_specific_trait_table["Reproductive Rate", "K"], terrestrial_specific_trait_table["Tail Length", "K"]),
                                                 "Estimate" = c(Terrestrial_Specific_Trait_Model_CVR$b[[1]], Terrestrial_Specific_Trait_Model_CVR$b[[2]], Terrestrial_Specific_Trait_Model_CVR$b[[3]], 
                                                                Terrestrial_Specific_Trait_Model_CVR$b[[4]], Terrestrial_Specific_Trait_Model_CVR$b[[5]], Terrestrial_Specific_Trait_Model_CVR$b[[6]], 
                                                                Terrestrial_Specific_Trait_Model_CVR$b[[7]], Terrestrial_Specific_Trait_Model_CVR$b[[8]], Terrestrial_Specific_Trait_Model_CVR$b[[9]], 
                                                                Terrestrial_Specific_Trait_Model_CVR$b[[10]], Terrestrial_Specific_Trait_Model_CVR$b[[11]], Terrestrial_Specific_Trait_Model_CVR$b[[12]]),
                                                 "CI Low" = c(Terrestrial_Specific_Trait_Model_CVR$ci.lb[1], Terrestrial_Specific_Trait_Model_CVR$ci.lb[2], Terrestrial_Specific_Trait_Model_CVR$ci.lb[3], 
                                                              Terrestrial_Specific_Trait_Model_CVR$ci.lb[4], Terrestrial_Specific_Trait_Model_CVR$ci.lb[5], Terrestrial_Specific_Trait_Model_CVR$ci.lb[6], 
                                                              Terrestrial_Specific_Trait_Model_CVR$ci.lb[7], Terrestrial_Specific_Trait_Model_CVR$ci.lb[8], Terrestrial_Specific_Trait_Model_CVR$ci.lb[9], 
                                                              Terrestrial_Specific_Trait_Model_CVR$ci.lb[10], Terrestrial_Specific_Trait_Model_CVR$ci.lb[11], Terrestrial_Specific_Trait_Model_CVR$ci.lb[12]), 
                                                 "CI High" = c(Terrestrial_Specific_Trait_Model_CVR$ci.ub[1], Terrestrial_Specific_Trait_Model_CVR$ci.ub[2], Terrestrial_Specific_Trait_Model_CVR$ci.ub[3], 
                                                               Terrestrial_Specific_Trait_Model_CVR$ci.ub[4], Terrestrial_Specific_Trait_Model_CVR$ci.ub[5], Terrestrial_Specific_Trait_Model_CVR$ci.ub[6], 
                                                               Terrestrial_Specific_Trait_Model_CVR$ci.ub[7], Terrestrial_Specific_Trait_Model_CVR$ci.ub[8], Terrestrial_Specific_Trait_Model_CVR$ci.ub[9], 
                                                               Terrestrial_Specific_Trait_Model_CVR$ci.ub[10], Terrestrial_Specific_Trait_Model_CVR$ci.ub[11], Terrestrial_Specific_Trait_Model_CVR$ci.ub[12]), 
                                                 "df" = c(Terrestrial_Specific_Trait_Model_CVR$ddf[[1]], Terrestrial_Specific_Trait_Model_CVR$ddf[[2]], Terrestrial_Specific_Trait_Model_CVR$ddf[[3]], 
                                                          Terrestrial_Specific_Trait_Model_CVR$ddf[[4]], Terrestrial_Specific_Trait_Model_CVR$ddf[[5]], Terrestrial_Specific_Trait_Model_CVR$ddf[[6]], 
                                                          Terrestrial_Specific_Trait_Model_CVR$ddf[[7]], Terrestrial_Specific_Trait_Model_CVR$ddf[[8]], Terrestrial_Specific_Trait_Model_CVR$ddf[[9]], 
                                                          Terrestrial_Specific_Trait_Model_CVR$ddf[[10]], Terrestrial_Specific_Trait_Model_CVR$ddf[[11]], Terrestrial_Specific_Trait_Model_CVR$ddf[[12]]), 
                                                 "p-value" = c(Terrestrial_Specific_Trait_Model_CVR$pval[1], Terrestrial_Specific_Trait_Model_CVR$pval[2], Terrestrial_Specific_Trait_Model_CVR$pval[3], 
                                                               Terrestrial_Specific_Trait_Model_CVR$pval[4], Terrestrial_Specific_Trait_Model_CVR$pval[5], Terrestrial_Specific_Trait_Model_CVR$pval[6], 
                                                               Terrestrial_Specific_Trait_Model_CVR$pval[7], Terrestrial_Specific_Trait_Model_CVR$pval[8], Terrestrial_Specific_Trait_Model_CVR$pval[9], 
                                                               Terrestrial_Specific_Trait_Model_CVR$pval[10], Terrestrial_Specific_Trait_Model_CVR$pval[11], Terrestrial_Specific_Trait_Model_CVR$pval[12]))

Raw_Acclimation_CVR <- data.frame("Acclimation Exposure" = c("MLMA"),
                                  "Studies" = c(intercept_study["Acclimation", "Study"]), 
                                  "Species" = c(intercept_table["Acclimation", "group_no"]), 
                                  "Effect Sizes" = c(intercept_table["Acclimation", "K"]),
                                  "Estimate" = c(Acclimation_Model_CVR$b[1]),
                                  "CI Low" = c(Acclimation_Model_CVR$ci.lb), 
                                  "CI High" = c(Acclimation_Model_CVR$ci.ub), 
                                  "df" = c(Acclimation_Model_CVR$ddf[[1]]), 
                                  "p-value" = c(Acclimation_Model_CVR$pval))

Raw_Acclimation_Amplitude_CVR <- data.frame("Acclimation Exposure" = c("Fluctuation Amplitude"),
                                            "Studies" = c(length(unique(Acclimation_Subset_Data$Study_ID))), 
                                            "Species" = c(length(unique(Acclimation_Subset_Data$Scientific_Name))), 
                                            "Effect Sizes" = c(length(Acclimation_Subset_Data$Effect_Size_ID)),
                                            "Estimate" = c(Acclimation_Amplitude_Model_CVR$b[1]),
                                            "CI Low" = c(Acclimation_Amplitude_Model_CVR$ci.lb), 
                                            "CI High" = c(Acclimation_Amplitude_Model_CVR$ci.ub), 
                                            "df" = c(Acclimation_Amplitude_Model_CVR$ddf[[1]]), 
                                            "p-value" = c(Acclimation_Amplitude_Model_CVR$pval))

Raw_Acclimation_Exposure_CVR <- data.frame("Acclimation Exposure" = c("Exposure Time"),
                                           "Studies" = c(length(unique(Acclimation_Subset_Data$Study_ID))), 
                                           "Species" = c(length(unique(Acclimation_Subset_Data$Scientific_Name))), 
                                           "Effect Sizes" = c(length(Acclimation_Subset_Data$Effect_Size_ID)),
                                           "Estimate" = c(Acclimation_Exposure_Model_CVR$b[1]),
                                           "CI Low" = c(Acclimation_Exposure_Model_CVR$ci.lb), 
                                           "CI High" = c(Acclimation_Exposure_Model_CVR$ci.ub), 
                                           "df" = c(Acclimation_Exposure_Model_CVR$ddf[[1]]), 
                                           "p-value" = c(Acclimation_Exposure_Model_CVR$pval))

Raw_Acclimation_Frequency_CVR <- data.frame("Acclimation Exposure" = c("Number of Fluctuations"),
                                            "Studies" = c(length(unique(Acclimation_Frequency_Data$Study_ID))), 
                                            "Species" = c(length(unique(Acclimation_Frequency_Data$Scientific_Name))), 
                                            "Effect Sizes" = c(length(Acclimation_Frequency_Data$Effect_Size_ID)),
                                            "Estimate" = c(Acclimation_Frequency_Model_CVR$b[1]),
                                            "CI Low" = c(Acclimation_Frequency_Model_CVR$ci.lb), 
                                            "CI High" = c(Acclimation_Frequency_Model_CVR$ci.ub), 
                                            "df" = c(Acclimation_Frequency_Model_CVR$ddf[[1]]), 
                                            "p-value" = c(Acclimation_Frequency_Model_CVR$pval))

Raw_Acclimation_Fluctuation_Type_CVR <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating", 
                                                                          "Stepwise"),
                                                   "Studies" = c(acclimation_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], acclimation_fluctuation_study["Alternating", "Study"], 
                                                                 acclimation_fluctuation_study["Stepwise", "Study"]), 
                                                   "Species" = c(acclimation_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], acclimation_fluctuation_table["Alternating", "group_no"], 
                                                                 acclimation_fluctuation_table["Stepwise", "group_no"]), 
                                                   "Effect Sizes" = c(acclimation_fluctuation_table["Sinusoidal (Sine Curve)", "K"], acclimation_fluctuation_table["Alternating", "K"], 
                                                                      acclimation_fluctuation_table["Stepwise", "K"]),
                                                   "Estimate" = c(Acclimation_Fluctuation_Model_CVR$b[[2]], Acclimation_Fluctuation_Model_CVR$b[[1]],
                                                                  Acclimation_Fluctuation_Model_CVR$b[[3]]),
                                                   "CI Low" = c(Acclimation_Fluctuation_Model_CVR$ci.lb[2], Acclimation_Fluctuation_Model_CVR$ci.lb[1], 
                                                                Acclimation_Fluctuation_Model_CVR$ci.lb[3]), 
                                                   "CI High" = c(Acclimation_Fluctuation_Model_CVR$ci.ub[2], Acclimation_Fluctuation_Model_CVR$ci.ub[1], 
                                                                 Acclimation_Fluctuation_Model_CVR$ci.ub[3]), 
                                                   "df" = c(Acclimation_Fluctuation_Model_CVR$ddf[[2]], Acclimation_Fluctuation_Model_CVR$ddf[[1]], 
                                                            Acclimation_Fluctuation_Model_CVR$ddf[[3]]), 
                                                   "p-value" = c(Acclimation_Fluctuation_Model_CVR$pval[2], Acclimation_Fluctuation_Model_CVR$pval[1], 
                                                                 Acclimation_Fluctuation_Model_CVR$pval[3]))

Raw_Acclimation_Trait_CVR <- data.frame("Phenotypic Trait Categories" = c("Behavioural", "Biochemical Assay", "Life-history Traits", "Physiological"),
                                        "Studies" = c(acclimation_trait_study["Behavioural", "Study"], acclimation_trait_study["Biochemical Assay", "Study"], 
                                                      acclimation_trait_study["Life-history Traits", "Study"], acclimation_trait_study["Physiological", "Study"]), 
                                        "Species" = c(acclimation_trait_table["Behavioural", "group_no"], acclimation_trait_table["Biochemical Assay", "group_no"], 
                                                      acclimation_trait_table["Life-history Traits", "group_no"], acclimation_trait_table["Physiological", "group_no"]), 
                                        "Effect Sizes" = c(acclimation_trait_table["Behavioural", "K"], acclimation_trait_table["Biochemical Assay", "K"], 
                                                           acclimation_trait_table["Life-history Traits", "K"], acclimation_trait_table["Physiological", "K"]),
                                        "Estimate" = c(Acclimation_Trait_Model_CVR$b[[1]], Acclimation_Trait_Model_CVR$b[[2]], Acclimation_Trait_Model_CVR$b[[3]], Acclimation_Trait_Model_CVR$b[[4]]),
                                        "CI Low" = c(Acclimation_Trait_Model_CVR$ci.lb[1], Acclimation_Trait_Model_CVR$ci.lb[2], Acclimation_Trait_Model_CVR$ci.lb[3], Acclimation_Trait_Model_CVR$ci.lb[4]), 
                                        "CI High" = c(Acclimation_Trait_Model_CVR$ci.ub[1], Acclimation_Trait_Model_CVR$ci.ub[2], Acclimation_Trait_Model_CVR$ci.ub[3], Acclimation_Trait_Model_CVR$ci.ub[4]), 
                                        "df" = c(Acclimation_Trait_Model_CVR$ddf[[1]], Acclimation_Trait_Model_CVR$ddf[[2]], Acclimation_Trait_Model_CVR$ddf[[3]], Acclimation_Trait_Model_CVR$ddf[[4]]), 
                                        "p-value" = c(Acclimation_Trait_Model_CVR$pval[1], Acclimation_Trait_Model_CVR$pval[2], Acclimation_Trait_Model_CVR$pval[3], Acclimation_Trait_Model_CVR$pval[4]))

Raw_Acclimation_Stage_CVR <- data.frame("Life-history Stages" = c("Adult", "Embryo", "Juvenile", "Larva"),
                                        "Studies" = c(acclimation_stage_study["Adult", "Study"], acclimation_stage_study["Embryo", "Study"], 
                                                      acclimation_stage_study["Juvenile", "Study"], acclimation_stage_study["Larva", "Study"]), 
                                        "Species" = c(acclimation_stage_table["Adult", "group_no"], acclimation_stage_table["Embryo", "group_no"], 
                                                      acclimation_stage_table["Juvenile", "group_no"], acclimation_stage_table["Larva", "group_no"]), 
                                        "Effect Sizes" = c(acclimation_stage_table["Adult", "K"], acclimation_stage_table["Embryo", "K"], 
                                                           acclimation_stage_table["Juvenile", "K"], acclimation_stage_table["Larva", "K"]),
                                        "Estimate" = c(Acclimation_Stage_Model_CVR$b[[1]], Acclimation_Stage_Model_CVR$b[[2]], Acclimation_Stage_Model_CVR$b[[3]], Acclimation_Stage_Model_CVR$b[[4]]),
                                        "CI Low" = c(Acclimation_Stage_Model_CVR$ci.lb[1], Acclimation_Stage_Model_CVR$ci.lb[2], Acclimation_Stage_Model_CVR$ci.lb[3], Acclimation_Stage_Model_CVR$ci.lb[4]), 
                                        "CI High" = c(Acclimation_Stage_Model_CVR$ci.ub[1], Acclimation_Stage_Model_CVR$ci.ub[2], Acclimation_Stage_Model_CVR$ci.ub[3], Acclimation_Stage_Model_CVR$ci.ub[4]), 
                                        "df" = c(Acclimation_Stage_Model_CVR$ddf[[1]], Acclimation_Stage_Model_CVR$ddf[[2]], Acclimation_Stage_Model_CVR$ddf[[3]], Acclimation_Stage_Model_CVR$ddf[[4]]), 
                                        "p-value" = c(Acclimation_Stage_Model_CVR$pval[1], Acclimation_Stage_Model_CVR$pval[2], Acclimation_Stage_Model_CVR$pval[3], Acclimation_Stage_Model_CVR$pval[4]))

Raw_Acclimation_Class_CVR <- data.frame("Taxonomic Class" = c("Actinopteri", "Bivalvia", "Gastropoda", "Holothuroidea", "Insecta", "Malacostraca"),
                                        "Studies" = c(acclimation_class_study["Actinopteri", "Study"], acclimation_class_study["Bivalvia", "Study"], acclimation_class_study["Gastropoda", "Study"], 
                                                      acclimation_class_study["Holothuroidea", "Study"], acclimation_class_study["Insecta", "Study"], acclimation_class_study["Malacostraca", "Study"]), 
                                        "Species" = c(acclimation_class_table["Actinopteri", "group_no"], acclimation_class_table["Bivalvia", "group_no"], acclimation_class_table["Gastropoda", "group_no"], 
                                                      acclimation_class_table["Holothuroidea", "group_no"], acclimation_class_table["Insecta", "group_no"], acclimation_class_table["Malacostraca", "group_no"]), 
                                        "Effect Sizes" = c(acclimation_class_table["Actinopteri", "K"], acclimation_class_table["Bivalvia", "K"], acclimation_class_table["Gastropoda", "K"], 
                                                           acclimation_class_table["Holothuroidea", "K"], acclimation_class_table["Insecta", "K"], acclimation_class_table["Malacostraca", "K"]),
                                        "Estimate" = c(Acclimation_Class_Model_CVR$b[[1]], Acclimation_Class_Model_CVR$b[[2]], Acclimation_Class_Model_CVR$b[[3]], 
                                                       Acclimation_Class_Model_CVR$b[[4]], Acclimation_Class_Model_CVR$b[[5]], Acclimation_Class_Model_CVR$b[[6]]),
                                        "CI Low" = c(Acclimation_Class_Model_CVR$ci.lb[1], Acclimation_Class_Model_CVR$ci.lb[2], Acclimation_Class_Model_CVR$ci.lb[3], 
                                                     Acclimation_Class_Model_CVR$ci.lb[4], Acclimation_Class_Model_CVR$ci.lb[5], Acclimation_Class_Model_CVR$ci.lb[6]), 
                                        "CI High" = c(Acclimation_Class_Model_CVR$ci.ub[1], Acclimation_Class_Model_CVR$ci.ub[2], Acclimation_Class_Model_CVR$ci.ub[3], 
                                                      Acclimation_Class_Model_CVR$ci.ub[4], Acclimation_Class_Model_CVR$ci.ub[5], Acclimation_Class_Model_CVR$ci.ub[6]), 
                                        "df" = c(Acclimation_Class_Model_CVR$ddf[[1]], Acclimation_Class_Model_CVR$ddf[[2]], Acclimation_Class_Model_CVR$ddf[[3]], 
                                                 Acclimation_Class_Model_CVR$ddf[[4]], Acclimation_Class_Model_CVR$ddf[[5]], Acclimation_Class_Model_CVR$ddf[[6]]), 
                                        "p-value" = c(Acclimation_Class_Model_CVR$pval[1], Acclimation_Class_Model_CVR$pval[2], Acclimation_Class_Model_CVR$pval[3], 
                                                      Acclimation_Class_Model_CVR$pval[4], Acclimation_Class_Model_CVR$pval[5], Acclimation_Class_Model_CVR$pval[6]))

Raw_Acclimation_Specific_Trait_CVR <- data.frame("Specific Phenotypic Traits" = c("Apparent Digestibility Coefficient", "Catalase Activity", "Cortisol Levels", "Food Consumption", 
                                                                                  "Immune Defense", "Metabolic Rate", "SOD Activity"),
                                                 "Studies" = c(acclimation_specific_trait_study["Apparent Digestability Coefficient", "Study"], acclimation_specific_trait_study["Catalase Activity", "Study"], acclimation_specific_trait_study["Cortisol", "Study"], 
                                                               acclimation_specific_trait_study["Food Consumption", "Study"], acclimation_specific_trait_study["Immune Defense", "Study"], acclimation_specific_trait_study["Metabolic Rate", "Study"], 
                                                               acclimation_specific_trait_study["SOD Activity", "Study"]), 
                                                 "Species" = c(acclimation_specific_trait_table["Apparent Digestability Coefficient", "group_no"], acclimation_specific_trait_table["Catalase Activity", "group_no"], acclimation_specific_trait_table["Cortisol", "group_no"], 
                                                               acclimation_specific_trait_table["Food Consumption", "group_no"], acclimation_specific_trait_table["Immune Defense", "group_no"], acclimation_specific_trait_table["Metabolic Rate", "group_no"], 
                                                               acclimation_specific_trait_table["SOD Activity", "group_no"]), 
                                                 "Effect Sizes" = c(acclimation_specific_trait_table["Apparent Digestability Coefficient", "K"], acclimation_specific_trait_table["Catalase Activity", "K"], acclimation_specific_trait_table["Cortisol", "K"], 
                                                                    acclimation_specific_trait_table["Food Consumption", "K"], acclimation_specific_trait_table["Immune Defense", "K"], acclimation_specific_trait_table["Metabolic Rate", "K"], 
                                                                    acclimation_specific_trait_table["SOD Activity", "K"]),
                                                 "Estimate" = c(Acclimation_Specific_Trait_Model_CVR$b[[1]], Acclimation_Specific_Trait_Model_CVR$b[[2]], Acclimation_Specific_Trait_Model_CVR$b[[3]], 
                                                                Acclimation_Specific_Trait_Model_CVR$b[[4]], Acclimation_Specific_Trait_Model_CVR$b[[5]], Acclimation_Specific_Trait_Model_CVR$b[[6]], 
                                                                Acclimation_Specific_Trait_Model_CVR$b[[7]]),
                                                 "CI Low" = c(Acclimation_Specific_Trait_Model_CVR$ci.lb[1], Acclimation_Specific_Trait_Model_CVR$ci.lb[2], Acclimation_Specific_Trait_Model_CVR$ci.lb[3], 
                                                              Acclimation_Specific_Trait_Model_CVR$ci.lb[4], Acclimation_Specific_Trait_Model_CVR$ci.lb[5], Acclimation_Specific_Trait_Model_CVR$ci.lb[6], 
                                                              Acclimation_Specific_Trait_Model_CVR$ci.lb[7]), 
                                                 "CI High" = c(Acclimation_Specific_Trait_Model_CVR$ci.ub[1], Acclimation_Specific_Trait_Model_CVR$ci.ub[2], Acclimation_Specific_Trait_Model_CVR$ci.ub[3], 
                                                               Acclimation_Specific_Trait_Model_CVR$ci.ub[4], Acclimation_Specific_Trait_Model_CVR$ci.ub[5], Acclimation_Specific_Trait_Model_CVR$ci.ub[6], 
                                                               Acclimation_Specific_Trait_Model_CVR$ci.ub[7]), 
                                                 "df" = c(Acclimation_Specific_Trait_Model_CVR$ddf[[1]], Acclimation_Specific_Trait_Model_CVR$ddf[[2]], Acclimation_Specific_Trait_Model_CVR$ddf[[3]], 
                                                          Acclimation_Specific_Trait_Model_CVR$ddf[[4]], Acclimation_Specific_Trait_Model_CVR$ddf[[5]], Acclimation_Specific_Trait_Model_CVR$ddf[[6]], 
                                                          Acclimation_Specific_Trait_Model_CVR$ddf[[7]]), 
                                                 "p-value" = c(Acclimation_Specific_Trait_Model_CVR$pval[1], Acclimation_Specific_Trait_Model_CVR$pval[2], Acclimation_Specific_Trait_Model_CVR$pval[3], 
                                                               Acclimation_Specific_Trait_Model_CVR$pval[4], Acclimation_Specific_Trait_Model_CVR$pval[5], Acclimation_Specific_Trait_Model_CVR$pval[6], 
                                                               Acclimation_Specific_Trait_Model_CVR$pval[7]))

Raw_Developmental_CVR <- data.frame("Developmental Exposure" = c("MLMA"),
                                    "Studies" = c(intercept_study["Developmental", "Study"]), 
                                    "Species" = c(intercept_table["Developmental", "group_no"]), 
                                    "Effect Sizes" = c(intercept_table["Developmental", "K"]),
                                    "Estimate" = c(Developmental_Model_CVR$b[1]),
                                    "CI Low" = c(Developmental_Model_CVR$ci.lb), 
                                    "CI High" = c(Developmental_Model_CVR$ci.ub), 
                                    "df" = c(Developmental_Model_CVR$ddf[[1]]), 
                                    "p-value" = c(Developmental_Model_CVR$pval))

Raw_Developmental_Amplitude_CVR <- data.frame("Developmental Exposure" = c("Fluctuation Amplitude"),
                                              "Studies" = c(length(unique(Developmental_Subset_Data$Study_ID))), 
                                              "Species" = c(length(unique(Developmental_Subset_Data$Scientific_Name))), 
                                              "Effect Sizes" = c(length(Developmental_Subset_Data$Effect_Size_ID)),
                                              "Estimate" = c(Developmental_Amplitude_Model_CVR$b[1]),
                                              "CI Low" = c(Developmental_Amplitude_Model_CVR$ci.lb), 
                                              "CI High" = c(Developmental_Amplitude_Model_CVR$ci.ub), 
                                              "df" = c(Developmental_Amplitude_Model_CVR$ddf[[1]]), 
                                              "p-value" = c(Developmental_Amplitude_Model_CVR$pval))

Raw_Developmental_Fluctuation_Type_CVR <- data.frame("Fluctuation Type" = c("Sinusoidal (Sine Curve)", "Alternating", 
                                                                            "Stepwise", "Stochastic"),
                                                     "Studies" = c(developmental_fluctuation_study["Sinusoidal (Sine Curve)", "Study"], developmental_fluctuation_study["Alternating", "Study"], 
                                                                   developmental_fluctuation_study["Stepwise", "Study"], developmental_fluctuation_study["Stochastic", "Study"]), 
                                                     "Species" = c(developmental_fluctuation_table["Sinusoidal (Sine Curve)", "group_no"], developmental_fluctuation_table["Alternating", "group_no"], 
                                                                   developmental_fluctuation_table["Stepwise", "group_no"], developmental_fluctuation_table["Stochastic", "group_no"]), 
                                                     "Effect Sizes" = c(developmental_fluctuation_table["Sinusoidal (Sine Curve)", "K"], developmental_fluctuation_table["Alternating", "K"], 
                                                                        developmental_fluctuation_table["Stepwise", "K"], developmental_fluctuation_table["Stochastic", "K"]),
                                                     "Estimate" = c(Developmental_Fluctuation_Model_CVR$b[[2]], Developmental_Fluctuation_Model_CVR$b[[1]],
                                                                    Developmental_Fluctuation_Model_CVR$b[[3]], Developmental_Fluctuation_Model_CVR$b[[4]]),
                                                     "CI Low" = c(Developmental_Fluctuation_Model_CVR$ci.lb[2], Developmental_Fluctuation_Model_CVR$ci.lb[1], 
                                                                  Developmental_Fluctuation_Model_CVR$ci.lb[3], Developmental_Fluctuation_Model_CVR$ci.lb[4]), 
                                                     "CI High" = c(Developmental_Fluctuation_Model_CVR$ci.ub[2], Developmental_Fluctuation_Model_CVR$ci.ub[1], 
                                                                   Developmental_Fluctuation_Model_CVR$ci.ub[3], Developmental_Fluctuation_Model_CVR$ci.ub[4]), 
                                                     "df" = c(Developmental_Fluctuation_Model_CVR$ddf[[2]], Developmental_Fluctuation_Model_CVR$ddf[[1]], 
                                                              Developmental_Fluctuation_Model_CVR$ddf[[3]], Developmental_Fluctuation_Model_CVR$ddf[[4]]), 
                                                     "p-value" = c(Developmental_Fluctuation_Model_CVR$pval[2], Developmental_Fluctuation_Model_CVR$pval[1], 
                                                                   Developmental_Fluctuation_Model_CVR$pval[3], Developmental_Fluctuation_Model_CVR$pval[4]))

Raw_Developmental_Trait_CVR <- data.frame("Phenotypic Trait Categories" = c("Behavioural", "Biochemical Assay", "Gene Expression", "Life-history Traits",  
                                                                            "Morphology", "Physiological"),
                                          "Studies" = c(developmental_trait_study["Behavioural", "Study"], developmental_trait_study["Biochemical Assay", "Study"], developmental_trait_study["Gene Expression", "Study"], developmental_trait_study["Life-history Traits", "Study"],
                                                        developmental_trait_study["Morphology", "Study"], developmental_trait_study["Physiological", "Study"]), 
                                          "Species" = c(developmental_trait_table["Behavioural", "group_no"], developmental_trait_table["Biochemical Assay", "group_no"], developmental_trait_table["Gene Expression", "group_no"], developmental_trait_table["Life-history Traits", "group_no"],
                                                        developmental_trait_table["Morphology", "group_no"], developmental_trait_table["Physiological", "group_no"]), 
                                          "Effect Sizes" = c(developmental_trait_table["Behavioural", "K"], developmental_trait_table["Biochemical Assay", "K"], developmental_trait_table["Gene Expression", "K"], developmental_trait_table["Life-history Traits", "K"],
                                                             developmental_trait_table["Morphology", "K"], developmental_trait_table["Physiological", "K"]),
                                          "Estimate" = c(Developmental_Trait_Model_CVR$b[[1]], Developmental_Trait_Model_CVR$b[[2]], Developmental_Trait_Model_CVR$b[[3]], Developmental_Trait_Model_CVR$b[[4]],
                                                         Developmental_Trait_Model_CVR$b[[5]], Developmental_Trait_Model_CVR$b[[6]]),
                                          "CI Low" = c(Developmental_Trait_Model_CVR$ci.lb[1], Developmental_Trait_Model_CVR$ci.lb[2], Developmental_Trait_Model_CVR$ci.lb[3], Developmental_Trait_Model_CVR$ci.lb[4], 
                                                       Developmental_Trait_Model_CVR$ci.lb[5], Developmental_Trait_Model_CVR$ci.lb[6]), 
                                          "CI High" = c(Developmental_Trait_Model_CVR$ci.ub[1], Developmental_Trait_Model_CVR$ci.ub[2], Developmental_Trait_Model_CVR$ci.ub[3], Developmental_Trait_Model_CVR$ci.ub[4], 
                                                        Developmental_Trait_Model_CVR$ci.ub[5], Developmental_Trait_Model_CVR$ci.ub[6]), 
                                          "df" = c(Developmental_Trait_Model_CVR$ddf[[1]], Developmental_Trait_Model_CVR$ddf[[2]], Developmental_Trait_Model_CVR$ddf[[3]], Developmental_Trait_Model_CVR$ddf[[4]], 
                                                   Developmental_Trait_Model_CVR$ddf[[5]], Developmental_Trait_Model_CVR$ddf[[6]]), 
                                          "p-value" = c(Developmental_Trait_Model_CVR$pval[1], Developmental_Trait_Model_CVR$pval[2], Developmental_Trait_Model_CVR$pval[3], Developmental_Trait_Model_CVR$pval[4], 
                                                        Developmental_Trait_Model_CVR$pval[5], Developmental_Trait_Model_CVR$pval[6]))

Raw_Developmental_Exposure_CVR <- data.frame("Exposure Time" = c("Embryo", "Juvenile", "Larva", "Pupa"),
                                             "Studies" = c(developmental_exposure_study["Embryo", "Study"], developmental_exposure_study["Juvenile", "Study"], 
                                                           developmental_exposure_study["Larva", "Study"], developmental_exposure_study["Pupa", "Study"]), 
                                             "Species" = c(developmental_exposure_table["Embryo", "group_no"], developmental_exposure_table["Juvenile", "group_no"], 
                                                           developmental_exposure_table["Larva", "group_no"], developmental_exposure_table["Pupa", "group_no"]), 
                                             "Effect Sizes" = c(developmental_exposure_table["Embryo", "K"], developmental_exposure_table["Juvenile", "K"], 
                                                                developmental_exposure_table["Larva", "K"], developmental_exposure_table["Pupa", "K"]),
                                             "Estimate" = c(Developmental_Exposure_Model_CVR$b[[1]], Developmental_Exposure_Model_CVR$b[[2]],
                                                            Developmental_Exposure_Model_CVR$b[[3]], Developmental_Exposure_Model_CVR$b[[4]]),
                                             "CI Low" = c(Developmental_Exposure_Model_CVR$ci.lb[1], Developmental_Exposure_Model_CVR$ci.lb[2], 
                                                          Developmental_Exposure_Model_CVR$ci.lb[3], Developmental_Exposure_Model_CVR$ci.lb[4]), 
                                             "CI High" = c(Developmental_Exposure_Model_CVR$ci.ub[1], Developmental_Exposure_Model_CVR$ci.ub[2], 
                                                           Developmental_Exposure_Model_CVR$ci.ub[3], Developmental_Exposure_Model_CVR$ci.ub[4]), 
                                             "df" = c(Developmental_Exposure_Model_CVR$ddf[[1]], Developmental_Exposure_Model_CVR$ddf[[2]], 
                                                      Developmental_Exposure_Model_CVR$ddf[[3]], Developmental_Exposure_Model_CVR$ddf[[4]]), 
                                             "p-value" = c(Developmental_Exposure_Model_CVR$pval[1], Developmental_Exposure_Model_CVR$pval[2], 
                                                           Developmental_Exposure_Model_CVR$pval[3], Developmental_Exposure_Model_CVR$pval[4]))

Raw_Developmental_Class_CVR <- data.frame("Taxonomic Class" = c("Actinopteri", "Amphibia", "Arachnida", "Branchiopoda", "Insecta"),
                                          "Studies" = c(developmental_class_study["Actinopteri", "Study"], developmental_class_study["Amphibia", "Study"], developmental_class_study["Arachnida", "Study"], 
                                                        developmental_class_study["Branchiopoda", "Study"], developmental_class_study["Insecta", "Study"]), 
                                          "Species" = c(developmental_class_table["Actinopteri", "group_no"], developmental_class_table["Amphibia", "group_no"], developmental_class_table["Arachnida", "group_no"], 
                                                        developmental_class_table["Branchiopoda", "group_no"], developmental_class_table["Insecta", "group_no"]), 
                                          "Effect Sizes" = c(developmental_class_table["Actinopteri", "K"], developmental_class_table["Amphibia", "K"], developmental_class_table["Arachnida", "K"], 
                                                             developmental_class_table["Branchiopoda", "K"], developmental_class_table["Insecta", "K"]),
                                          "Estimate" = c(Developmental_Class_Model_CVR$b[[1]], Developmental_Class_Model_CVR$b[[2]], Developmental_Class_Model_CVR$b[[3]], 
                                                         Developmental_Class_Model_CVR$b[[4]], Developmental_Class_Model_CVR$b[[5]]),
                                          "CI Low" = c(Developmental_Class_Model_CVR$ci.lb[1], Developmental_Class_Model_CVR$ci.lb[2], Developmental_Class_Model_CVR$ci.lb[3], 
                                                       Developmental_Class_Model_CVR$ci.lb[4], Developmental_Class_Model_CVR$ci.lb[5]), 
                                          "CI High" = c(Developmental_Class_Model_CVR$ci.ub[1], Developmental_Class_Model_CVR$ci.ub[2], Developmental_Class_Model_CVR$ci.ub[3], 
                                                        Developmental_Class_Model_CVR$ci.ub[4], Developmental_Class_Model_CVR$ci.ub[5]), 
                                          "df" = c(Developmental_Class_Model_CVR$ddf[[1]], Developmental_Class_Model_CVR$ddf[[2]], Developmental_Class_Model_CVR$ddf[[3]], 
                                                   Developmental_Class_Model_CVR$ddf[[4]], Developmental_Class_Model_CVR$ddf[[5]]), 
                                          "p-value" = c(Developmental_Class_Model_CVR$pval[1], Developmental_Class_Model_CVR$pval[2], Developmental_Class_Model_CVR$pval[3], 
                                                        Developmental_Class_Model_CVR$pval[4], Developmental_Class_Model_CVR$pval[5]))

Raw_Developmental_Specific_Trait_CVR <- data.frame("Specific Phenotypic Traits" = c("Development Time", "Fecundity", "Food Consumption", 
                                                                                    "Head Width", "Length", "Locomotor Performance", 
                                                                                    "Longevity", "Mass", "Tail Length"),
                                                   "Studies" = c(developmental_specific_trait_study["Development Time", "Study"], developmental_specific_trait_study["Fecundity", "Study"], developmental_specific_trait_study["Food Consumption", "Study"], 
                                                                 developmental_specific_trait_study["Head Width", "Study"], developmental_specific_trait_study["Length", "Study"], developmental_specific_trait_study["Locomotor Performance", "Study"], 
                                                                 developmental_specific_trait_study["Longevity", "Study"], developmental_specific_trait_study["Mass", "Study"], developmental_specific_trait_study["Tail Length", "Study"]), 
                                                   "Species" = c(developmental_specific_trait_table["Development Time", "group_no"], developmental_specific_trait_table["Fecundity", "group_no"], developmental_specific_trait_table["Food Consumption", "group_no"], 
                                                                 developmental_specific_trait_table["Head Width", "group_no"], developmental_specific_trait_table["Length", "group_no"], developmental_specific_trait_table["Locomotor Performance", "group_no"], 
                                                                 developmental_specific_trait_table["Longevity", "group_no"], developmental_specific_trait_table["Mass", "group_no"], developmental_specific_trait_table["Tail Length", "group_no"]), 
                                                   "Effect Sizes" = c(developmental_specific_trait_table["Development Time", "K"], developmental_specific_trait_table["Fecundity", "K"], developmental_specific_trait_table["Food Consumption", "K"], 
                                                                      developmental_specific_trait_table["Head Width", "K"], developmental_specific_trait_table["Length", "K"], developmental_specific_trait_table["Locomotor Performance", "K"], 
                                                                      developmental_specific_trait_table["Longevity", "K"], developmental_specific_trait_table["Mass", "K"], developmental_specific_trait_table["Tail Length", "K"]),
                                                   "Estimate" = c(Developmental_Specific_Trait_Model_CVR$b[[1]], Developmental_Specific_Trait_Model_CVR$b[[2]], Developmental_Specific_Trait_Model_CVR$b[[3]], 
                                                                  Developmental_Specific_Trait_Model_CVR$b[[4]], Developmental_Specific_Trait_Model_CVR$b[[5]], Developmental_Specific_Trait_Model_CVR$b[[6]], 
                                                                  Developmental_Specific_Trait_Model_CVR$b[[7]], Developmental_Specific_Trait_Model_CVR$b[[8]], Developmental_Specific_Trait_Model_CVR$b[[9]]),
                                                   "CI Low" = c(Developmental_Specific_Trait_Model_CVR$ci.lb[1], Developmental_Specific_Trait_Model_CVR$ci.lb[2], Developmental_Specific_Trait_Model_CVR$ci.lb[3], 
                                                                Developmental_Specific_Trait_Model_CVR$ci.lb[4], Developmental_Specific_Trait_Model_CVR$ci.lb[5], Developmental_Specific_Trait_Model_CVR$ci.lb[6], 
                                                                Developmental_Specific_Trait_Model_CVR$ci.lb[7], Developmental_Specific_Trait_Model_CVR$ci.lb[8], Developmental_Specific_Trait_Model_CVR$ci.lb[9]), 
                                                   "CI High" = c(Developmental_Specific_Trait_Model_CVR$ci.ub[1], Developmental_Specific_Trait_Model_CVR$ci.ub[2], Developmental_Specific_Trait_Model_CVR$ci.ub[3], 
                                                                 Developmental_Specific_Trait_Model_CVR$ci.ub[4], Developmental_Specific_Trait_Model_CVR$ci.ub[5], Developmental_Specific_Trait_Model_CVR$ci.ub[6], 
                                                                 Developmental_Specific_Trait_Model_CVR$ci.ub[7], Developmental_Specific_Trait_Model_CVR$ci.ub[8], Developmental_Specific_Trait_Model_CVR$ci.ub[9]), 
                                                   "df" = c(Developmental_Specific_Trait_Model_CVR$ddf[[1]], Developmental_Specific_Trait_Model_CVR$ddf[[2]], Developmental_Specific_Trait_Model_CVR$ddf[[3]], 
                                                            Developmental_Specific_Trait_Model_CVR$ddf[[4]], Developmental_Specific_Trait_Model_CVR$ddf[[5]], Developmental_Specific_Trait_Model_CVR$ddf[[6]], 
                                                            Developmental_Specific_Trait_Model_CVR$ddf[[7]], Developmental_Specific_Trait_Model_CVR$ddf[[8]], Developmental_Specific_Trait_Model_CVR$ddf[[9]]), 
                                                   "p-value" = c(Developmental_Specific_Trait_Model_CVR$pval[1], Developmental_Specific_Trait_Model_CVR$pval[2], Developmental_Specific_Trait_Model_CVR$pval[3], 
                                                                 Developmental_Specific_Trait_Model_CVR$pval[4], Developmental_Specific_Trait_Model_CVR$pval[5], Developmental_Specific_Trait_Model_CVR$pval[6], 
                                                                 Developmental_Specific_Trait_Model_CVR$pval[7], Developmental_Specific_Trait_Model_CVR$pval[8], Developmental_Specific_Trait_Model_CVR$pval[9]))

write.csv(Raw_Overall_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Overall_CVR.csv", row.names = FALSE)
write.csv(Raw_Amplitude_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Amplitude_CVR.csv", row.names = FALSE)
write.csv(Raw_Fluctuation_Type_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Fluctuation_Type_CVR.csv", row.names = FALSE)
write.csv(Raw_Trait_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Trait_CVR.csv", row.names = FALSE)
write.csv(Raw_Class_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Class_CVR.csv", row.names = FALSE)
write.csv(Raw_Specific_Trait_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Specific_Trait_CVR.csv", row.names = FALSE)
write.csv(Raw_Population_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Population_CVR.csv", row.names = FALSE)
write.csv(Raw_Population_Amplitude_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Population_Amplitude_CVR.csv", row.names = FALSE)
write.csv(Raw_Population_Fluctuation_Type_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Population_Fluctuation_Type_CVR.csv", row.names = FALSE)
write.csv(Raw_Population_Class_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Population_Class_CVR.csv", row.names = FALSE)
write.csv(Raw_Individual_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Individual_CVR.csv", row.names = FALSE)
write.csv(Raw_Individual_Amplitude_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Individual_Amplitude_CVR.csv", row.names = FALSE)
write.csv(Raw_Individual_Fluctuation_Type_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Individual_Fluctuation_Type_CVR.csv", row.names = FALSE)
write.csv(Raw_Individual_Class_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Individual_Class_CVR.csv", row.names = FALSE)
write.csv(Raw_Aquatic_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Aquatic_CVR.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Amplitude_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Aquatic_Amplitude_CVR.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Fluctuation_Type_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Aquatic_Fluctuation_Type_CVR.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Trait_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Aquatic_Trait_CVR.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Plasticity_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Aquatic_Plasticity_CVR.csv", row.names = FALSE)
write.csv(Raw_Aquatic_Specific_Trait_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Aquatic_Specific_Trait_CVR.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Terrestrial_CVR.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Amplitude_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Terrestrial_Amplitude_CVR.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Fluctuation_Type_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Terrestrial_Fluctuation_Type_CVR.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Trait_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Terrestrial_Trait_CVR.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Plasticity_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Terrestrial_Plasticity_CVR.csv", row.names = FALSE)
write.csv(Raw_Terrestrial_Specific_Trait_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Terrestrial_Specific_Trait_CVR.csv", row.names = FALSE)
write.csv(Raw_Acclimation_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Acclimation_CVR.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Amplitude_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Acclimation_Amplitude_CVR.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Exposure_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Acclimation_Exposure_CVR.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Frequency_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Acclimation_Frequency_CVR.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Fluctuation_Type_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Acclimation_Fluctuation_Type_CVR.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Trait_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Acclimation_Trait_CVR.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Stage_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Acclimation_Stage_CVR.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Class_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Acclimation_Class_CVR.csv", row.names = FALSE)
write.csv(Raw_Acclimation_Specific_Trait_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Acclimation_Specific_Trait_CVR.csv", row.names = FALSE)
write.csv(Raw_Developmental_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Developmental_CVR.csv", row.names = FALSE)
write.csv(Raw_Developmental_Amplitude_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Developmental_Amplitude_CVR.csv", row.names = FALSE)
write.csv(Raw_Developmental_Fluctuation_Type_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Developmental_Fluctuation_Type_CVR.csv", row.names = FALSE)
write.csv(Raw_Developmental_Trait_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Developmental_Trait_CVR.csv", row.names = FALSE)
write.csv(Raw_Developmental_Exposure_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Developmental_Exposure_CVR.csv", row.names = FALSE)
write.csv(Raw_Developmental_Class_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Developmental_Class_CVR.csv", row.names = FALSE)
write.csv(Raw_Developmental_Specific_Trait_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Raw_Developmental_Specific_Trait_CVR.csv", row.names = FALSE)

# Heterogeneity Table

Heterogeneity_Overall_CVR <- data.frame("Models" = c("Overall", "Fluctuation Amplitude", "Fluctuation Type", 
                                                     "Phenotypic Trait Categories", "Specific Phenotypic Traits", "Taxonomic Class"), 
                                        "Shared Animal" = c(Overall_Model_CVR_i2[6, 1], Amplitude_Model_CVR_i2[6, 1], Fluctuation_Model_CVR_i2[6, 1], 
                                                            Trait_Model_CVR_i2[6, 1], Specific_Trait_Model_CVR_i2[6, 1], Class_Model_CVR_i2[6, 1]), 
                                        "Measurement" = c(Overall_Model_CVR_i2[7, 1], Amplitude_Model_CVR_i2[7, 1], Fluctuation_Model_CVR_i2[7, 1], 
                                                          Trait_Model_CVR_i2[7, 1], NA, Class_Model_CVR_i2[7, 1]),
                                        "Observational" = c(Overall_Model_CVR_i2[4, 1], Amplitude_Model_CVR_i2[4, 1], Fluctuation_Model_CVR_i2[4, 1], 
                                                            Trait_Model_CVR_i2[4, 1], Specific_Trait_Model_CVR_i2[4, 1], Class_Model_CVR_i2[4, 1]), 
                                        "Phylogenetic Relatedness" = c(Overall_Model_CVR_i2[2, 1], Amplitude_Model_CVR_i2[2, 1], Fluctuation_Model_CVR_i2[2, 1], 
                                                                       Trait_Model_CVR_i2[2, 1], Specific_Trait_Model_CVR_i2[2, 1], Class_Model_CVR_i2[2, 1]), 
                                        "Species" = c(Overall_Model_CVR_i2[5, 1], Amplitude_Model_CVR_i2[5, 1], Fluctuation_Model_CVR_i2[5, 1], 
                                                      Trait_Model_CVR_i2[5, 1], Specific_Trait_Model_CVR_i2[5, 1], Class_Model_CVR_i2[5, 1]),
                                        "Study" = c(Overall_Model_CVR_i2[3, 1], Amplitude_Model_CVR_i2[3, 1], Fluctuation_Model_CVR_i2[3, 1], 
                                                    Trait_Model_CVR_i2[3, 1], Specific_Trait_Model_CVR_i2[3, 1], Class_Model_CVR_i2[3, 1]), 
                                        "Total" = c(Overall_Model_CVR_i2[1, 1], Amplitude_Model_CVR_i2[1, 1], Fluctuation_Model_CVR_i2[1, 1], 
                                                    Trait_Model_CVR_i2[1, 1], Specific_Trait_Model_CVR_i2[1, 1], Class_Model_CVR_i2[1, 1]))

Heterogeneity_Population_CVR <- data.frame("Models" = c("Population-level Traits", "Fluctuation Amplitude", 
                                                        "Fluctuation Type", "Taxonomic Class"), 
                                           "Shared Animal" = c(Population_Model_CVR_i2[6, 1], Population_Amplitude_Model_CVR_i2[6, 1], 
                                                               Population_Fluctuation_Model_CVR_i2[6, 1], Population_Class_Model_CVR_i2[6, 1]), 
                                           "Measurement" = c(Population_Model_CVR_i2[7, 1], Population_Amplitude_Model_CVR_i2[7, 1], 
                                                             Population_Fluctuation_Model_CVR_i2[7, 1], Population_Class_Model_CVR_i2[7, 1]),
                                           "Observational" = c(Population_Model_CVR_i2[4, 1], Population_Amplitude_Model_CVR_i2[4, 1], 
                                                               Population_Fluctuation_Model_CVR_i2[4, 1], Population_Class_Model_CVR_i2[4, 1]), 
                                           "Phylogenetic Relatedness" = c(Population_Model_CVR_i2[2, 1], Population_Amplitude_Model_CVR_i2[2, 1], 
                                                                          Population_Fluctuation_Model_CVR_i2[2, 1], Population_Class_Model_CVR_i2[2, 1]), 
                                           "Species" = c(Population_Model_CVR_i2[5, 1], Population_Amplitude_Model_CVR_i2[5, 1], 
                                                         Population_Fluctuation_Model_CVR_i2[5, 1], Population_Class_Model_CVR_i2[5, 1]),
                                           "Study" = c(Population_Model_CVR_i2[3, 1], Population_Amplitude_Model_CVR_i2[3, 1], 
                                                       Population_Fluctuation_Model_CVR_i2[3, 1], Population_Class_Model_CVR_i2[3, 1]), 
                                           "Total" = c(Population_Model_CVR_i2[1, 1], Population_Amplitude_Model_CVR_i2[1, 1], 
                                                       Population_Fluctuation_Model_CVR_i2[1, 1], Population_Class_Model_CVR_i2[1, 1]))

Heterogeneity_Individual_CVR <- data.frame("Models" = c("Individual-level Traits", "Fluctuation Amplitude", 
                                                        "Fluctuation Type", "Taxonomic Class"), 
                                           "Shared Animal" = c(Individual_Model_CVR_i2[6, 1], Individual_Amplitude_Model_CVR_i2[6, 1], 
                                                               Individual_Fluctuation_Model_CVR_i2[6, 1], Individual_Class_Model_CVR_i2[6, 1]), 
                                           "Measurement" = c(Individual_Model_CVR_i2[7, 1], Individual_Amplitude_Model_CVR_i2[7, 1], 
                                                             Individual_Fluctuation_Model_CVR_i2[7, 1], Individual_Class_Model_CVR_i2[7, 1]),
                                           "Observational" = c(Individual_Model_CVR_i2[4, 1], Individual_Amplitude_Model_CVR_i2[4, 1], 
                                                               Individual_Fluctuation_Model_CVR_i2[4, 1], Individual_Class_Model_CVR_i2[4, 1]), 
                                           "Phylogenetic Relatedness" = c(Individual_Model_CVR_i2[2, 1], Individual_Amplitude_Model_CVR_i2[2, 1], 
                                                                          Individual_Fluctuation_Model_CVR_i2[2, 1], Individual_Class_Model_CVR_i2[2, 1]), 
                                           "Species" = c(Individual_Model_CVR_i2[5, 1], Individual_Amplitude_Model_CVR_i2[5, 1], 
                                                         Individual_Fluctuation_Model_CVR_i2[5, 1], Individual_Class_Model_CVR_i2[5, 1]),
                                           "Study" = c(Individual_Model_CVR_i2[3, 1], Individual_Amplitude_Model_CVR_i2[3, 1], 
                                                       Individual_Fluctuation_Model_CVR_i2[3, 1], Individual_Class_Model_CVR_i2[3, 1]), 
                                           "Total" = c(Individual_Model_CVR_i2[1, 1], Individual_Amplitude_Model_CVR_i2[1, 1], 
                                                       Individual_Fluctuation_Model_CVR_i2[1, 1], Individual_Class_Model_CVR_i2[1, 1]))

Heterogeneity_Aquatic_CVR <- data.frame("Models" = c("Aquatic Organisms", "Exposure Type", "Fluctuation Amplitude", 
                                                     "Fluctuation Type", "Phenotypic Trait Categories", "Specific Phenotypic Traits"), 
                                        "Shared Animal" = c(Aquatic_Model_CVR_i2[6, 1], Aquatic_Plasticity_Model_CVR_i2[6, 1], Aquatic_Amplitude_Model_CVR_i2[6, 1], 
                                                            Aquatic_Fluctuation_Model_CVR_i2[6, 1], Aquatic_Trait_Model_CVR_i2[6, 1], Aquatic_Specific_Trait_Model_CVR_i2[6, 1]), 
                                        "Measurement" = c(Aquatic_Model_CVR_i2[7, 1], Aquatic_Plasticity_Model_CVR_i2[7, 1], Aquatic_Amplitude_Model_CVR_i2[7, 1], 
                                                          Aquatic_Fluctuation_Model_CVR_i2[7, 1], Aquatic_Trait_Model_CVR_i2[7, 1], NA),
                                        "Observational" = c(Aquatic_Model_CVR_i2[4, 1], Aquatic_Plasticity_Model_CVR_i2[4, 1], Aquatic_Amplitude_Model_CVR_i2[4, 1], 
                                                            Aquatic_Fluctuation_Model_CVR_i2[4, 1], Aquatic_Trait_Model_CVR_i2[4, 1], Aquatic_Specific_Trait_Model_CVR_i2[4, 1]), 
                                        "Phylogenetic Relatedness" = c(Aquatic_Model_CVR_i2[2, 1], Aquatic_Plasticity_Model_CVR_i2[2, 1], Aquatic_Amplitude_Model_CVR_i2[2, 1], 
                                                                       Aquatic_Fluctuation_Model_CVR_i2[2, 1], Aquatic_Trait_Model_CVR_i2[2, 1], Aquatic_Specific_Trait_Model_CVR_i2[2, 1]), 
                                        "Species" = c(Aquatic_Model_CVR_i2[5, 1], Aquatic_Plasticity_Model_CVR_i2[5, 1], Aquatic_Amplitude_Model_CVR_i2[5, 1], 
                                                      Aquatic_Fluctuation_Model_CVR_i2[5, 1], Aquatic_Trait_Model_CVR_i2[5, 1], Aquatic_Specific_Trait_Model_CVR_i2[5, 1]),
                                        "Study" = c(Aquatic_Model_CVR_i2[3, 1], Aquatic_Plasticity_Model_CVR_i2[3, 1], Aquatic_Amplitude_Model_CVR_i2[3, 1], 
                                                    Aquatic_Fluctuation_Model_CVR_i2[3, 1], Aquatic_Trait_Model_CVR_i2[3, 1], Aquatic_Specific_Trait_Model_CVR_i2[3, 1]), 
                                        "Total" = c(Aquatic_Model_CVR_i2[1, 1], Aquatic_Plasticity_Model_CVR_i2[1, 1], Aquatic_Amplitude_Model_CVR_i2[1, 1], 
                                                    Aquatic_Fluctuation_Model_CVR_i2[1, 1], Aquatic_Trait_Model_CVR_i2[1, 1], Aquatic_Specific_Trait_Model_CVR_i2[1, 1]))

Heterogeneity_Terrestrial_CVR <- data.frame("Models" = c("Terrestrial Organisms", "Exposure Type", "Fluctuation Amplitude", 
                                                         "Fluctuation Type", "Phenotypic Trait Categories", "Specific Phenotypic Traits"), 
                                            "Shared Animal" = c(Terrestrial_Model_CVR_i2[6, 1], Terrestrial_Plasticity_Model_CVR_i2[6, 1], Terrestrial_Amplitude_Model_CVR_i2[6, 1], 
                                                                Terrestrial_Fluctuation_Model_CVR_i2[6, 1], Terrestrial_Trait_Model_CVR_i2[6, 1], Terrestrial_Specific_Trait_Model_CVR_i2[6, 1]),
                                            "Measurement" = c(Terrestrial_Model_CVR_i2[7, 1], Terrestrial_Plasticity_Model_CVR_i2[7, 1], Terrestrial_Amplitude_Model_CVR_i2[7, 1], 
                                                              Terrestrial_Fluctuation_Model_CVR_i2[7, 1], Terrestrial_Trait_Model_CVR_i2[7, 1], NA),
                                            "Observational" = c(Terrestrial_Model_CVR_i2[4, 1], Terrestrial_Plasticity_Model_CVR_i2[4, 1], Terrestrial_Amplitude_Model_CVR_i2[4, 1], 
                                                                Terrestrial_Fluctuation_Model_CVR_i2[4, 1], Terrestrial_Trait_Model_CVR_i2[4, 1], Terrestrial_Specific_Trait_Model_CVR_i2[4, 1]), 
                                            "Phylogenetic Relatedness" = c(Terrestrial_Model_CVR_i2[2, 1], Terrestrial_Plasticity_Model_CVR_i2[2, 1], Terrestrial_Amplitude_Model_CVR_i2[2, 1], 
                                                                           Terrestrial_Fluctuation_Model_CVR_i2[2, 1], Terrestrial_Trait_Model_CVR_i2[2, 1], Terrestrial_Specific_Trait_Model_CVR_i2[2, 1]), 
                                            "Species" = c(Terrestrial_Model_CVR_i2[5, 1], Terrestrial_Plasticity_Model_CVR_i2[5, 1], Terrestrial_Amplitude_Model_CVR_i2[5, 1], 
                                                          Terrestrial_Fluctuation_Model_CVR_i2[5, 1], Terrestrial_Trait_Model_CVR_i2[5, 1], Terrestrial_Specific_Trait_Model_CVR_i2[5, 1]),
                                            "Study" = c(Terrestrial_Model_CVR_i2[3, 1], Terrestrial_Plasticity_Model_CVR_i2[3, 1], Terrestrial_Amplitude_Model_CVR_i2[3, 1], 
                                                        Terrestrial_Fluctuation_Model_CVR_i2[3, 1], Terrestrial_Trait_Model_CVR_i2[3, 1], Terrestrial_Specific_Trait_Model_CVR_i2[3, 1]), 
                                            "Total" = c(Terrestrial_Model_CVR_i2[1, 1], Terrestrial_Plasticity_Model_CVR_i2[1, 1], Terrestrial_Amplitude_Model_CVR_i2[1, 1], 
                                                        Terrestrial_Fluctuation_Model_CVR_i2[1, 1], Terrestrial_Trait_Model_CVR_i2[1, 1], Terrestrial_Specific_Trait_Model_CVR_i2[1, 1]))

Heterogeneity_Acclimation_CVR <- data.frame("Models" = c("Acclimation", "Exposure Time", "Fluctuation Amplitude", "Fluctuation Type", "Life-history Stage", 
                                                         "Number of Fluctuations", "Phenotypic Trait Categories", "Specific Phenotypic Traits", "Taxonomic Class"), 
                                            "Shared Animal" = c(Acclimation_Model_CVR_i2[6, 1], Acclimation_Exposure_Model_CVR_i2[6, 1], Acclimation_Amplitude_Model_CVR_i2[6, 1], Acclimation_Fluctuation_Model_CVR_i2[6, 1], Acclimation_Stage_Model_CVR_i2[6, 1], 
                                                                Acclimation_Frequency_Model_CVR_i2[6, 1], Acclimation_Trait_Model_CVR_i2[6, 1], Acclimation_Specific_Trait_Model_CVR_i2[6, 1], Acclimation_Class_Model_CVR_i2[6, 1]), 
                                            "Measurement" = c(Acclimation_Model_CVR_i2[7, 1], Acclimation_Exposure_Model_CVR_i2[7, 1], Acclimation_Amplitude_Model_CVR_i2[7, 1], Acclimation_Fluctuation_Model_CVR_i2[7, 1], Acclimation_Stage_Model_CVR_i2[7, 1], 
                                                              Acclimation_Frequency_Model_CVR_i2[7, 1], Acclimation_Trait_Model_CVR_i2[7, 1], NA, Acclimation_Class_Model_CVR_i2[7, 1]), 
                                            "Observational" = c(Acclimation_Model_CVR_i2[4, 1], Acclimation_Exposure_Model_CVR_i2[4, 1], Acclimation_Amplitude_Model_CVR_i2[4, 1], Acclimation_Fluctuation_Model_CVR_i2[4, 1], Acclimation_Stage_Model_CVR_i2[4, 1], 
                                                                Acclimation_Frequency_Model_CVR_i2[4, 1], Acclimation_Trait_Model_CVR_i2[4, 1], Acclimation_Specific_Trait_Model_CVR_i2[4, 1], Acclimation_Class_Model_CVR_i2[4, 1]), 
                                            "Phylogenetic Relatedness" = c(Acclimation_Model_CVR_i2[2, 1], Acclimation_Exposure_Model_CVR_i2[2, 1], Acclimation_Amplitude_Model_CVR_i2[2, 1], Acclimation_Fluctuation_Model_CVR_i2[2, 1], Acclimation_Stage_Model_CVR_i2[2, 1], 
                                                                           Acclimation_Frequency_Model_CVR_i2[2, 1], Acclimation_Trait_Model_CVR_i2[2, 1], Acclimation_Specific_Trait_Model_CVR_i2[2, 1], Acclimation_Class_Model_CVR_i2[2, 1]), 
                                            "Species" = c(Acclimation_Model_CVR_i2[5, 1], Acclimation_Exposure_Model_CVR_i2[5, 1], Acclimation_Amplitude_Model_CVR_i2[5, 1], Acclimation_Fluctuation_Model_CVR_i2[5, 1], Acclimation_Stage_Model_CVR_i2[5, 1], 
                                                          Acclimation_Frequency_Model_CVR_i2[5, 1], Acclimation_Trait_Model_CVR_i2[5, 1], Acclimation_Specific_Trait_Model_CVR_i2[5, 1], Acclimation_Class_Model_CVR_i2[5, 1]),
                                            "Study" = c(Acclimation_Model_CVR_i2[3, 1], Acclimation_Exposure_Model_CVR_i2[3, 1], Acclimation_Amplitude_Model_CVR_i2[3, 1], Acclimation_Fluctuation_Model_CVR_i2[3, 1], Acclimation_Stage_Model_CVR_i2[3, 1], 
                                                        Acclimation_Frequency_Model_CVR_i2[3, 1], Acclimation_Trait_Model_CVR_i2[3, 1], Acclimation_Specific_Trait_Model_CVR_i2[3, 1], Acclimation_Class_Model_CVR_i2[3, 1]), 
                                            "Total" = c(Acclimation_Model_CVR_i2[1, 1], Acclimation_Exposure_Model_CVR_i2[1, 1], Acclimation_Amplitude_Model_CVR_i2[1, 1], Acclimation_Fluctuation_Model_CVR_i2[1, 1], Acclimation_Stage_Model_CVR_i2[1, 1], 
                                                        Acclimation_Frequency_Model_CVR_i2[1, 1], Acclimation_Trait_Model_CVR_i2[1, 1], Acclimation_Specific_Trait_Model_CVR_i2[1, 1], Acclimation_Class_Model_CVR_i2[1, 1]))

Heterogeneity_Developmental_CVR <- data.frame("Models" = c("Developmental", "Exposure Time", "Fluctuation Amplitude", "Fluctuation Type", 
                                                           "Phenotypic Trait Categories", "Specific Phenotypic Traits", "Taxonomic Class"), 
                                              "Shared Animal" = c(Developmental_Model_CVR_i2[6, 1], Developmental_Exposure_Model_CVR_i2[6, 1], Developmental_Amplitude_Model_CVR_i2[6, 1], Developmental_Fluctuation_Model_CVR_i2[6, 1], 
                                                                  Developmental_Trait_Model_CVR_i2[6, 1], Developmental_Specific_Trait_Model_CVR_i2[6, 1], Developmental_Class_Model_CVR_i2[6, 1]), 
                                              "Measurement" = c(Developmental_Model_CVR_i2[7, 1], Developmental_Exposure_Model_CVR_i2[7, 1], Developmental_Amplitude_Model_CVR_i2[7, 1], Developmental_Fluctuation_Model_CVR_i2[7, 1], 
                                                                Developmental_Trait_Model_CVR_i2[7, 1], NA, Developmental_Class_Model_CVR_i2[7, 1]),
                                              "Observational" = c(Developmental_Model_CVR_i2[4, 1], Developmental_Exposure_Model_CVR_i2[4, 1], Developmental_Amplitude_Model_CVR_i2[4, 1], Developmental_Fluctuation_Model_CVR_i2[4, 1], 
                                                                  Developmental_Trait_Model_CVR_i2[4, 1], Developmental_Specific_Trait_Model_CVR_i2[4, 1], Developmental_Class_Model_CVR_i2[4, 1]), 
                                              "Phylogenetic Relatedness" = c(Developmental_Model_CVR_i2[2, 1], Developmental_Exposure_Model_CVR_i2[2, 1], Developmental_Amplitude_Model_CVR_i2[2, 1], Developmental_Fluctuation_Model_CVR_i2[2, 1], 
                                                                             Developmental_Trait_Model_CVR_i2[2, 1], Developmental_Specific_Trait_Model_CVR_i2[2, 1], Developmental_Class_Model_CVR_i2[2, 1]), 
                                              "Species" = c(Developmental_Model_CVR_i2[5, 1], Developmental_Exposure_Model_CVR_i2[5, 1], Developmental_Amplitude_Model_CVR_i2[5, 1], Developmental_Fluctuation_Model_CVR_i2[5, 1], 
                                                            Developmental_Trait_Model_CVR_i2[5, 1], Developmental_Specific_Trait_Model_CVR_i2[5, 1], Developmental_Class_Model_CVR_i2[5, 1]),
                                              "Study" = c(Developmental_Model_CVR_i2[3, 1], Developmental_Exposure_Model_CVR_i2[3, 1], Developmental_Amplitude_Model_CVR_i2[3, 1], Developmental_Fluctuation_Model_CVR_i2[3, 1], 
                                                          Developmental_Trait_Model_CVR_i2[3, 1], Developmental_Specific_Trait_Model_CVR_i2[3, 1], Developmental_Class_Model_CVR_i2[3, 1]), 
                                              "Total" = c(Developmental_Model_CVR_i2[1, 1], Developmental_Exposure_Model_CVR_i2[1, 1], Developmental_Amplitude_Model_CVR_i2[1, 1], Developmental_Fluctuation_Model_CVR_i2[1, 1], 
                                                          Developmental_Trait_Model_CVR_i2[1, 1], Developmental_Specific_Trait_Model_CVR_i2[1, 1], Developmental_Class_Model_CVR_i2[1, 1]))

write.csv(Heterogeneity_Overall_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Heterogeneity_Overall_CVR.csv", row.names = FALSE)
write.csv(Heterogeneity_Population_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Heterogeneity_Population_CVR.csv", row.names = FALSE)
write.csv(Heterogeneity_Individual_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Heterogeneity_Individual_CVR.csv", row.names = FALSE)
write.csv(Heterogeneity_Aquatic_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Heterogeneity_Aquatic_CVR.csv", row.names = FALSE)
write.csv(Heterogeneity_Terrestrial_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Heterogeneity_Terrestrial_CVR.csv", row.names = FALSE)
write.csv(Heterogeneity_Acclimation_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Heterogeneity_Acclimation_CVR.csv", row.names = FALSE)
write.csv(Heterogeneity_Developmental_CVR, file = "./3.Data_Analysis/2.Outputs/Supplementary_Material/Heterogeneity_Developmental_CVR.csv", row.names = FALSE)
