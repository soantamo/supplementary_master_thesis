###
# code for plots in method figure 

library(tidyverse)
library(here)
library(readr)
library(viridis)
library(readxl)
library(forcats)

# read in total GAM predictions
model_predictions <- readRDS("total_models/df_pred_all.rds")

model_predictions$species <- as.factor(model_predictions$species)
levels(model_predictions$species)

# plot predictions: Fig. 1a
method_df <- model_predictions |> 
  filter(species %in% c("Alburnus_alburnus", "Lota_lota", "Cottus_gobio_Aare_littoral"))

method_df_rescaled <- method_df |> 
  group_by(species) |> 
  mutate(fit_rescaled = (fit - min(fit)) / (max(fit) - min(fit)))

mycolors <-  c("Alburnus_alburnus" = "#FE6949", "Lota_lota" = "#5ADFFE", "Cottus_gobio_Aare_littoral" = "#FAB0B6")

method_fig_pred <- method_df_rescaled |> 
  ggplot(aes(temp, fit_rescaled)) +
  geom_line(aes(color = factor(species)), linewidth = 1.3) +
  theme_classic(base_size = 25) +
  scale_color_manual(values = mycolors, guide = NULL) +
  xlab("") +
  ylab("") +
  xlim(5.2, 20)
method_fig_pred

# save as TIFF
# tiff(paste("total_models/plots/2_method_figure_predictions.tiff", sep = ""), units="in", width = 12, height=8, res=300)
# plot(method_fig_pred)
# 
# # Closing the graphical device
# dev.off()


# plot derivatives of the species: Fig. 1b
# derivatives data

all_models_derivatives <- readRDS("total_models/df_deriv_all.rds")

all_deriv <- as_tibble(all_models_derivatives)

# select the species we need
method_deriv <- all_deriv |> 
  filter(species %in% c("Alburnus_alburnus", "Lota_lota", "Cottus_gobio_Aare_littoral"))

# plot derivatives of the three species -> Fig. 1b
method_fig_deriv <- method_deriv |> 
  ggplot(aes(temp, derivative)) +
  geom_line(aes(color = factor(species)), linewidth = 1.3) +
  theme_classic(base_size = 25) +
  scale_color_manual(values = mycolors, guide = NULL) +
  xlab("") +
  ylab("") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlim(5.2, 20)

# save as TIFF
# tiff(paste("total_models/plots/2_method_figure_derivatives.tiff", sep = ""), units="in", width = 12, height=8, res=300)
# 
# plot(method_fig_deriv)
# 
# # Closing the graphical device
# dev.off()


# Examples of lake communities

# Divergence: Fig. 1c
lake1_deriv <- all_deriv |> 
  filter(species %in% c("Alburnus_alburnus", "Coregonus_sp_felchen", "Salmo_trutta"))

mycolors2 <-  c("Alburnus_alburnus" = "black", "Coregonus_sp_felchen" = "black", "Salmo_trutta" = "black")

# plot Fig. 1c
lake1_fig_deriv <- lake1_deriv |> 
  ggplot(aes(temp, derivative)) +
  geom_line(aes(color = factor(species)), linewidth = 1.3) +
  theme_classic(base_size = 25) +
  scale_color_manual(values = mycolors2, guide = NULL) +
  xlab("") +
  ylab("") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlim(4.8, 20)

lake1_fig_deriv
# save the plot
# tiff(paste("total_models/plots/lake1_figure_derivatives.tiff", sep = ""), units="in", width = 12, height=8, res=300)
# 
# plot(lake1_fig_deriv)
# 
# # Closing the graphical device
# dev.off()

# calculate response diversity metrics for the three species
source(here("functions.R"))

lake1_deriv$fLake <- as.character(lake1_deriv$fLake)

# select the columns we need
lake1_data <- lake1_deriv |>
  select(temp, derivative, species, fLake)

# wider format to calculate 
df_resp_div_lake1 <- lake1_data |>
  pivot_wider(
    names_from = species,
    values_from = derivative) |> 
  drop_na()  |> 
  select(-fLake) #we do not need lake information here

# calculate
# dissimilarity
df_resp_div_lake1$rdiv <- apply(df_resp_div_lake1[,-1, drop = FALSE], 1, resp_div, sign_sens = F)
# divergence
df_resp_div_lake1$sign <- apply(df_resp_div_lake1[,-1, drop = FALSE], 1, resp_div, sign_sens = T)
# median dissimilarity
df_resp_div_lake1$Med <- median(df_resp_div_lake1$rdiv)

str(df_resp_div)

lake1_fig_div <- df_resp_div_lake1 |> 
  ggplot(aes(temp, sign)) +
  geom_line( linewidth = 1) +
  theme_classic(base_size = 25) +
  xlab("") +
  ylab("") +
  ylim(0,1)

lake1_fig_div

# save as tiff
# tiff(paste("total_models/plots/lake1_figure_div.tiff", sep = ""), units="in", width = 12, height=8, res=300)
# 
# plot(lake1_fig_div)
# 
# dev.off()

# Lake 2
# Dissimilarity: Fig. 1d
# Derivatives in this lake community
method_deriv <- all_deriv |> 
  filter(species %in% c("Alburnus_alburnus", "Lota_lota", "Cottus_gobio_Aare_littoral"))

mycolors2 <-  c("Alburnus_alburnus" = "black", "Lota_lota" = "black", "Cottus_gobio_Aare_littoral" = "black")

method_fig_deriv_black_2 <- method_deriv |> 
  ggplot(aes(temp, derivative)) +
  geom_line(aes(color = factor(species)), linewidth = 1.3) +
  theme_classic(base_size = 25) +
  scale_color_manual(values = mycolors2, guide = NULL) +
  xlab("") +
  ylab("") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlim(4.8, 20)

method_fig_deriv_black_2

# tiff(paste("total_models/plots/black_2_method_figure_derivatives.tiff", sep = ""), units="in", width = 12, height=8, res=300)
# # plot(ggarrange(depth1, depth2, ncol = 2))
# # plot science discussion
# plot(method_fig_deriv_black_2)
# 
# # Closing the graphical device
# dev.off()


# calculate response diversity metrics with the other lake community
# Dissimilarity metric

source(here("functions.R"))

method_deriv$fLake <- as.character(method_deriv$fLake)

data <- method_deriv |>
  select(temp, derivative, species, fLake)

df_resp_div <- data |>
  pivot_wider(
    names_from = species,
    values_from = derivative) |> 
  drop_na() |> 
  select(-fLake)

df_resp_div$rdiv <- apply(df_resp_div[,-1, drop = FALSE], 1, resp_div, sign_sens = F)
df_resp_div$sign <- apply(df_resp_div[,-1, drop = FALSE], 1, resp_div, sign_sens = T)
df_resp_div$Med <- median(df_resp_div$rdiv)

str(df_resp_div)

method_fig_diss <- df_resp_div |> 
  ggplot(aes(temp, rdiv)) +
  geom_line( linewidth = 1) +
  theme_classic(base_size = 25) +
  xlab("") +
  ylab("") +
  ylim(1,3) +
  xlim(5.2, 20)

method_fig_diss

# tiff(paste("total_models/plots/2_method_figure_diss.tiff", sep = ""), units="in", width = 12, height=8, res=300)
# 
# plot(method_fig_diss)
# 
# dev.off()


