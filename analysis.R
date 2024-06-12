#library

library(ggplot2)
library(tidyverse)
library(lme4)
library(lmerTest)
library(adegenet)
library(hierfstat)
library(interactions)
library(dplyr)
library(emmeans)


getwd()
setwd("C:/Uni/Master/Macroecology/project")

source('../genetic_biodiv_practical-main/other_functions.R')

#read and explore data

data <- read.csv("macrocourse_EU_data_2024.csv")
taxa <- read.csv2("taxonomy.csv")

data <- merge(data, taxa)

data$bats <- ifelse(data$order == "chiroptera", 1, 0)
data$bats <- as.factor(data$bats)

head(data)
str(data)

levels(as.factor(data$species))

data$species <- as.factor(data$species)
data$genus <- as.factor(data$genus)

hist(data$altitude_10)
hist(data$annual_precip_10)
hist(data$roughness_10)

bats <- data[data$bats== 1,]
non_bats <- data[data$bats == 0,]

#transform variables

data$log_alt_5 <- log(data$altitude_5+0.0000000001)
data$log_alt_10 <- log(data$altitude_10)
data$log_alt_50 <- log(data$altitude_50)
data$log_alt_100 <- log(data$altitude_100)

data$log_roughness_10 <- log(data$roughness_10)

hist(data$log_alt_100)


#modelling altitude effect

gd_environment_5 <- lmer(gene_diversity ~ log_alt_5 + 
                            (1|genus/species), data = data)
summary(gd_environment_5)

gd_environment_10 <- lmer(gene_diversity ~  log_alt_10 + 
                          (1|genus/species), data = data)
summary(gd_environment_10)

ar_environment_10 <- lmer(allelic_richness ~ log_alt_10 + (1|genus/species), data = data)
summary(ar_environment_10) # same trend


res <- residuals(gd_environment_10)
plot(res)
qqnorm(res)
qqline(res)

gd_environment_50 <- lmer(gene_diversity ~ log_alt_50 +
                         (1|genus/species), data = data)
summary(gd_environment_50)



gd_environment_100 <- lmer(gene_diversity ~ log_alt_100 + 
                            (1|genus/species), data = data)
summary(gd_environment_100)

summary(gd_environment_5)
summary(gd_environment_10)
summary(gd_environment_50)
summary(gd_environment_100)



ggplot(data = data, aes(x = log_alt_10, y = gene_diversity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()

interact_plot(gd_environment_50, pred = PET_50, modx = log_alt_50)+
  theme_bw()


model_bats <- lmer(gene_diversity ~ log_alt_10 * bats + (1|genus/species), data = data)
summary(model_bats)
anova(model_bats)

model_order <- lmer(gene_diversity ~ log_alt_10 * order + (1|genus/species), data = data)
summary(model_order)
anova(model_order)

model_order2 <- lmer(gene_diversity ~ log_alt_10 + (log_alt_10|order/species), data = data)
summary(model_order2)
coef(model_order2)

coef_data <- coef(model_order2)$order

coef_data <- as.data.frame(coef_data)

coef_data$factor_level <- rownames(coef_data)


emtrends(model_order, pairwise ~ order,var = "log_alt_10")

# modelling impact of altitude on allelic richness
model_ar_order <- lmer(allelic_richness ~ log_alt_10 * order + (1|genus/species), data = data)
summary(model_ar_order)
anova(model_ar_order) #similiar trend to gene diversity but no order significance (altitude + interaction significant)



ggplot(data = data, aes(x = log_alt_10, y = gene_diversity, color = order)) +
  geom_point() +
  geom_smooth(method = "lm", aes(col = order), se = FALSE)


#modelling impact of hfi

model_hfi <- lmer(gene_diversity ~ hfi_10 * order + (1|genus/species), data = data)
summary(model_hfi)
anova(model_hfi)

model_hfi_r <- lmer(gene_diversity ~ hfi_10 + (1|genus/species), data = rodents)
summary(model_hfi_r)

model_hfi_b <- lmer(gene_diversity ~ hfi_10 + (1|genus/species), data = bats)
summary(model_hfi_b)

emtrends(model_hfi, pairwise ~ order,var = "hfi_10")

# modelling allelic richness

model_ar_hfi <- lmer(allelic_richness ~ hfi_10 * order + (1|genus/species), data = data)
summary(model_ar_hfi)
anova(model_ar_hfi) # shows no trend (not like gene diversity) 


ggplot(data = data, aes(x = hfi_10, y = gene_diversity, color = order)) +
  geom_point() +
  geom_smooth(method = "lm", aes(col = order), se = FALSE)



#plots for rodents and bats (hfi)

rodents <- data[data$order == "rodentia",]

kein_lemming <- rodents[rodents$species!= "Lemmus_lemmus",]

ggplot(data = rodents, aes(x = log_alt_10, y = gene_diversity, color = species)) +
  geom_point() +
  geom_smooth(method = "lm", aes(col = NULL), se = FALSE)

ggplot(data = bats, aes(x = log_alt_10, y = gene_diversity, color = species)) +
  geom_point() +
  geom_smooth(method = "lm", aes(col = NULL), se = FALSE)

ggplot(data = kein_lemming, aes(x = hfi_10, y = gene_diversity, color = species)) +
  geom_point() +
  geom_smooth(method = "lm", aes(col = NULL), se = FALSE)


#models altitude + hfi

model_alt_hfi <- lmer(gene_diversity ~ log_alt_10 + hfi_10 + (1|genus/species), data = data)
summary(model_alt_hfi)


model_ar_both <- lmer(allelic_richness ~ log_alt_10 + hfi_10 + (1|genus/species), data = data)
summary(model_ar_both) # similar trend to gene diversity but hfi slighty above significance

#boxplot

ggplot(data= data, aes(x= order, y= gene_diversity, color = order)) +
  geom_boxplot()

