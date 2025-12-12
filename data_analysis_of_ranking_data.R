################################################################################
# Ranking data analysis of on-farm trials  
# 
# Authority: Wageningen university
# Author: Joost van Heerwaarden and Hugo Dorado
# Institution/Project:Wageningen university
# Date:24/10/2025
# Description: 
#   This script estimates Thurstonian model parameters for maize yield ranking 
#   data collected from Rwanda (1000 Farms Project). It includes data preparation, 
#   model fitting, and computation of means and standard errors.
# Dependencies: tidyverse, PlackettLuce, gosset, qvcalc
################################################################################

#----------------------------- Load Libraries ----------------------------------
library(tidyverse)     # Data manipulation and visualization
library(PlackettLuce)  # Ranking model functions
library(gosset)        # Ranking data structure
library(qvcalc)        # Quasi-variance and uncertainty estimation
library(asreml)        # Fitting mixed-effect models
library(here)
#----------------------------- Set Working Directory ---------------------------
# Ensure the working directory points to the folder containing required scripts
setwd(here())

#----------------------------- Load Functions ----------------------------------
# Source custom functions used for Thurstonian model fitting
source('master_funs.R')

#----------------------------- Load Dataset ------------------------------------
# Read maize yield ranking data for Rwanda
dataset <- read.csv('Data/yield_maize_rwanda.csv')

#--------------------- Thurstonian Model Estimation ----------------------------

# Example of only one zone

# Filter dataset to include only Zone 2
dataset_zone2 <- dataset %>% 
  filter(zone == 2)

# Convert numeric yield rankings to ranking data format
# 'id' identifies individual surveys/farms, 
# 'items' are genotypes, and 'input' is the ranking variable
yield_rankings <- rank_numeric(
  data  = dataset_zone2,
  id    = 'registration_surveyid',
  items = 'genotype',
  input = 'yield_ranked'
)

# Inspect formatted rankings
yield_rankings

# Fit the Thurstonian model
# 'ref' indicates the reference genotype against which others are compared
thurtonian_model <- thurstonianMod(
  rnks = yield_rankings,
  ref  = 'wh 507'
)

# Compute Thurstonian means, standard errors, and quasi-standard errors
th_means_se(thurtonian_model)

#---------------------- GxE Interaction Analysis -------------------------------

# Environment is defined as a combination between zone and year

dataset$environment <- paste("Z",dataset$zone,"Y",dataset$year,sep='')

# This function fits the Thurstonian model for each environment
# and returns a table of genotypic means and standard errors.
# 'ref' specifies the reference genotype used for identification.

gmeans_table <- genotypic_means_table(data=dataset,env = 'environment',
                      id    = 'registration_surveyid',items = 'genotype',
                      input = 'yield_ranked',ref   ='wh 507')

# Optional: Save and reload the results

# write.csv(gmeans_table,'Results/gmeans_table.csv',row.names=F)
# gmeans_table <- read.csv('Results/gmeans_table.csv')

# The FA(2) model allows for modeling the genotype × environment
# interaction structure through two latent factors.

gmeans_table$wt <- 1/gmeans_table$quasiVar # define weights in the model

gmeans_table$environment  <- factor(gmeans_table$environment )

gmeans_table$genotype  <- factor(gmeans_table$genotype )

fa2_mod <- asreml( fixed = estimate~environment,random =~ fa(environment, 2):genotype, 
                  family = asreml::asr_gaussian(dispersion = 1), weights = wt,
                  na.action = list(x = "exclude", y = "include"), 
                  trace = 0, maxiter = 2000,data = gmeans_table)


summary(fa2_mod)

