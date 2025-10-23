# Longitudinal development of ITS - Tau values ---------------------------------
library(readr)
library(readxl)
library(tidyverse)
library(dplyr)
library(sjPlot)
library(multcomp)
library(AICcmodavg)
library(flextable)
library(rempsyc)
library(lme4)
library(lmerTest)

set.seed(42) #Set seed for reproducibility 

table_glht <- function(x) { # Function to extract results from glht objects
  pq <- summary(x)$test
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  error <- attr(pq$pvalues, "error")
  pname <- switch(x$alternativ, less = paste("Pr(<", ifelse(x$df ==0, "z", "t"), ")", sep = ""), 
                  greater = paste("Pr(>", ifelse(x$df == 0, "z", "t"), ")", sep = ""), two.sided = paste("Pr(>|",ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
  colnames(mtests) <- c("Estimate", "Std. Error", ifelse(x$df ==0, "z value", "t value"), pname)
  return(mtests)
  
}

# Where the data is
tablepath = "~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/"

## Covariates preparation -------------------------------
#### This part of the code should be removed and just load the covariates_its.csv file provided
ses <- readxl::read_excel("~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/bexat_data_quest.xlsx", sheet = "ses")
ses <- ses%>%
  mutate(MeanEduM = mean(ZSES_EdMadre_6m, na.rm = T), 
         MeanEduF = mean(ZSES_EdPadre_6m, na.rm = T),
         MeanOccM = mean(ZSES_LabMadre_6m, na.rm = T),
         MeanOccF = mean(ZSES_LabPadre_6m, na.rm = T),
         MeanI2N  = mean(IncNeed_6m, na.rm = T))%>%
  mutate(MEdu = if_else(!is.na(ZSES_EdMadre_6m), ZSES_EdMadre_6m, MeanEduM), 
         FEdu = if_else(!is.na(ZSES_EdPadre_6m), ZSES_EdPadre_6m, MeanEduF),
         MOcc = if_else(!is.na(ZSES_LabMadre_6m), ZSES_LabMadre_6m, if_else(!is.na(ZSES_LabMadre_16m), ZSES_LabMadre_16m,  MeanOccM)), 
         FOcc = if_else(!is.na(ZSES_LabPadre_6m), ZSES_LabPadre_6m, if_else(!is.na(ZSES_LabPadre_16m), ZSES_LabPadre_16m,  MeanOccF)),
         I2N  = if_else(!is.na(IncNeed_6m), IncNeed_6m, if_else(!is.na(IncNeed_16m), IncNeed_16m, MeanI2N)),
         by = suj)%>%dplyr::select(suj, MEdu, FEdu, MOcc, FOcc, I2N)

basicdemos <- read_delim("~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Results/Chapter_3/Figures_and_Tables/basicdemos.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)%>%filter(sesage==6)%>%dplyr::select(-1, -sesage)%>%
  mutate(meanbw = mean(bw, na.rm = T), 
         meangw = mean(gw, na.rm = T))%>%dplyr::select(-sex)

basicdemos$bw[is.na(basicdemos$bw)] = mean(basicdemos$meanbw)
basicdemos$gw[is.na(basicdemos$gw)] = mean(basicdemos$meangw)


bdi <- readxl::read_excel("~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/bexat_data_quest.xlsx", sheet = "bdi")
bdi <- bdi%>%
  mutate(Mbdi = mean(BDI_6m, na.rm = T))%>%
  mutate(BDI  = if_else(!is.na(BDI_6m), BDI_6m, if_else(!is.na(BDI_16m), BDI_16m, Mbdi)))%>%dplyr::select(suj, BDI)

covs = merge(ses, basicdemos, by = "suj")
covs = merge(covs, bdi, by = "suj")%>%dplyr::select(-value, -meanbw, -meangw)

colnames(covs) = c("suj", "MEdu", "FEdu", "MOcc", "FOcc", "I2N", "bw", "gw", "BDI")

write_csv(covs, file = "~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/covariates_its.csv")

covs = covs%>%mutate(bw = scale(bw, center = T, scale = T), 
                     gw = scale(gw, center = T, scale = T), 
                     BDI = scale (BDI, center = T, scale = T))

 # Cohorts 
cohorts = c("exploratory", "validation")

## Run the analysis ----------------------------------------------------------
for (cidx in cohorts){
  
  # Load data 
  if (cidx == "validation"){
    taudata = read_csv("Desktop/Tau_data/tau_development_validation_v4.csv")
  } else {
    taudata = read_csv("Desktop/Tau_data/tau_development_exploratory_v4.csv")
  }
  
  dataset <- taudata%>%filter(pernan < .25, event != 9999)%>% # Filter only the good epochs and remove bad events
    mutate(trials = mean(n())/length(unique(channel)), .by = c(subject, ses_age))%>%
    group_by(subject, Gender, ses_age, event_n, channel, Area)%>%
    summarise(tau    = mean(tauimp, na.rm = T), 
              trials = mean(trials, na.rm= T),
              pernan = mean(pernan, na.rm = T))
  
  
  dataset             <- dplyr::select(dataset, subject, Gender, ses_age, channel, Area, tau, event_n, trials, pernan)%>%
    group_by(subject, Gender, ses_age, event_n, Area)%>%
    summarise(tau = mean(tau, na.rm = TRUE),
              trials = mean(trials, na.rm = T),
              pernan = mean(pernan, na.rm = T))
  
  colnames(dataset)   <- c("suj", "Gender", "sesage","event_n", "Area", "tau", "trials", "pernan") 
  dataset_Area        <- dataset%>%group_by(suj, Area)%>%summarise(tau = mean(tau, na.rm = T), 
                                                                   trials = mean(trials,  na.rm = T),
                                                                   pernan = mean(pernan, na.rm = T)) # Summarise tau by Area including the percentage of non-convergent models
  
  dataset = merge(dataset, covs, by = "suj", all.x = T) # Merge the dataset with the covariates
  
  
  dataset$time    <- (dataset$sesage-6)/12 # Center the age for the intercept set to 0
  dataset$time_sq <- dataset$time^2 # time squared
  
  model_trials <- lmer(trials ~ time + time_sq + (1|suj), data = dataset%>%group_by(suj, time, time_sq)%>%summarise(trials = mean(trials,na.rm=T))) # Determine if there is an impact of age on the number of trials
  
  model.offset <- lmer(tau ~ pernan  + trials + 1  + (1|suj),    data = dataset, REML = TRUE) # Random intercept model
  model.slope  <- lmer(tau ~ pernan + trials + 1  + (time|suj), data = dataset, REML = TRUE) # Random slope model
  
  aic_compt <- aictab(c(model.offset, model.slope), modnames = c("Random Offset", "Random Slope")) # Compare models
  nametable <- paste(tablepath, "random_effects_devits_", cidx, ".html", sep = "") #Save the model comparison
  aictest      <- nice_table(aic_compt)
  flextable::save_as_html(aictest, path = nametable)
  
  if (aic_compt$Modnames[1] == "Random Slope") { # If the best model included random slope create the rest of models with that random structure 
    
    model.Area         <- lmer(tau ~ pernan + trials + Area  + (time|suj), data = dataset, REML = FALSE)
    model.time         <- lmer(tau ~ pernan + trials + time  + (time|suj), data = dataset, REML = FALSE)
    model.time_sq      <- lmer(tau ~ pernan + trials + time  + time_sq + (time|suj), data = dataset, REML = FALSE)
    model.time_sq_Area <- lmer(tau ~ pernan + trials + time  + time_sq + Area + (time|suj), data = dataset, REML = FALSE)
    model.inter_l      <- lmer(tau ~ pernan + trials + time*Area  + (time|suj), data = dataset, REML = FALSE)
    model.inter_sq     <- lmer(tau ~ pernan + trials + time*Area + time_sq  + (time|suj), data = dataset, REML = FALSE)
    
  } else { #Otherwise go for random offset 
    
    model.Area         <- lmer(tau ~ pernan + trials +Area  + (1|suj), data = dataset, REML = FALSE)
    model.time         <- lmer(tau ~ pernan + trials +time  + (1|suj), data = dataset, REML = FALSE)
    model.time_sq      <- lmer(tau ~ pernan + trials +time  + time_sq + (1|suj), data = dataset, REML = FALSE)
    model.time_sq_Area <- lmer(tau ~ pernan + trials +time  + time_sq + Area + (1|suj), data = dataset, REML = FALSE)
    model.inter_l      <- lmer(tau ~ pernan + trials +time*Area  + (1|suj), data = dataset, REML = FALSE)
    model.inter_sq     <- lmer(tau ~ pernan + trials +time*Area  + time_sq + (1|suj), data = dataset, REML = FALSE)
    
  }
  
  aic_compt <- aictab(c(model.Area, model.time, model.time_sq, model.time_sq_Area, model.inter_l, model.inter_sq), modnames = c("Area", "Time", "Time SQ", "Time SQ + Area", "Interaction", "Interaction SQ"))
  nametable <- paste(tablepath, "fixed_effects_devits_", cidx, ".html", sep = "")
  aictest      <- nice_table(aic_compt)
  flextable::save_as_html(aictest, path = nametable) #Compare all models and save the table
  
  if (aic_compt$Modnames[1] == "Interaction SQ"){ # Find the best model and save the Area post-hoc tests in case Area effect was significant
    model.final <- model.inter_sq
    nametable   <- paste(tablepath, "Area_effect_devits_", cidx, ".html", sep = "")
    model.Area           <- lmer(tau ~ pernan + trials +Area + (1|suj), data = dataset_Area, REML = FALSE) # Area model
    Area_test            <- glht(model.Area, linfct = mcp(Area = "Tukey"), test = adjusted("hommel")) # Extract the post hoc test
    Area_test            <-  as.data.frame(table_glht(Area_test))
    Area_test$comparison <- rownames(Area_test)
    flextable::save_as_html(nice_table(Area_test), path = nametable)
    
    six_to_nine     <- lmer(tau ~ time*Area + pernan + trials + (1|suj), data = filter(dataset, sesage<16), REML = TRUE)
    nine_to_sixteen <- lmer(tau ~ I(time - min(time))*Area + pernan + trials + (1|suj), data = filter(dataset, sesage>6), REML = TRUE)
    
  } else if (aic_compt$Modnames[1] == "Interaction"){
    model.final <- model.inter_l
    nametable <- paste(tablepath, "Area_effect_devits_", cidx, ".html", sep = "")
    model.Area           <- lmer(tau ~ pernan + trials +Area + (1|suj), data = dataset_Area, REML = FALSE)
    Area_test            <- glht(model.Area, linfct = mcp(Area = "Tukey"), test = adjusted("hommel"))
    Area_test            <-  as.data.frame(table_glht(Area_test))
    Area_test$comparison <- rownames(Area_test)
    flextable::save_as_html(nice_table(Area_test), path = nametable)
  } else if (aic_compt$Modnames[1] == "Time SQ + Area"){
    model.final <- model.time_sq_Area
    nametable <- paste(tablepath, "Area_effect_devits_", cidx, ".html", sep = "")
    model.Area           <- lmer(tau ~ pernan + trials +Area + (1|suj), data = dataset_Area, REML = FALSE)
    Area_test            <- glht(model.Area, linfct = mcp(Area = "Tukey"), test = adjusted("hommel"))
    Area_test            <-  as.data.frame(table_glht(Area_test))
    Area_test$comparison <- rownames(Area_test)
    flextable::save_as_html(nice_table(Area_test), path = nametable)
    six_to_nine     <- lmer(tau ~ time + pernan + trials + Area + (1|suj), data = filter(dataset, sesage<16), REML = TRUE)
    nine_to_sixteen <- lmer(tau ~ time + pernan + trials + Area + (1|suj), data = filter(dataset, sesage>6), REML = TRUE)
  } else if (aic_compt$Modnames[1] == "Time SQ") {
    model.final <- model.time_sq
    six_to_nine     <- lmer(tau ~ time + pernan + trials + (1|suj), data = filter(dataset, sesage<16), REML = TRUE)
    nine_to_sixteen <- lmer(tau ~ I(time - min(time)) + pernan + trials + (1|suj), data = filter(dataset, sesage>6), REML = TRUE)
    
  } else if (aic_compt$Modnames[1] == "Time") {
    model.final <- model.time
  } else {
    model.final    <- model.Area
    model.finalcov <- model.Areacov
    nametable <- paste(tablepath, "Area_effect_devits_", cidx, ".html", sep = "")
    model.Area           <- lmer(tau ~ pernan + trials +Area + (1|suj), data = dataset_Area, REML = FALSE)
    Area_test            <- glht(model.Area, linfct = mcp(Area = "Tukey"), test = adjusted("hommel"))
    Area_test            <-  as.data.frame(table_glht(Area_test))
    Area_test$comparison <- rownames(Area_test)
    flextable::save_as_html(nice_table(Area_test), path = nametable)
  }
  
  
  # This lines of code were added after the initial analysis knowing which was the best model to run the final models with covariates
  
  if (cidx == "exploratory") {
    model.final    <- lmer(tau ~ pernan + trials + time*Area  + time_sq + (time|suj), data = dataset, REML = TRUE)
    model.finalcov <- lmer(tau ~ pernan + trials + time*Area  + time_sq + MEdu + FEdu + MOcc + FOcc + I2N + Gender + gw + bw + BDI + (time|suj), data = dataset, REML = TRUE)
    
  } else { 
    model.final    <- lmer(tau ~ pernan + trials + time + Area  + time_sq + (time|suj), data = dataset, REML = TRUE)
    model.finalcov <- lmer(tau ~ pernan + trials + time + Area  + time_sq + MEdu + FEdu + MOcc + FOcc + I2N + Gender + gw + bw + BDI + (time|suj), data = dataset, REML = TRUE)
    }
  
  tabletitle = paste("Final Model - Tau Development")
  tablepathf  = paste(tablepath, "final_model_devits_", cidx, ".html", sep = "")
  print(tab_model(model.final, show.df = TRUE, show.intercept = TRUE, show.est = TRUE, show.se= T, show.p = TRUE, show.stat = TRUE, show.ngroups = TRUE, show.aic = TRUE, show.obs = TRUE, show.std = T, p.style = "scientific_stars", p.val = "satterthwaite", title = tabletitle, file = tablepathf ))
  
  tabletitle = paste("Final Model W/Covs - Tau Development")
  tablepathf  = paste(tablepath, "final_model_devits", cidx, "_covs.html", sep = "")
  print(tab_model(model.finalcov, show.df = TRUE, show.intercept = TRUE, show.est = TRUE, show.se= T, show.p = TRUE, show.stat = TRUE, show.ngroups = TRUE, show.aic = TRUE, show.obs = TRUE, show.std = T, p.style = "scientific_stars", p.val = "satterthwaite", title = tabletitle, file = tablepathf ))
  
  
  tabletitle = paste("CI Boots of the Final Model W/Covs - Tau Development")
  tablepathf  = paste(tablepath, "final_model_devits_ciboot_", cidx, ".html", sep = "")
  print(tab_model(model.final, iterations = 1000, bootstrap = TRUE, show.df = TRUE, show.intercept = TRUE, show.est = TRUE, show.p = TRUE, show.stat = TRUE, show.ngroups = TRUE, show.aic = TRUE, show.obs = TRUE, p.style = "scientific_stars", p.val = "satterthwaite", title = tabletitle, file = tablepathf))
  
  
  tabletitle = paste("CI Boots of the Final Model - Tau Development")
  tablepathf  = paste(tablepath, "final_model_devits_ciboot_", cidx, "_covs.html", sep = "")
  print(tab_model(model.finalcov, iterations = 1000, bootstrap = TRUE, show.intercept = TRUE, show.est = TRUE, show.aic = TRUE, show.se = T, p.style = "scientific_stars", p.val = "satterthwaite", file = tablepathf))
  
  tabletitle = paste("Final Model Wo/Covs - Tau Development 6 to 9 mo")
  tablepathf  = paste(tablepath, "final_model_devits", cidx, "_wocovs_6_to_9.html", sep = "")
  print(tab_model(six_to_nine, show.df = TRUE, show.intercept = TRUE, show.est = TRUE, show.se= T, show.p = TRUE, show.stat = TRUE, show.ngroups = TRUE, show.aic = TRUE, show.obs = TRUE, show.std = T, p.style = "scientific_stars", p.val = "satterthwaite", title = tabletitle, file = tablepathf ))
  
  tabletitle = paste("Final Model Wo/Covs - Tau Development 6 to 9 mo")
  tablepathf  = paste(tablepath, "final_model_devits", cidx, "_wocovs_6_to_9_boots.html", sep = "")
  print(tab_model(six_to_nine, iterations = 1000, bootstrap = TRUE, show.df = TRUE, show.intercept = TRUE, show.est = TRUE, show.p = TRUE, show.stat = TRUE, show.ngroups = TRUE, show.aic = TRUE, show.obs = TRUE, p.style = "scientific_stars", p.val = "satterthwaite", title = tabletitle, file = tablepathf))
  
  tabletitle = paste("Final Model Wo/Covs - Tau Development 9 to 16 mo")
  tablepathf  = paste(tablepath, "final_model_devits", cidx, "_wocovs_9_to_16.html", sep = "")
  print(tab_model(nine_to_sixteen, show.df = TRUE, show.intercept = TRUE, show.est = TRUE, show.se= T, show.p = TRUE, show.stat = TRUE, show.ngroups = TRUE, show.aic = TRUE, show.obs = TRUE, show.std = T, p.style = "scientific_stars", p.val = "satterthwaite", title = tabletitle, file = tablepathf ))
  
  tabletitle = paste("Final Model Wo/Covs - Tau Development 9 to 16 mo")
  tablepathf  = paste(tablepath, "final_model_devits", cidx, "_wocovs_9_to_16_boots.html", sep = "")
  print(tab_model(nine_to_sixteen, iterations = 1000, bootstrap = TRUE, show.df = TRUE, show.intercept = TRUE, show.est = TRUE, show.p = TRUE, show.stat = TRUE, show.ngroups = TRUE, show.aic = TRUE, show.obs = TRUE, p.style = "scientific_stars", p.val = "satterthwaite", title = tabletitle, file = tablepathf))
  
  
}