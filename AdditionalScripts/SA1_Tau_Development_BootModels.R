library(dplyr)
library(lme4)
library(lmerTest)
library(readr)
library(effectsize)
library(reshape2)
library(ggplot2)
library(ggpubr)

# Set the seed for reproducibility
set.seed(42)
nboots <- 5000

# Function to extract statistics from the model
stats_models <- function(model) {
  
  fixed_effects <- coef(summary(model))
  
  # Intercept
  t_values_int       <- fixed_effects[1, "t value"]
  beta_estimates_int <- fixed_effects[1, "Estimate"]
  p_values_int       <- fixed_effects[1, "Pr(>|t|)"]
  se_values_int      <- fixed_effects[1, "Std. Error"]
  
  # Age
  t_values_age <- fixed_effects[2, "t value"]
  beta_estimates_age <- fixed_effects[2, "Estimate"]
  p_values_age <- fixed_effects[2, "Pr(>|t|)"]
  se_values_age <- fixed_effects[2, "Std. Error"]
  
  # Age^2
  t_values_age2 <- fixed_effects[3, "t value"]
  beta_estimates_age2 <- fixed_effects[3, "Estimate"]
  p_values_age2 <- fixed_effects[3, "Pr(>|t|)"]
  se_values_age2 <- fixed_effects[3, "Std. Error"]
  
  stats <- data.frame(beta_estimates_int = beta_estimates_int, 
                      se_values_int = se_values_int,
                      t_values_int = t_values_int, 
                      p_values_int = p_values_int, 
                      beta_estimates_age = beta_estimates_age, 
                      se_values_age = se_values_age,
                      t_values_age = t_values_age, 
                      p_values_age = p_values_age, 
                      beta_estimates_age2 = beta_estimates_age2, 
                      se_values_age2 = se_values_age2,
                      t_values_age2 = t_values_age2, 
                      p_values_age2 = p_values_age2)
  return(stats)
}

# Function to sampling the trials for the bootstrapping
sample_trials <- function(data) {
  data %>%
    group_by(subject, ses_age, area) %>%
    sample_n(3, replace = T) %>%
    ungroup()
}

# Where is the data
tablepath <- "~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/"

# Load the covariates and center them
covdata   <- read_csv(paste(tablepath, "covariates_its.csv", sep = ""))%>%
  mutate(bw = scale(bw, center = T, scale = T), 
         gw = scale(gw, center = T, scale = T), 
         BDI = scale (BDI, center = T, scale = T)) # BDI = Beck Depression Inventory

# The cohorts and covariates to loop through
cohorts = c("Exploratory",  "Validation")
covs   = c("Basic", "Complete")
anl    = c("No Boot", "Boot")

loopidx = 1 
for (cohort in cohorts){
  
  if (cohort == "Exploratory") {
    taudata <- read_csv("Desktop/Tau_data/tau_development_exploratory_v4.csv")%>%
      rename(sex = Gender, area = Area, channel_number = channel, channel = Label)
    
  } else {
    taudata <- read_csv("Desktop/Tau_data/tau_development_validation_v4.csv")%>%
      rename(sex = Gender, area = Area, channel_number = channel, channel = Label)}
  
  for (covidx in covs){

    for (aidx in anl){
      # Average the dataset excluding epochs with low convergence and artifacted epochs
      dataset <- taudata %>%
        filter(pernan < .25, event != 9999) %>%
        mutate(trials = mean(n()) / length(unique(channel)), .by = c(subject, ses_age)) %>%
        group_by(subject, sex, ses_age, event_n, area) %>%
        summarise(tau = mean(tauimp, na.rm = TRUE), 
                  trials = mean(trials, na.rm = TRUE), 
                  pernan = mean(pernan, na.rm = TRUE)) %>%
        mutate(time = (ses_age - 6) / 12) # Center the age for the lmm models so the intercept corresponds with age 6 months
      
      if (covidx == "Complete"){
        dataset  = merge(dataset, covdata, by.x = "subject", by.y = "suj") # In case we selected completed covariates, we merge the data with the covariates list
      }
      
      if (aidx == "Boot"){
        
        
        if (covidx == "Basic"){ # The function depends on the covariates 
          fit_lmer <- function(data) {
            lmer(tau ~ time + I((time^2)) + area + pernan + (time|subject), data = data, REML = TRUE)
          }
        } else {
          fit_lmer <- function(data){
            lmer(tau ~ time + I((time^2)) + area + pernan + MEdu + FEdu + MOcc + FOcc + I2N + sex + gw + bw + BDI  + (time|subject), data = data, REML = TRUE)}
        }
        # Perform bootstrapping
        bootstrap_samples <- replicate(nboots, sample_trials(dataset), simplify = FALSE) # 5000 boots to sample the trials
        
        # Fit the model on each bootstrapped sample
        bootstrap_models <- lapply(bootstrap_samples, fit_lmer)
        
        
        # Apply the function to all bootstrapped models
        bootstraps_stats <- lapply(bootstrap_models, stats_models)
        
        # Combine the results into one data frame
        combined_stats <- bind_rows(bootstraps_stats, .id = "bootstrap_id") %>%
          reshape2::melt(id.vars = "bootstrap_id") %>%
          mutate(fixed_effect = if_else(grepl("int", as.character(variable)), "Intercept", 
                                        if_else(grepl("age2", as.character(variable)), "AgeSQ", "Age")), 
                 stat = if_else(grepl("t_value", as.character(variable)), "t", 
                                if_else(grepl("beta", as.character(variable)), "B", if_else(grepl("se", as.character(variable)), "SE", "p-val")))) %>%
          dplyr::select(-variable)
        
        # Plot the distribution of the parameters for each fixed effect and statistic
        plot_data <- combined_stats %>%
          filter(fixed_effect != "Intercept") %>%
          ungroup() %>%
          group_by(bootstrap_id, fixed_effect) %>%
          mutate(genp = value[stat == "p-val"], 
                 signif = if_else(genp < 0.001, 0, if_else(genp < .01, 1, if_else(genp < .05, 2, 3))))
        
        # Plot the distribution of the parameters for each fixed effect and statistic adding the p values 
        parameters_plot = ggplot(plot_data, aes(x = value, fill = factor(signif, levels = c(3,2,1,0), labels = c("P>.05", "P<.05", "P<.01", "P<.001")))) +
          geom_histogram(alpha = 0.3, bins = 25) + 
          facet_grid(fixed_effect ~ stat, scales = "free") +
          ggpubr::theme_pubr() + labs(y = "Frequency", 
                                      x = "Parameter Value", 
                                      fill = "", 
                                      title = "LMM Parameters Bootstrapping the Trials", 
                                      subtitle = paste(covidx, " covariates / ", cohort, " Cohort", sep = "")) + theme(legend.position = "bottom")
        ggsave(parameters_plot, filename = paste(tablepath, "descriptive_stat_boottrials_", cohort, "_", covidx, ".svg", sep = ""))
        rm(parameters_plot)
        rm(plot_data)
        
        # Summarise the results
        model_stats <- combined_stats %>%
          dplyr::select(-bootstrap_id) %>%
          group_by(fixed_effect) %>%summarise(B = mean(value[stat=="B"]),
                                              SE = mean(value[stat=="SE"]), 
                                              Cil = quantile(value[stat=="B"], .025, na.rm = T),
                                            Ciu = quantile(value[stat=="B"], .975, na.rm = T))%>%
          mutate(cohortname   = cohort, 
                 covsinclu    = covidx, 
                 boottype     = aidx)
        
      } else {
        
        if (covidx == "Basic"){
          fit_lmer <- function(data) {
            lmer(tau ~ time + I((time^2)) + area + pernan +  trials + (time|subject), data = data, REML = TRUE)
          }
        } else {
          fit_lmer <- function(data){
            lmer(tau ~ time + I((time^2)) + area + pernan + MEdu + FEdu + MOcc + FOcc + I2N + sex + gw + bw + BDI + trials + (time|subject), data = data, REML = TRUE)}
        }
        
        model = fit_lmer(dataset)
        
        # Define the function to extract the fixed effects
        fixed_effects <- function(model) {
          fixef(model)
        }
        
        boot_results <- bootMer(model, FUN = fixed_effects, nsim = nboots, use.u = TRUE)
        
        boot_ci <- function(boot_results, index) {
          resampled_stats <- boot_results$t[, index]
          original_stat <- boot_results$t0[index]
          quantile(resampled_stats, probs = c(0.025, 0.975)) - original_stat
        }
        
        ci_df <- data.frame(
          Parameter = names(fixef(model)),
          CI_lower = sapply(1:length(fixef(model)), function(i) fixef(model)[i] + boot_ci(boot_results, i)[1]),
          CI_upper = sapply(1:length(fixef(model)), function(i) fixef(model)[i] + boot_ci(boot_results, i)[2])
        )
        
        
        model_stats = stats_models(model)%>%
          melt()%>%mutate(fixed_effect = if_else(grepl("int", as.character(variable)), "Intercept", 
                                                 if_else(grepl("age2", as.character(variable)), "AgeSQ", "Age")), 
                          stat = if_else(grepl("t_value", as.character(variable)), "t", 
                                         if_else(grepl("beta", as.character(variable)), "B", if_else(grepl("se", as.character(variable)), "SE", "p-val")))) %>%
          dplyr::select(-variable)%>%
          group_by(fixed_effect) %>%summarise(B = mean(value[stat=="B"]),
                                              SE = mean(value[stat=="SE"]), 
                                              Cil = NA,
                                              Ciu = NA)
        
        model_stats$Cil[model_stats$fixed_effect=="Intercept"] = ci_df$CI_lower[ci_df$Parameter=="(Intercept)"]
        model_stats$Ciu[model_stats$fixed_effect=="Intercept"] = ci_df$CI_upper[ci_df$Parameter=="(Intercept)"]
        
        model_stats$Cil[model_stats$fixed_effect=="Age"] = ci_df$CI_lower[ci_df$Parameter=="time"]
        model_stats$Ciu[model_stats$fixed_effect=="Age"] = ci_df$CI_upper[ci_df$Parameter=="time"]
        
        model_stats$Cil[model_stats$fixed_effect=="AgeSQ"] = ci_df$CI_lower[ci_df$Parameter=="I((time^2))"]
        model_stats$Ciu[model_stats$fixed_effect=="AgeSQ"] = ci_df$CI_upper[ci_df$Parameter=="I((time^2))"]
        
        model_stats = model_stats%>%
          mutate(cohortname   = cohort, 
                 covsinclu    = covidx, 
                 boottype     = aidx)
      }
      
      
      if (loopidx ==1) {
        
        complete_summary = model_stats
      } else {
        
        complete_summary = rbind(complete_summary, model_stats)
      }
      
      loopidx = loopidx + 1 
      
      
    } 
  }
}

estimate_tables = nice_table(complete_summary)
path2table      = paste(tablepath, "boots_results.docx", sep="")
save_as_docx(estimate_tables,path = path2table)

complete_summary$boottype  = factor(complete_summary$boottype, levels = c("No Boot", "Boot"), labels = c("All Trials", "Trials Boot"))
complete_summary$covsinclu = factor(complete_summary$covsinclu, levels = c("Basic", "Complete"), labels = c("w/o Covs", "w/ Covs"))

ggplot(complete_summary, aes(x = B, y = interaction(covsinclu, boottype))) + geom_point(size = 3, shape = 23, fill = "black") + 
  geom_segment(aes(x = B - SE, xend = B + SE, y = interaction(covsinclu, boottype), yend = interaction(covsinclu, boottype)), size = 1.5) +
  geom_segment(aes(x = Cil, xend = Ciu, y = interaction(covsinclu, boottype), yend = interaction(covsinclu, boottype)), size = .75) + facet_grid(fixed_effect~cohortname) +
  geom_vline(xintercept = 0) + ggpubr::theme_pubr() + 
  labs(title = "Beta estimates", 
       y = "Model")
ggsave(filename = paste(tablepath, "beta_estimates_across_models.svg", sep = ""), width = 8, height = 4)


complete_summary2 <- complete_summary%>%mutate(across(c(B, SE, Cil, Ciu), ~round(.x, 2), .names = "{.col}"))%>%
  mutate(cov_boot = paste(covsinclu, ".", boottype, sep=""), 
         estats   = paste(B, " (", SE, ")", " [", Cil, " - ", Ciu,"]", sep = ""))%>%
  dplyr::select(cohortname, fixed_effect, cov_boot, estats)%>%
  dcast(cohortname + fixed_effect ~ cov_boot)

estimates_table = nice_table(complete_summary2, separate.header = T)
path2table      = paste(tablepath, "boots_results_sorted.docx", sep="")
save_as_docx(estimates_table,path = path2table)