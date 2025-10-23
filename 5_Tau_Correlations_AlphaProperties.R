# ITS Tau correlation with the oscillatory and rhythmic alpha variables ---------------------------

## Load the required packages -------------------------------------------------
required_packages <- c(
  "tidyverse", # Load this first because it includes ggplot2 and dplyr
  "ggpubr",
  "psych",
  "reshape2",
  "flextable",
  "rempsyc",
  "ppcor",
  "boot")

new_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(new_packages) > 0) {
  install.packages(new_packages)
}

sapply(required_packages, library, character.only = TRUE)


## Functions and plot settings ------------------------------------------------

## Plot themes -----------------
settheme =  theme(strip.text.x = element_text(size = 9.5, face = "bold"), strip.text.y = element_text(size = 9.5, face = "bold")) + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), title = element_text(size = 10), axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold")) 

## Partial correlation function for bootstrapping-----------------------
partial_corr <- function(data, indices) {
  cormat <- data[indices, ] # Resample with replacement
  pcor_results <- ppcor::pcor.test(cormat$tau, cormat$voi, cormat[c("pernan", "rsq")], method = "spearman") # alpha oscillatory power also controls for R2
  return(pcor_results$estimate)}  

partial_corr_lag <- function(data, indices) {
  cormat <- data[indices, ] # Resample with replacement
  pcor_results <- ppcor::pcor.test(cormat$tau, cormat$voi, cormat[c("pernan")], method = "spearman")
  return(pcor_results$estimate)} 

#set seed for reproducibility 
set.seed(42)

## Load the data ---------------------------------------------------------------
## Data parameters -------------------------------------------------------------
cohorts = c('exploratory', 'validation')

# Where do you want to save the results?
tablepath = "~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/"

## Tau data load, preparation, and analysis ------------------------------------
for (cohort in cohorts) {
  if (cohort == "exploratory") {
    taudata <- read_csv("Desktop/Tau_data/tau_development_exploratory_v4.csv")%>%
      rename(sex = Gender, area = Area, channel_number = channel, channel = Label)

    } else {
    taudata <- read_csv("Desktop/Tau_data/tau_development_validation_v4.csv")%>%
      rename(sex = Gender, area = Area, channel_number = channel, channel = Label)}

# Exclude 9999 events (social/high movement) and the epochs with low convergence
taudata_space <- taudata%>%filter(pernan < .25, event != 9999)%>%
  group_by(subject, sex, ses_age, channel, area)%>% # Average the remaining epochs
  summarise(tau    = mean(tauimp, na.rm = T),
            pernan = mean(is.na(tauorg), na.rm = T))%>%filter(pernan < .25)

taudata_space    <- dplyr::select(taudata_space, subject, sex, ses_age, channel, area, tau, pernan) #select only the relevant variables

taudata <- taudata%>%filter(pernan < .25, event != 9999)%>% # Create a dataset with one value per electrode and area 
  mutate(trials = mean(n()) / length(unique(channel)), .by = c(subject, ses_age))%>%
  group_by(subject, sex, ses_age, channel, area)%>%
  summarise(tau    = mean(tauimp, na.rm = T),
            trials = mean(trials, na.rm = T),
            pernan = mean(pernan, na.rm = T))

colnames(taudata)       <- c("suj", "sex", "sesage", "elect", "area", "tau","trials", "pernan") 
colnames(taudata_space) <- c("suj", "sex", "sesage", "elect", "area", "tau","pernan")

## 1) Oscillatory alpha activity ----------------------------------------------
###   Load the oscillatory alpha activity data --------------------------------
pow                 <- read_csv("~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/pow_20Hz_6m_16m_ITS.csv")

pow                 <- filter(pow, band == "Alpha" & rsquared > .949) # Select only the alpha band and the epochs with good R2


fit_descriptives <- pow%>%filter(incluelect == 1, !is.na(trials), develop_f == 1, suj %in% taudata$suj)%>%
  group_by(sesage)%>%
  summarise(mtrials = mean(trials, na.rm = T),
            sdtrials = sd(trials, na.rm = T),
            mrsquared = mean(rsquared, na.rm = T),
            sdrsquared = sd(rsquared, na.rm = T))

pow               <- pow%>%dplyr::select(suj, sex, sesage, elect, area, oscpow, abspow, relpow, slope, offset, rsquared, peak, freq)%>%
  group_by(suj, sex, sesage, elect, area)%>%
  summarise(osc    = mean(oscpow, na.rm= TRUE), 
            peak   = mean(peak, na.rm = TRUE),
            freq   = mean(freq, na.rm = TRUE),
            rsq    = mean(rsquared, na.rm = T))

descriptives        <- pow%>%filter(suj %in% taudata$suj) # model fit, power, and peak frequency descriptives
descriptives        <- dplyr::select(descriptives, suj, sex, sesage, area, osc, peak, freq, rsq)%>%
  group_by(suj, sex, sesage, area)%>%
  summarise(osc = mean(osc, na.rm = T), 
            freq= mean(freq[peak==1], na.rm = T), 
            rsq = mean(rsq, na.rm = T))

descriptives <- melt(descriptives, id.vars = c("suj", "sex", "sesage", "area"))%>% # parse the descriptives to more standardized name
  group_by(sesage, sex, variable, area)%>%
  summarise(m  = mean(as.double(value), na.rm = T), 
            sd =sd(value,na.rm = T))

descriptives$m  <- round(descriptives$m, 2) # round descriptives to the second decimal point
descriptives$sd <- round(descriptives$sd, 2)

descriptives$stat <- paste(descriptives$m, " (", descriptives$sd, ")", sep ="") # create a statistic column with the mean and sd as m (sd)

descriptives      <- descriptives%>%dplyr::select(-m,-sd) # remove the mean and sd columns

desctable         <- dcast(descriptives, variable + sesage + sex ~  area) # reshape the descriptive table

desctable <- nice_table(desctable, separate.header = F) # Formated table

path2table = paste(tablepath, "oscpow_descriptives_", cohort, ".docx", sep="")
save_as_docx(desctable,path = path2table) # save the table

pow <- melt(pow, id.vars = c("suj", "sex", "sesage", "elect", "area", "peak"))

varnames <- as.vector(unique(pow$variable))[1:2] # frequency and oscillatory power
ses      <- c(6,9,16)

### Within-participant correlations between alpha oscillatory activity and tau --------------------
session       = NA
variable      = NA
rs            = NA
n             = NA
ci_min        = NA
ci_max        = NA
pval          = NA


sesidx  = 0 
loopidx = 1 
for (age in ses){
  sesidx = sesidx + 1 
  varidx = 0 
  
  data1 <- taudata%>%filter(sesage == age)%>% # prepare tau data filtering only the relevant age
    group_by(suj, sex)%>%summarise(tau    = mean(tau, na.rm = TRUE),
                                   trials = mean(trials, na.rm = TRUE),
                                   pernan = mean(trials, na.rm = TRUE))
  
  for (items in varnames) {
    varidx = varidx + 1  
    
    if (varnames[varidx] != "freq"){ # we divide into frequency and oscillatory because frequency can only be computed in electrodes that had an oscillatory peak (i.e., peak == 1)
      data2  <- pow%>%filter(sesage == ses[sesidx] & (variable == varnames[varidx] | variable == "rsq"))%>% # select the variable of interest and the rsq
        group_by(suj, variable)%>%summarise(M = mean(value, na.rm = T))%>%dcast(suj ~ variable)
      colnames(data2) <- c("suj", "voi", "rsq")  # rename the variable of interest to "voi" so it matches the general function below
    } else {
      data2  <- pow%>%filter(sesage == ses[sesidx] & (variable == varnames[varidx] | variable == "rsq") & peak == 1)%>%
        group_by(suj, variable)%>%summarise(M = mean(value, na.rm = T))%>%dcast(suj ~ variable)
      colnames(data2) <- c("suj", "voi", "rsq")  
      
    }
    
    cormat <-merge(data1,data2, by = "suj") #merge tau and oscpow with tau data
    cor_results <- pcor.test(cormat$tau, cormat$voi, cormat[c("pernan", "rsq")], method = "spearman") #partial correlation controlling for rsq and percentage of nan electrodes
    
    boot_results <- boot(cormat, partial_corr, R = 5000) # bootstrap the partial correlation to create the confidence intervals
    ci           <- boot.ci(boot_results, type = "perc") # save the CI
    
    #store the results
    session[loopidx]       = age
    variable[loopidx]      = items
    rs[loopidx]            = cor_results$estimate
    n[loopidx]             = cor_results$n
    ci_min[loopidx]        = ci$percent[4]
    ci_max[loopidx]        = ci$percent[5]
    pval[loopidx]          = cor_results$p
    loopidx = loopidx + 1 
  }
}

session  <- c(session)
variable <- c(variable)
rs       <- c(rs)
n        <- c(n)
ci_min   <- c(ci_min)
ci_max   <- c(ci_max)
pval     <- c(pval)

pval <- p.adjust(pval, method = "fdr") # correct with FDR
results_correlation <- data.frame(n, session, variable, rs, pval, ci_min, ci_max, stringsAsFactors = TRUE) # create data frame with the results


results_correlation$significance[results_correlation$pval < .001] <- "***" # significance levels after FDR
results_correlation$significance[results_correlation$pval < .01 & results_correlation$pval > .001] <- "**"
results_correlation$significance[results_correlation$pval < .05 & results_correlation$pval > .01]  <- "*"
results_correlation$significance[results_correlation$pval < .1 & results_correlation$pval > .05]  <- ""
results_correlation$significance[results_correlation$pval > .1]  <- ""

results_correlation$session[results_correlation$session==6]  <- "6-mo." #rename the ages
results_correlation$session[results_correlation$session==9]  <- "9-mo."
results_correlation$session[results_correlation$session==16] <- "16-mo."

results_correlation$session <- factor(results_correlation$session, levels = c("6-mo.", "9-mo.", "16-mo.")) #reorder the ages


results_table        <- dplyr::select(results_correlation, c("session", "n", "variable", "rs", "ci_min", "ci_max", "significance")) #select significant variables
results_table$rs     <- round(results_table$rs, 2) #round parameters
results_table$ci_min <- round(results_table$ci_min,2)
results_table$ci_max <- round(results_table$ci_max,2)

results_table$stats  <- paste(results_table$rs, results_table$significance, " [", results_table$ci_min, " - ", results_table$ci_max, "]", sep ="") #create a single correlation value as r-significance [ci]
results_table        <- dplyr::select(results_table, c("session", "variable", "n", "stats")) # Keep only the relevant columns

results_table                  <- dcast(results_table, session + n ~ variable, value.var = "stats") #format the table

path2table = paste(tablepath, "tau_correlation_powvars_", cohort, ".docx", sep ="") 
save_as_docx(nice_table(results_table),path = path2table) #save a formated table as docx


#### Within-participant correlation plots --------------------------------------
powplot <- pow%>% filter(variable != 'rsq')%>%
  group_by(suj, sesage)%>%
  summarise(freq = mean(value[variable == "freq" & peak == 1], na.rm = TRUE),
            osc  = mean(value[variable == "osc"], na.rm = TRUE))

powplot <- melt(powplot, id.vars = c("suj", "sesage")) # reshape to long format

tauplot <- taudata%>% # prepare tau data filtering only the relevant age
  group_by(suj, sesage)%>%summarise(tau    = mean(tau, na.rm = TRUE),
                                 trials = mean(trials, na.rm = TRUE),
                                 pernan = mean(trials, na.rm = TRUE))

plotdataset <- merge(powplot, tauplot, by = c("suj", "sesage"))

plotdataset$sesage[plotdataset$sesage == 6]  <- "6-mo." 
plotdataset$sesage[plotdataset$sesage == 9]  <- "9-mo." 
plotdataset$sesage[plotdataset$sesage == 16] <- "16-mo." 

path2plot <- paste(tablepath, "individualcor_its_osc_", cohort, ".svg", sep ="")
plot_osc = ggplot(filter(plotdataset, variable == "osc"), aes(y=tau, x=value)) +
  geom_smooth(method = "lm", se = T, color = "black", alpha = .30, fill = "grey75") + ggpubr::theme_pubr() +
  geom_point(aes(color = "lightpink")) + 
  facet_grid(.~factor(sesage, levels = c("6-mo.", "9-mo.", "16-mo.")), scales = "free") + 
  theme(legend.position = "none") +
  settheme +
  xlab("Oscillatory Power") + 
  ylab("Tau (s)") 

plot_freq = ggplot(filter(plotdataset, variable == "freq"), aes(y=tau, x=value)) +
  geom_smooth(method = "lm", se = T, color = "black", alpha = .30, fill = "grey75") + ggpubr::theme_pubr() +
  geom_point(aes(color = "lightpink")) + 
  facet_grid(.~factor(sesage, levels = c("6-mo.", "9-mo.", "16-mo.")), scales = "free") + 
  theme(legend.position = "none") +
  settheme +
  xlab("Peak Frequency (Hz)") + 
  ylab("Tau (s)") 

ggpubr::ggarrange(plot_osc, plot_freq, nrow = 2)
ggsave(path2plot, width = 14.2, height = 11, units = "cm")

### Space correlations between alpha oscillatory actiivty and tau --------------
varnames <- as.vector(unique(pow$variable))[1:2]
ses      <- c(6,9,16)

session       = NA
variable      = NA
rs            = NA
n             = NA
ci_min        = NA
ci_max        = NA
pval          = NA

loopidx = 1 
for (sesi in ses){
  
  data1 <- taudata_space%>%filter(sesage == sesi)%>% #Now we reorder by electrode
    group_by(elect)%>%summarise(tau    = mean(tau, na.rm = TRUE),
                                pernan = mean(pernan, na.rm = TRUE))
  
  for (vari in varnames) {
    varidx = varidx + 1  
    
    if (vari != "freq"){
      data2  <- pow%>%filter(sesage == sesi & (variable == vari | variable == "rsq"))%>%
        group_by(elect, variable)%>%summarise(M = mean(value, na.rm = T))%>%dcast(elect ~ variable)
      colnames(data2) <- c("elect", "voi", "rsq")  
    } else {
      data2  <- pow%>%filter(sesage == sesi & (variable == vari | variable == "rsq") & peak == 1)%>%
        group_by(elect, variable)%>%summarise(M = mean(value, na.rm = T))%>%dcast(elect ~ variable)
      colnames(data2) <- c("elect", "voi", "rsq")  
      
    }
    
    cormat <-merge(data1,data2, by = "elect")
    cor_results <- pcor.test(cormat$tau, cormat$voi, cormat[c("pernan", "rsq")], method = "spearman")
    
    boot_results <- boot(cormat, partial_corr, R = 5000)
    ci           <- boot.ci(boot_results, type = "perc")
    
    session[loopidx]       = sesi
    variable[loopidx]      = vari
    rs[loopidx]            = cor_results$estimate
    n[loopidx]             = cor_results$n
    ci_min[loopidx]        = ci$percent[4]
    ci_max[loopidx]        = ci$percent[5]
    pval[loopidx]          = cor_results$p
    loopidx = loopidx + 1 
    
  }
}

session  <- c(session)
varaible <- c(variable)
rs       <- c(rs)
n        <- c(n)
ci_min   <- c(ci_min)
ci_max   <- c(ci_max)
pval     <- c(pval)

pval <- p.adjust(pval, method = "fdr")
results_correlation <- data.frame(n, session, variable, rs, pval, ci_min, ci_max, stringsAsFactors = TRUE)


results_correlation$significance[results_correlation$pval < .001] <- "***"
results_correlation$significance[results_correlation$pval < .01 & results_correlation$pval > .001] <- "**"
results_correlation$significance[results_correlation$pval < .05 & results_correlation$pval > .01]  <- "*"
results_correlation$significance[results_correlation$pval < .1  & results_correlation$pval > .05]  <- ""
results_correlation$significance[results_correlation$pval > .1]  <- ""

results_correlation$session[results_correlation$session==6]  <- "6-mo."
results_correlation$session[results_correlation$session==9]  <- "9-mo."
results_correlation$session[results_correlation$session==16] <- "16-mo."

results_correlation$session <- factor(results_correlation$session, levels = c("6-mo.", "9-mo.", "16-mo."))


results_table        <- dplyr::select(results_correlation, c("session", "n", "variable", "rs", "ci_min", "ci_max", "significance"))
results_table$rs     <- round(results_table$rs, 2)
results_table$ci_min <- round(results_table$ci_min,2)
results_table$ci_max <- round(results_table$ci_max,2)

results_table$stats  <- paste(results_table$rs, results_table$significance, " [", results_table$ci_min, " - ", results_table$ci_max, "]", sep ="")
results_table        <- dplyr::select(results_table, c("session", "variable", "n", "stats"))

results_table                  <- dcast(results_table, session + n ~ variable, value.var = "stats")

path2table = paste(tablepath, "tau_correlation_powvars_space_", cohort, ".docx", sep ="")
save_as_docx(nice_table(results_table),path = path2table)


#### Space correlation plots  --------------------------------------------------
powplot <- pow%>% filter(variable != 'rsq')%>%
  group_by(sesage, elect, area)%>%
  summarise(freq = mean(value[variable == "freq" & peak == 1], na.rm = TRUE),
            osc  = mean(value[variable == "osc"], na.rm = TRUE))

powplot <- melt(powplot, id.vars = c("sesage", "elect", "area")) # reshape to long format

tauplot <- taudata_space%>% 
  group_by(sesage, elect, area)%>%summarise(tau    = mean(tau, na.rm = TRUE))

plotdataset <- merge(powplot, tauplot, by = c("sesage", "elect", "area"))

plotdataset$sesage[plotdataset$sesage == 6]  <- "6-mo." 
plotdataset$sesage[plotdataset$sesage == 9]  <- "9-mo." 
plotdataset$sesage[plotdataset$sesage == 16] <- "16-mo." 

plotdataset$sesage <- factor(plotdataset$sesage, levels = c("6-mo.", "9-mo.", "16-mo."))

path2plot <- paste(tablepath, "spatial_stability_its_osc_", cohort, ".svg", sep ="")
plot_osc = ggplot(filter(plotdataset, variable == "osc"), aes(y=tau, x=value)) + geom_point(aes(color = area)) +
  geom_smooth(method = "lm", color = "black", se = T, alpha = .30, fill = "grey75") +
  scale_color_manual(values =c("#C4961A","#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")) + 
  theme_pubr() + facet_grid(.~factor(sesage, levels = c("6-mo.", "9-mo.", "16-mo.")), scales = "free") +
  settheme +
  labs(x = "Oscillatory Power", 
       y = "Tau (s)", 
       color = "Area")

plot_freq = ggplot(filter(plotdataset, variable == "freq"), aes(y=tau, x=value)) + geom_point(aes(color = area)) + 
  geom_smooth(method = "lm", color = "black", se = T, alpha = .30, fill = "grey75") +
  scale_color_manual(values =c("#C4961A","#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")) + 
  theme_pubr() + facet_grid(.~factor(sesage, levels = c("6-mo.", "9-mo.", "16-mo.")), scales = "free") +
  settheme +
  labs(x = "Frequency (Hz)", 
       y = "Tau (s)", 
       color = "Area")

ggarrange(plot_osc, plot_freq, nrow = 2, common.legend = T, legend = "bottom")
ggsave(path2plot, width = 16, height = 12, units = "cm")


## 2) Lagged coherence correlations --------------------------------------------
#Load the raw lagged coherence data and compute the metrics of interest
#pow                 <- read_csv("~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/pow_burst_rythm_6m_36m_ITS.csv") #Load the lagged coherence data
#pow                 <- filter(pow, band == "Alpha") # We keep only the alpha band
#pow                 <- pow%>%dplyr::select(suj, sex, sesage, elect, incluelect, area, lag, lag_pow)%>% # Relevant variables
  #group_by(suj, sex, sesage, elect, area)%>%
  #summarise(rhythmlag    = mean(lag_pow[lag>2.5&lag<3.5], na.rm= TRUE), # Burst lags
            #burstlag     = mean(lag_pow[lag<2.5], na.rm = TRUE)) # Rhythm lags

### Save for future use 
#write.csv(pow, "~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/lagged_coherence_alpha_6_16m_metrics.csv", row.names = F)

# Load the dataset with the lagged coherence data
pow <- read_csv("~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/lagged_coherence_alpha_6_16m_metrics.csv")

descriptives        <- pow%>%filter(suj %in% taudata$suj) 

descriptives <- melt(descriptives, id.vars = c("suj", "sex", "sesage", "elect", "area"))%>% #Descriptive table
  group_by(sesage, sex, variable, area)%>%
  summarise(m  = mean(as.double(value), na.rm = T), 
            sd =sd(value,na.rm = T))

descriptives$m  <- round(descriptives$m, 2)
descriptives$sd <- round(descriptives$sd, 2)

descriptives$stat <- paste(descriptives$m, " (", descriptives$sd, ")", sep ="")

descriptives      <- descriptives%>%dplyr::select(-m,-sd)

desctable         <- dcast(descriptives, variable + sesage + sex ~  area)


desctable <- nice_table(desctable, separate.header = F)

path2table = paste(tablepath, "lagcoh_descriptives_", cohort, ".docx", sep="")
save_as_docx(desctable,path = path2table)


pow <- melt(pow, id.vars = c("suj", "sex", "sesage", "elect", "area"))

varnames <- as.vector(unique(pow$variable))
ses      <- c(6,9,16)


#### Within-subject correlations between lagged coherence and tau  -------------
session       = NA
variable      = NA
rs            = NA
n             = NA
ci_min        = NA
ci_max        = NA
pval          = NA

loopidx = 1 
for (sesi in ses){
  
  data1 <- taudata%>%filter(sesage == sesi)%>%
    group_by(suj)%>%summarise(tau    = mean(tau, na.rm = TRUE),
                              pernan = mean(pernan, na.rm = TRUE))
  
  for (vari in varnames) {
    
    if (vari != "freq"){
      data2  <- pow%>%filter(sesage == sesi & variable == vari)%>%
        group_by(suj)%>%summarise(voi = mean(value, na.rm = TRUE))}
    else {
      data2 <- filter(pow, peak == 1)
      data2 <- data2%>%filter(sesage == sesi & variable == vari)%>%
        group_by(suj)%>%summarise(voi = mean(value, na.rm = TRUE))
    }
    
    cormat <-merge(data1,data2, by = "suj")
    
    cor_results <- pcor.test(cormat$tau, cormat$voi, cormat[c("pernan")], method = "spearman")
    
    
    boot_results <- boot(cormat, partial_corr_lag, R = 5000) # Different correlation function because we cannot control for R2 
    ci           <- boot.ci(boot_results, type = "perc")
    
    session[loopidx]       = sesi
    variable[loopidx]      = vari
    rs[loopidx]            = cor_results$estimate
    n[loopidx]             = cor_results$n
    ci_min[loopidx]        = ci$percent[4]
    ci_max[loopidx]        = ci$percent[5]
    pval[loopidx]          = cor_results$p
    loopidx = loopidx + 1 
  }
}

session  <- c(session)
varaible <- c(variable)
rs       <- c(rs)
n        <- c(n)
ci_min   <- c(ci_min)
ci_max   <- c(ci_max)
pval     <- c(pval)

pval <- p.adjust(pval, method = "fdr")
results_correlation <- data.frame(n, session, variable, rs, pval, ci_min, ci_max, stringsAsFactors = TRUE)


results_correlation$significance[results_correlation$pval < .001] <- "***"
results_correlation$significance[results_correlation$pval < .01 & results_correlation$pval > .001] <- "**"
results_correlation$significance[results_correlation$pval < .05 & results_correlation$pval > .01]  <- "*"
results_correlation$significance[results_correlation$pval < .1 & results_correlation$pval > .05]  <- ""
results_correlation$significance[results_correlation$pval > .1]  <- ""

results_correlation$session[results_correlation$session==6]  <- "6-mo."
results_correlation$session[results_correlation$session==9]  <- "9-mo."
results_correlation$session[results_correlation$session==16] <- "16-mo."

results_correlation$session <- factor(results_correlation$session, levels = c("6-mo.", "9-mo.", "16-mo."))


results_table        <- dplyr::select(results_correlation, c("session", "n", "variable", "rs", "ci_min", "ci_max", "significance"))
results_table$rs     <- round(results_table$rs, 2)
results_table$ci_min <- round(results_table$ci_min,2)
results_table$ci_max <- round(results_table$ci_max,2)

results_table$stats  <- paste(results_table$rs, results_table$significance, " [", results_table$ci_min, " - ", results_table$ci_max, "]", sep ="")
results_table        <- dplyr::select(results_table, c("session", "variable", "n", "stats"))

results_table                  <- dcast(results_table, session + n ~ variable, value.var = "stats")

path2table = paste(tablepath, "tau_correlation_brythm_", cohort, ".docx", sep ="")
save_as_docx(nice_table(results_table),path = path2table)

#### Within participants correlation plots  ------------------------------------

powplot <- pow%>%group_by(suj, sesage, variable)%>%
  summarise(value = mean(value, na.rm = TRUE))

tauplot <- taudata%>% # prepare tau data filtering only the relevant age
  group_by(suj, sesage)%>%summarise(tau    = mean(tau, na.rm = TRUE))

plotdataset <- merge(powplot, tauplot, by = c("suj", "sesage"))

plotdataset <- filter(plotdataset, variable == "rhythmlag" | variable == "burstlag")

plotdataset$sesage[plotdataset$sesage == 6]  <- "6-mo." 
plotdataset$sesage[plotdataset$sesage == 9]  <- "9-mo." 
plotdataset$sesage[plotdataset$sesage == 16] <- "16-mo." 

plotdataset$sesage <- factor(plotdataset$sesage, levels = c("6-mo.", "9-mo.", "16-mo."))

path2plot <- paste(tablepath, "individualcor_its_lcoh_", cohort, ".svg", sep ="")
plot_osc = ggplot(filter(plotdataset, variable == "burstlag"), aes(y=tau, x=value)) +
  geom_smooth(method = "lm", se = T, color = "black", alpha = .30, fill = "grey75") + ggpubr::theme_pubr() +
  geom_point(aes(color = "lightpink")) + 
  facet_grid(.~factor(sesage, levels = c("6-mo.", "9-mo.", "16-mo.")), scales = "free") + 
  theme(legend.position = "none") +
  settheme +
  xlab("Lagged Coh. Burst") + 
  ylab("Tau (s)") 

plot_freq = ggplot(filter(plotdataset, variable == "rhythmlag"), aes(y=tau, x=value)) +
  geom_smooth(method = "lm", se = T, color = "black", alpha = .30, fill = "grey75") + ggpubr::theme_pubr() +
  geom_point(aes(color = "lightpink")) + 
  facet_grid(.~factor(sesage, levels = c("6-mo.", "9-mo.", "16-mo.")), scales = "free") + 
  theme(legend.position = "none") +
  settheme +
  xlab("Lagged Coh. Rhythm") + 
  ylab("Tau (s)") 

ggpubr::ggarrange(plot_osc, plot_freq, nrow = 2)
ggsave(path2plot, width = 14.2, height = 11, units = "cm")


#### Space correlations between lagged coherence and tau  ----------------------

varnames <- as.vector(unique(pow$variable))
ses      <- c(6,9,16)


session       = NA
variable      = NA
rs            = NA
n             = NA
ci_min        = NA
ci_max        = NA
pval          = NA

loopidx = 1 
for (sesi in ses){
  
  data1 <- taudata_space%>%filter(sesage == sesi)%>%
    group_by(elect)%>%summarise(tau    = mean(tau, na.rm = TRUE), 
                                pernan = mean(pernan, na.rm = TRUE))
  
  for (vari in varnames) {
    
    data2  <- pow%>%filter(sesage == sesi & variable == vari)%>%
      group_by(elect)%>%summarise(voi = mean(value, na.rm = TRUE))
    
    cormat <-merge(data1,data2, by = "elect")
    
    cor_results <- pcor.test(cormat$tau, cormat$voi, cormat[c("pernan")], method = "spearman")
    
    boot_results <- boot(cormat, partial_corr_lag, R = 5000)
    ci           <- boot.ci(boot_results, type = "perc")
    
    session[loopidx]       = sesi
    variable[loopidx]      = vari
    rs[loopidx]            = cor_results$estimate
    n[loopidx]             = cor_results$n
    ci_min[loopidx]        = ci$percent[4]
    ci_max[loopidx]        = ci$percent[5]
    pval[loopidx]          = cor_results$p
    loopidx = loopidx + 1 
    
  }
}

session  <- c(session)
varaible <- c(variable)
rs       <- c(rs)
n        <- c(n)
ci_min   <- c(ci_min)
ci_max   <- c(ci_max)
pval     <- c(pval)

pval <- p.adjust(pval, method = "fdr")
results_correlation <- data.frame(n, session, variable, rs, pval, ci_min, ci_max, stringsAsFactors = TRUE)


results_correlation$significance[results_correlation$pval < .001] <- "***"
results_correlation$significance[results_correlation$pval < .01 & results_correlation$pval > .001] <- "**"
results_correlation$significance[results_correlation$pval < .05 & results_correlation$pval > .01]  <- "*"
results_correlation$significance[results_correlation$pval < .1 & results_correlation$pval > .05]  <- ""
results_correlation$significance[results_correlation$pval > .1]  <- ""

results_correlation$session[results_correlation$session==6]  <- "6-mo."
results_correlation$session[results_correlation$session==9]  <- "9-mo."
results_correlation$session[results_correlation$session==16] <- "16-mo."

results_correlation$session <- factor(results_correlation$session, levels = c("6-mo.", "9-mo.", "16-mo."))


results_table        <- dplyr::select(results_correlation, c("session", "n", "variable", "rs", "ci_min", "ci_max", "significance"))
results_table$rs     <- round(results_table$rs, 2)
results_table$ci_min <- round(results_table$ci_min,2)
results_table$ci_max <- round(results_table$ci_max,2)

results_table$stats  <- paste(results_table$rs, results_table$significance, " [", results_table$ci_min, " - ", results_table$ci_max, "]", sep ="")
results_table        <- dplyr::select(results_table, c("session", "variable", "n", "stats"))

results_table                  <- dcast(results_table, session + n ~ variable, value.var = "stats")

path2table = paste(tablepath, "tau_correlation_brythm_space_", cohort, ".docx", sep ="")
save_as_docx(nice_table(results_table),path = path2table)

#### Space correlation plots  --------------------------------------------------
plotdataset <- merge(pow, taudata_space, by = c("suj", "sex", "sesage", "elect", "area"))

plotdataset <- plotdataset%>%group_by(elect, sesage, area, variable)%>%
  summarise(value = mean(value, na.rm = TRUE),
            tau   = mean(tau, na.rm = TRUE))

plotdataset <- filter(plotdataset, variable == "rhythmlag" | variable == "burstlag")


plotdataset$sesage[plotdataset$sesage == 6]  <- "6-mo." 
plotdataset$sesage[plotdataset$sesage == 9]  <- "9-mo." 
plotdataset$sesage[plotdataset$sesage == 16] <- "16-mo." 

plotdataset$sesage <- factor(plotdataset$sesage, levels = c("6-mo.", "9-mo.", "16-mo."))

path2plot <- paste(tablepath, "spatial_stability_its_burstlag_", cohort, ".svg", sep ="")
plot_osc = ggplot(filter(plotdataset, variable == "burstlag"), aes(y=tau, x=value)) + geom_point(aes(color = area)) +
  geom_smooth(method = "lm", color = "black", se = T, alpha = .30, fill = "grey75") +
  scale_color_manual(values =c("#C4961A","#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")) + 
  theme_pubr() + facet_grid(.~factor(sesage, levels = c("6-mo.", "9-mo.", "16-mo.")), scales = "free") +
  settheme +
  labs(x = "Lagged Coh. Burst", 
       y = "Tau (s)", 
       color = "Area")

plot_freq = ggplot(filter(plotdataset, variable == "rhythmlag"), aes(y=tau, x=value)) + geom_point(aes(color = area)) + 
  geom_smooth(method = "lm", color = "black", se = T, alpha = .30, fill = "grey75") +
  scale_color_manual(values =c("#C4961A","#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")) + 
  theme_pubr() + facet_grid(.~factor(sesage, levels = c("6-mo.", "9-mo.", "16-mo.")), scales = "free") +
  settheme +
  labs(x = "Lagged Coh. Rhythm", 
       y = "Tau (s)", 
       color = "Area")

ggarrange(plot_osc, plot_freq, nrow = 2, common.legend = T, legend = "bottom")
ggsave(path2plot, width = 16, height = 12, units = "cm")
}
