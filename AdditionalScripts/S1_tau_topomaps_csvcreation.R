## TAU - Topomaps 

## Load the files --------------------------------------------------------------
tablepath = "Desktop/Tau_data/"

cohorts = c('exploratory', 'validation', 'adults')

for (cohort in cohorts) {
  
  if (cohort != "adults") {
    
    if (cohort == 'exploratory') {
      taudata = readr::read_csv("Desktop/Tau_data/tau_development_exploratory_v4.csv")
    } else if (cohort == 'validation') {
      taudata = readr::read_csv("Desktop/Tau_data/tau_development_validation_v4.csv") } 
    
    tau = taudata %>%
      group_by(ses_age, channel, Label) %>%
      summarise(t = mean(tauimp, na.rm = TRUE))
    write.csv(tau, paste0(tablepath, "tau_", cohort, "_v4.csv"))
    rm(tau)
    
    
    taudata$dummynan[is.na(taudata$tauorg)    & taudata$pernan < .25 & taudata$event != 9999]  <- 1 
    taudata$dummynan[is.na(taudata$tauorg)==0 & taudata$pernan < .25 & taudata$event != 9999]  <- 0
    
    nan_avg      <- taudata%>%group_by(ses_age, channel, Label)%>%
      summarise(t = mean(dummynan, na.rm = T))
    
    write.csv(nan_avg, paste0(tablepath, "tau_", cohort, "_v4_nanper.csv")) 
  } else { 
    taudata = readr::read_csv("Desktop/Tau_data/tau_development_adults_v4.csv")
    
    taudata = taudata%>%mutate(event = if_else(event == 'None', 'Video', 
                                               if_else(event == 'eyes_open', 'EO', 
                                                 if_else(event == 'eyes_closed', 'EC', as.character(event)))))
    
    tau = taudata %>%
      group_by(event, channel, Label) %>%
      summarise(t = mean(tauimp, na.rm = TRUE))
    
    write.csv(tau, paste0(tablepath, "tau_", cohort, "_v4.csv"))
    rm(tau)
    
    taudata$dummynan[is.na(taudata$tauorg)    & taudata$pernan < .25 & taudata$event != 9999]  <- 1
    taudata$dummynan[is.na(taudata$tauorg)==0 & taudata$pernan < .25 & taudata$event != 9999]  <- 0
    
    nan_avg      <- taudata%>%group_by(event, channel, Label)%>%
      summarise(t = mean(dummynan, na.rm = T))
    
    write.csv(nan_avg, paste0(tablepath, "tau_", cohort, "_v4_nanper.csv"))
  }
}

