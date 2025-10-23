library(tidyverse)
library(readr)
library(stringr)
library(reshape2)

elec             = readxl::read_excel("Desktop/elec_info_bexat.xlsx")
path2data        = "Desktop/"
path2save        = "Desktop/Tau_data/"
basefile         = "EEGdata_tau_"
sufixfile        = "_10s_numpy1000_longformat_seconds_nooutliers_lowmovement_nohighnans"
bexat_desc       = readr::read_csv("Desktop/tau_descriptives.csv")
exact_ages       = bexat_desc%>%dplyr::select(subject, age6m, age9m, age16m)
exact_ages       = reshape2::melt(exact_ages, id.var = 'subject')%>%
  mutate(ses_age = as.numeric(stringr::str_extract(variable, '\\d+')))%>%
           dplyr::select(-variable)


cohorts = c('exploratory', 'validation', 'adults')
version = 4

if (!dir.exists(path2save)) {
  dir.create(path2save)
  message(paste("Directory created:", path2save))
} 

for (cohort in cohorts) {
  if (cohort == 'adults') {
    ages = c('adults')
    valsufix = '_envelope.csv'
  } else { 
    ages    = c('6mo', '9mo', '16mo')
    if (cohort == 'validation') {
      valsufix = '_validation_envelope.csv'
    } else { 
      valsufix = '_envelope.csv'
    }
  }
for (age in ages) {
  filename         = paste0(path2data, basefile, age, sufixfile, valsufix)
  tau              = read_csv(filename)%>%
    mutate(cond_name= if_else((cond_name == 'None' | cond_name == 'video'), 'Video', 
                               if_else(cond_name == 'eyes_open', 'EO', 
                                 if_else(cond_name == 'eyes_closed', 'EC', as.character(cond_name)))))
  colnames(tau)[1] = 'RowNumber'
  
  tau              = tau%>%arrange(sub, channel, RowNumber)%>%
    mutate(RowNumber = RowNumber + 1)
  
  for (s in unique(tau$sub)) {
    
    subtau = filter(tau, sub == s)%>%mutate(event_n = 0)
    subrow = filter(subtau, channel == 0)
    
    if (length(subrow$channel) == 1) {
      subtau$event_n = 1
    } else { 
      l = 1 
      e = 1
      for (r in subtau$RowNumber) { 
        if (r == subtau$RowNumber[1] || subtau$channel[l] != subtau$channel[l-1]) { 
          subtau$event_n[l] = 1 
          l = l + 1
          e = 2
        } else {
          subtau$event_n[l] = e 
          e = e + 1 
          l = l + 1
        }
      } 
    }
    
    if (s == unique(tau$sub)[1]) {
      ses_tau = subtau
    } else { 
      ses_tau = rbind(ses_tau, subtau)}
    rm(subtau, subrow)
  }
  rm(tau)
  ses_tau          = merge(ses_tau, elec, by.x = 'channel', by.y = 'Number')%>%
    mutate(pernan  = sum(is.na(tau))/n(), .by = sub)%>%
    mutate(tauorg  = tau)%>%dplyr::select(-tau)%>%
    mutate(tauimp  = if_else(!is.na(tauorg), tauorg, 
                             mean(tauorg, na.rm = T)), .by = c(sub, event_n, Subarea))%>%
    mutate(tauimp  = if_else(!is.na(tauimp), tauimp, mean(tauorg, na.rm = T)), .by = c(sub, event_n, Area))%>%
    mutate(event   = cond_name, 
           subject = sub)%>%dplyr::select(-cond_name, -sub)%>%arrange(subject, channel, RowNumber) 
  
  if (cohort != 'adults') {
    ses_tau     = mutate(ses_tau, ses_age = as.numeric(stringr::str_extract(age, '\\d+')), 
                         age = (ses_age-6)/12, 
                         age_sq = age^2)
    if (age == ages[1]) { 
      all_tau = ses_tau
    } else {
      all_tau = rbind(all_tau, ses_tau)
      }
    rm(ses_tau)
  } 
}
  
  if (cohort != 'adults') { 
    descs = dplyr::select(bexat_desc, subject, Gender, Pesogramos, Semas)%>%
      rename(bweigth = Pesogramos, gweeks = Semas)
    all_tau = merge(all_tau, descs, by = 'subject')
    all_tau = merge(all_tau, exact_ages, by = c('subject', 'ses_age'))
  } else { 
    all_tau = ses_tau
    rm(ses_tau)}
  
  merged_filename = paste0(path2save, 'tau_development_', cohort, '_v', version, '.csv')
  write_csv(all_tau, merged_filename)
  rm(all_tau)
}




