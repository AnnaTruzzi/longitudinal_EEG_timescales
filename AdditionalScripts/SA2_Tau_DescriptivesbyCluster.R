# Supplementary plots TAU per region
library(readr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggrain)
library(rempsyc)

# Set the plot theme
settheme =  theme(strip.text.x = element_text(size = 10, face = "bold"), strip.text.y = element_text(size = 10, face = "bold")) + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), title = element_text(size = 10), axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold")) 

tablepath = "~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/" # General path for the data
cohorts   = c('exploratory', 'validation', 'adults')

# Load the demographic information 
desc    <- read_delim("~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/basicdemos.csv",
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)%>%filter(sesage < 36) # Exclude the 36mo information

for (cohort in cohorts){ # Loop for each developmental cohort 
  
if (cohort != "adults") {
  if (cohort == "exploratory") {
    taudata <- read_csv("Desktop/Tau_data/tau_development_exploratory_v4.csv")%>%
      rename(sex = Gender, area = Area, channel_number = channel, channel = Label)
    
  } else {
    taudata <- read_csv("Desktop/Tau_data/tau_development_validation_v4.csv")%>%
      rename(sex = Gender, area = Area, channel_number = channel, channel = Label)}

## Tau Nans --------------------------------------------------------------------
taudata$DummyNan[is.na(taudata$tauorg)] <- 1 # Set one to non-convergent epochs
taudata$DummyNan[!is.na(taudata$tauorg)] <- 0

data <- taudata%>%filter(pernan < .25, event != 9999)%>%group_by(ses_age, sex, area)%>%
  summarise(nanper = mean(DummyNan, na.rm = T), 
            nansd  = sd(DummyNan, na.rm = T)) # Create the percentage of non convergent epochs and its standard deviation

data       <- melt(data, id.vars = c("ses_age", "sex", "area"))
data$value <- round(data$value, 2) # Keep this for event and area

data      <- dcast(data, ses_age + sex + area ~ variable)
data$stat <- paste(data$nanper, " (", data$nansd, ")", sep = "")
data      <- dplyr::select(data, -nanper, -nansd) # Create a variable as M (SD) and remove the mean and sd 

desctable <- dcast(data, ses_age + sex ~ area)

desctable <- nice_table(desctable, separate.header = F)

path2table = paste(tablepath, "pernan_its_dev_", cohort, ".docx", sep="") # Name to store the descriptives
save_as_docx(desctable,path = path2table) # Save the descriptive table


data <- taudata%>%filter(pernan < .25, event != 9999)%>%group_by(subject, sex, ses_age, channel, area)%>%
  summarise(tau = mean(tauimp, na.rm = T)) # Keep only the events with good convergence and that are not artifacted epochs, average tau per channel keeping the information of the area

## Demographic info ------------------------------------------------------------
datadesc                              <- taudata%>%group_by(subject, sex, ses_age)%>%summarise(epochs = length(unique(event_n)))
datadesc$epochs[datadesc$epochs == 1] <- NA 
colnames(datadesc)                    <- c("suj", "sex", "sesage", "epochs")

descset                               <- merge(desc, datadesc, by = c("suj",  "sesage"), all.x = F, all.y = T)

descset <- dplyr::select(descset, -sex.x) # Remove sex.x
descset <- melt(descset, id.vars = c("suj", "sesage", "sex.y"))

descset <- descset%>%group_by(sesage, sex.y, variable)%>% # For each variable of interest compute the mean and sd (age, birth weigth, gestation weeks, and sex)
  summarise(m       = mean(value, na.rm = T), 
            sd      = sd(value, na.rm = T))

descset$m  <- round(descset$m, 2) # Round to the second decimal point
descset$sd <- round(descset$sd, 2)

descset$stats <- paste(descset$m, " (", descset$sd, ")", sep = "") # Create a variable as M (SD)
descset       <- dplyr::select(descset, -m, -sd) # Remove the mean and sd variables

colnames(descset) <- c('sesage', 'sex', 'descvar', 'stats')
desctable     <- reshape2::dcast(descset, sesage + sex ~ descvar, value.var = 'stats')

desctable <- nice_table(desctable, separate.header = F)

path2table = paste(tablepath, "demographics_", cohort, ".docx", sep="")
save_as_docx(desctable,path = path2table)

## Tau descriptives and plots --------------------------------------------------
### Avg. Tau --------------------------------------------------------------------
colnames(data) <- c("suj", "sex", "sesage", "channel", "area", "tau")

datatable <- data%>%group_by(sex, sesage)%>%
  summarise(n = length(unique(suj[tau>0])),
    mtau   = mean(tau, na.rm = T), 
            sdtau  = sd(tau, na.rm = T))%>%
  melt(id.vars = c("sesage", "sex"))

datatable$value <- round(datatable$value,2)
datatable       <- dcast(datatable, sesage + sex ~ variable)
datatable$tau <- paste(datatable$mtau, " (", datatable$sdtau, ")", sep = "")
desctable       <- dplyr::select(datatable, -mtau, -sdtau)

desctable <- nice_table(desctable, separate.header = F)

path2table = paste(tablepath, "desctiptives_itsdev_avg_", cohort, ".docx", sep="")
save_as_docx(desctable,path = path2table)

### Develop by Area ------------------------------------------------------------
datatable <- data%>%group_by(sex, sesage, area)%>%
  summarise(mtau   = mean(tau, na.rm = T), 
            sdtau  = sd(tau, na.rm = T))%>%
  melt(id.vars = c("sesage", "sex", "area"))

datatable$value <- round(datatable$value,2)
datatable       <- dcast(datatable, sesage + sex + area ~ variable)
datatable$stats <- paste(datatable$mtau, " (", datatable$sdtau, ")", sep = "")
desctable       <- dplyr::select(datatable, -mtau, -sdtau)
desctable       <- dcast(desctable, sesage + sex ~ area)

desctable <- nice_table(desctable, separate.header = F)

path2table = paste(tablepath, "desctiptives_itsdev_area_", cohort, ".docx", sep="")
save_as_docx(desctable,path = path2table)

taudata_plot <- data%>%group_by(suj, sesage, area)%>%
  summarise(tau = mean(tau, na.rm = TRUE))%>%ungroup()

plot_reg = ggplot(taudata_plot, aes(x = sesage, y = tau, color = factor(sesage))) + 
  geom_jitter(alpha = .75, size = 0.75, width = .25) + 
  geom_line(aes(group = suj), alpha = .25, size = .25, color = "grey65") +
  geom_smooth(method = "lm", formula = y ~ poly(x,2)  ,aes(group = 1), color = "black", size = 1, linewidth = .5, alpha = .33) +
  scale_color_manual("Session", values=c("grey20", "#5C3030", "#BC7575")) + 
  xlab("Session (months)") + ylab("Tau (s)") + settheme + ggpubr::theme_pubr() + theme(legend.position = "none") +
  scale_x_continuous(breaks = c(6, 9, 16)) + facet_wrap(factor(taudata_plot$area))
path2plot <- paste(tablepath, "taudevelop_area_curve_", cohort, ".svg")
ggsave(path2plot, width = 18, height = 10, units = "cm")


taudata_plot$sesage <- factor(taudata_plot$sesage, levels = c(6, 9, 16))
path2plot <- paste(tablepath, "taudevelop_area_", cohort, ".svg")
baseplot <- ggplot(taudata_plot, aes(x = sesage, y = tau, fill = sesage)) + 
  geom_rain(boxplot.args = list(outlier.shape = NA, alpha = .75), violin.args = list(color = "black", alpha = .5), point.args = list(colour = "black", fill = "grey75", alpha = .5, size = .33), rain.side = 'l', id.long.var = "suj", line.args = list(alpha = .40, color = "grey90")) + facet_wrap(factor(taudata_plot$area))
aestplot <- baseplot +  scale_fill_manual("Session", values=c("grey20", "#5C3030", "#BC7575")) +  stat_summary(mapping = aes(x = factor(sesage), y =  tau), fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}, fun = median, color = "red", size = .33, alpha = 1)
aestplot <- aestplot + stat_summary(fun = median, color = "red", geom = "line", group = 1, alpha = 1, linewidth = .5)
aestplot  + xlab("Session (months)") + ylab("Tau (s)") + settheme + ggpubr::theme_pubr() + theme(legend.position = "none")
ggsave(path2plot, width = 18, height = 10, units = "cm")


} else { 
  taudata <- read_csv("Desktop/Tau_data/tau_development_adults_v4.csv")%>%
    rename(area = Area, channel_number = channel, channel = Label)%>%
    mutate(event = if_else(event == "None", "Video", 
                               if_else(event == "eyes_open", "EO", 
                                       if_else(event == "eyes_closed", "EC", event))))%>%
    mutate(event = factor(event, levels = c("Video", "EO", "EC")))
  
  taudata$DummyNan[is.na(taudata$tauorg)] <- 1 # Set one to non-convergent epochs
  taudata$DummyNan[!is.na(taudata$tauorg)] <- 0
  
  data <- taudata%>%filter(pernan < .25, event != 9999)%>%group_by(event, area)%>%
    summarise(nanper = mean(DummyNan, na.rm = T), 
              nansd  = sd(DummyNan, na.rm = T)) # Create the percentage of non convergent epochs and its standard deviation
  
  data       <- melt(data, id.vars = c("event", "area"))
  data$value <- round(data$value, 2) # Keep this for event and area
  
  data      <- dcast(data, event + area ~ variable)
  data$stat <- paste(data$nanper, " (", data$nansd, ")", sep = "")
  data      <- dplyr::select(data, -nanper, -nansd) # Create a variable as M (SD) and remove the mean and sd 
  
  desctable <- dcast(data, event ~ area)
  
  desctable <- nice_table(desctable, separate.header = F)
  
  path2table = paste(tablepath, "pernan_its_dev_", cohort, ".docx", sep="") # Name to store the descriptives
  save_as_docx(desctable,path = path2table) # Save the descriptive table
  
  data <- taudata%>%filter(pernan < .25, event != 9999)%>%group_by(subject, event, channel, area)%>%
    summarise(tau = mean(tauimp, na.rm = T)) # Keep only the events with good convergence and that are not artifacted epochs, average tau per channel keeping the information of the area
  
   # Average tau
  colnames(data) <- c("suj", "event", "channel", "area", "tau")
  
  datatable <- data%>%group_by(event)%>%
    summarise(n = length(unique(suj[tau>0])),
              mtau   = mean(tau, na.rm = T), 
              sdtau  = sd(tau, na.rm = T))%>%
    melt(id.vars = c("event"))
  
  datatable$value <- round(datatable$value,2)
  datatable       <- dcast(datatable, event  ~ variable)
  datatable$tau <- paste(datatable$mtau, " (", datatable$sdtau, ")", sep = "")
  desctable       <- dplyr::select(datatable, -mtau, -sdtau)
  
  desctable <- nice_table(desctable, separate.header = F)
  
  path2table = paste(tablepath, "desctiptives_itsdev_avg_", cohort, ".docx", sep="")
  save_as_docx(desctable,path = path2table)
  
  # Average tau per area 
  datatable <- data%>%group_by(event, area)%>%
    summarise(mtau   = mean(tau, na.rm = T), 
              sdtau  = sd(tau, na.rm = T))%>%
    melt(id.vars = c("event", "area"))
  
  datatable$value <- round(datatable$value,2)
  datatable       <- dcast(datatable, event + area ~ variable)
  datatable$stats <- paste(datatable$mtau, " (", datatable$sdtau, ")", sep = "")
  desctable       <- dplyr::select(datatable, -mtau, -sdtau)
  
  desctable       <- dcast(desctable, event ~ area)
  desctable <- nice_table(desctable, separate.header = F)
  path2table = paste(tablepath, "desctiptives_itsdev_area_", cohort, ".docx", sep="")
  save_as_docx(desctable,path = path2table)
  
}
}


