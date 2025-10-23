## Correlations with peak frequency and oscillatory power in alpha, and other bands 
library(readr)
library(ggrain)
library(ggplot2) 
library(tidyverse)
library(ggpubr)

tablepath = "~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/"
tauorign   <- read_csv("Desktop/Tau_data/tau_development_exploratory_v4.csv") #Load exploratory cohort data - this data has been already averaged and excluded pernans
tauvalid   <- read_csv("Desktop/Tau_data/tau_development_validation_v4.csv") #Load validation cohort data - this data has been already averaged and excluded pernans
tauadult   <- read_csv("Desktop/Tau_data/tau_development_adults_v4.csv") #Load adult cohort data

# Set the theme for the text
settheme =  theme(strip.text.x = element_text(size = 10, face = "bold"), strip.text.y = element_text(size = 10, face = "bold")) + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), title = element_text(size = 10), axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold")) 

dataorig    <- tauorign%>%group_by(subject, ses_age)%>%
  summarise(tau = mean(tauimp, na.rm = T)) # Compute the mean tau for each subject and session age in the exploratory cohort

dataval     <- tauvalid%>%group_by(subject, ses_age)%>%
  summarise(tau = mean(tauimp, na.rm = T)) # Compute the mean tau for each subject and session age in the validation cohort

datadult    <- tauadult%>%group_by(subject, event)%>%filter(pernan<.25)%>%
  summarise(tau = mean(tauimp, na.rm = T)) # Compute the mean tau for each subject and event in the adult cohort

tauvals <- c(dataorig$tau, dataval$tau, datadult$tau) # Combine all tau values from the exploratory, validation, and adult cohorts

miny    <- round(min(tauvals), 1) # Compute the maximum and minimum to have the same Y scale in the three groups
maxy    <- round(max(tauvals), 1) + .1

dataorig = dataorig%>%mutate(cohort = "Exploratory") # Add a cohort variable to the exploratory cohort
dataval  = dataval%>%mutate(cohort  = "Validation") # Add a cohort variable to the validation cohort
datadev  = rbind(dataorig, dataval) # Combine the exploratory and validation cohorts into a single dataset for development


plot_reg = ggplot(datadev, aes(x = ses_age, y = tau, color = factor(ses_age))) + 
  geom_jitter(alpha = .75, size = 0.75, width = .25) + 
  geom_line(aes(group = subject), alpha = .25, size = .25, color = "grey65") +
  geom_smooth(method = "lm", formula = y ~ poly(x,2)  ,aes(group = 1), color = "black", size = 1, linewidth = .5, alpha = .33) +
  scale_color_manual("Session", values=c("grey20", "#5C3030", "#BC7575")) + 
  xlab("Session (months)") + ylab("Tau (s)") + settheme + ggpubr::theme_pubr() + theme(legend.position = "none") +
  scale_x_continuous(breaks = c(6, 9, 16)) + facet_wrap(factor(datadev$cohort))
path2plot <- paste(tablepath, "tau_development_curve.svg")
ggsave(path2plot, plot = plot_reg, width = 10, height = 5, units = "cm") 


datadev$ses_age <- factor(datadev$ses_age, levels = c(6,9,16), labels = c("6", "9", "16")) #Order the ages
datadult$event   <- factor(datadult$event,    levels = c("Video", "EO", "EC"), labels = c("Video", "EO", "EC")) #Order adult-cohort condition and change the labels


#Create the plot for the developmental sample
path2plot  <- paste(tablepath, "tau_development.svg", sep ="")
baseplot <- ggplot(datadev, aes(x = ses_age, y = tau, fill = ses_age)) + geom_rain(boxplot.args = list(outlier.shape = NA, alpha = .75), violin.args = list(color = "black", alpha = .5), point.args = list(colour = "black", fill = "grey75", alpha = .5, size = .33), rain.side = 'l', id.long.var = "subject", line.args = list(alpha = .40, color = "grey90")) + theme_pubr()
aestplot <- baseplot +  scale_fill_manual("Session", values=c("grey20", "#5C3030", "#BC7575")) +  stat_summary(mapping = aes(x = factor(ses_age), y =  tau), fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}, fun = median, color = "red", size = .33, alpha = 1)
aestplot <- aestplot + stat_summary(fun = median, color = "red", geom = "line", group = 1, alpha = 1, linewidth = .5) + ylim(c(0,.4)) + facet_grid(.~cohort)
plotdev <- aestplot + settheme + xlab("Session (months)") + ylab("Tau (s)") + theme(legend.position = "none") 
ggsave(path2plot, plot = plotdev, width = 10, height = 5, units = "cm") 

# Create the plot for the adult sample
path2plot  <- paste(tablepath, "tau_adults.svg", sep ="")
baseplot <- ggplot(datadult, aes(x = event, y = tau, fill = event)) + geom_rain(boxplot.args = list(outlier.shape = NA, alpha = .75), violin.args = list(color = "black", alpha = .5), point.args = list(colour = "black", fill = "grey75", alpha = .5, size = .33), rain.side = 'l', id.long.var = "subject", line.args = list(alpha = .40, color = "grey90")) 
aestplot <- baseplot +  scale_fill_manual("Condition", values=c("beige", "beige", "beige")) +  stat_summary(mapping = aes(x = event, y =  tau), fun.min = function(z) {quantile(z,0.25)}, fun.max = function(z) {quantile(z,0.75)}, fun = median, color = "red", size = .33, alpha = 1) + ylim(c(0,.4)) + ggpubr::theme_pubr() + settheme
plot_tauadult  <- aestplot + ylab("Tau (s)") + xlab("Block") + labs(subtitle ='') + theme(legend.position = "none")
ggsave(path2plot, plot = plot_tauadult, width = 5, height = 5, units = "cm") 

# Combine both plots and save it
path2plot  <- paste(tablepath, "tau_comparison.svg", sep ="")
plotmerged <- ggarrange(plotdev, plot_tauadult, nrow = 1, widths = c(2, 1)) 
ggsave(path2plot, plot = plotmerged, width = 15, height = 6, units = "cm")


path2plot  <- paste(tablepath, "tau_comparison_curve.svg", sep ="")
plotmerged <- ggarrange(plot_reg, plot_tauadult, nrow = 1, widths = c(2, 1)) 
ggsave(path2plot, plot = plotmerged, width = 15, height = 6, units = "cm")

  


