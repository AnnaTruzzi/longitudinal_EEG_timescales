## Correlations with peak frequency and oscillatory power in alpha, and other bands 
library(ggplot2)
library(ggpubr)
library(readr)
library(tidyverse)
library(psych)
library(reshape2)
library(flextable)
library(boot)
library(rempsyc)
library(ppcor)

settheme =  theme(strip.text.x = element_text(size = 9.5, face = "bold"), strip.text.y = element_text(size = 9.5, face = "bold")) + 
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), title = element_text(size = 10), axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold")) 

taudata      <- read_delim("~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/develop_tau_validation_v3.csv", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

validationds <- unique(taudata$subject)
rm(taudata)

taudata <- read_delim("~/Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/develop_tau_v3.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)
mainids <- unique(taudata$subject)
rm(taudata)


powsets <- list.files(path = "Library/CloudStorage/Dropbox/ITS_Project/PSD_Data/",
                              pattern = ".*\\.csv$",
                              recursive = TRUE,
                              full.names = TRUE)

data <- lapply(powsets, read.csv, header = TRUE)%>%bind_rows()


data <- filter(data, suj %in% c(validationds, mainids))%>%
  mutate(Group = ifelse(suj %in% validationds, "Validation", "Exploratory"))%>%
  filter(rsq > .949)

fit_descriptives = data %>%
  group_by(suj, session, Group) %>%
  summarise(rsq = mean(rsq, na.rm =T))%>%
  ungroup() %>%
  group_by(session, Group) %>%
  summarise(mrsq = round(mean(rsq, na.rm =T),3),
            sdrsq = round(sd(rsq, na.rm =T),3),
            minrsq = round(min(rsq, na.rm =T),3),
            maxrsq = round(max(rsq, na.rm =T),3))

psd_plots = data %>%
  filter(rsq > .949) %>%
  group_by(session, Group, hz)%>%
  summarise(moscpow = mean(oscpow, na.rm = T),
            sdoscpow = sd(oscpow, na.rm = T)/sqrt(length(unique(suj))), 
            merror = mean(error, na.rm = T),
            sderror = sd(error, na.rm = T)/sqrt(length(unique(suj)))) 

oscpow_plot = ggplot(psd_plots, aes(x = hz, y = moscpow, colour = Group, fill = after_scale(color))) + 
  geom_ribbon(aes(ymin = moscpow - sdoscpow, ymax = moscpow + sdoscpow), alpha = 0.33) +geom_line(linewidth = 1.2) +
  facet_wrap(~session, scales = "free_x") + 
  scale_color_manual(values = c("Exploratory" = "steelblue", "Validation" = "coral")) +
  ggpubr::theme_pubr() + settheme + 
  labs(x = "Frequency (Hz)", y = "Oscillatory Power", color = "Group")
ggsave("Library/CloudStorage/Dropbox/ITS_Project/PSD_Data/oscpow_plot.jpg", oscpow_plot, width = 8, height = 4)
  

error_plot = ggplot(psd_plots, aes(x = hz, y = merror, color = Group, fill = after_scale(color))) + geom_line(linewidth = 1.2) +
  facet_wrap(~session, scales = "free_x") + 
  scale_color_manual(values = c("Exploratory" = "steelblue", "Validation" = "coral")) +
  ggpubr::theme_pubr() + settheme + geom_ribbon(aes(ymin = merror - sderror, ymax = merror + sderror), alpha = 0.33) +
  labs(x = "Frequency (Hz)", y = "Absolute Error", color = "Group")
ggsave("Library/CloudStorage/Dropbox/ITS_Project/PSD_Data/error_plot.jpg", error_plot, width = 8, height = 4)



