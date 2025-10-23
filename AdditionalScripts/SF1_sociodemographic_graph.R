ses_theme <-   theme(axis.line = element_line(linewidth = 0, colour = "black"), strip.background = element_rect(color = "black", fill = NA, size =1)) + 
  theme(panel.margin=unit(.1, "lines"), panel.border = element_rect(color = "black", fill = NA, size = 1), strip.background = element_rect(color = "black", fill = NA, size = 1)) + 
  theme(strip.text.x = element_text(size = 10, face = "bold"), strip.text.y = element_text(size = 10, face = "bold")) + theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11), title = element_text(size = 11), axis.title.x = element_text(size = 11, face = "bold"), axis.title.y = element_text(size = 11, face = "bold")) + 
  theme(legend.position = "none")
ggplot2::theme_set(ses_theme)

validation      <- read_csv("Desktop/Tau_data/tau_development_validation_v4.csv")
validation_sujs <- unique(validation$subject)

exploratory      <- read_csv("Desktop/Tau_data/tau_development_exploratory_v4.csv")
exploratory_sujs <- unique(exploratory$subject)

bdi        <- readxl::read_excel("Desktop/Tau_data/bexat_data_quest.xlsx", 
                            sheet = "bdi")%>%
  mutate(group= ifelse(suj %in% validation_sujs, "Validation",
                        ifelse(suj %in% exploratory_sujs, "Exploratory", NA)))%>%
  filter(!is.na(group))%>%
  mutate(BDI = if_else(is.na(BDI_6m), 
                       if_else(is.na(BDI_16m), BDI_36m, BDI_16m), BDI_6m))

ses   <- readxl::read_excel("Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/bexat_data_quest.xlsx", 
                  sheet = "ses")

chaos <- readxl::read_excel("Library/CloudStorage/OneDrive-Personal/Documentos/Papers_msi/Thesis/Shared_Code_and_Data/bexat_data_quest.xlsx", 
                           sheet = "chaos ")

data <- merge(bdi, ses, by.x = "suj", by.y = "suj")
data <- merge(data, chaos, by.x = "suj", by.y = "suj")
           
bdi = ggplot(data, aes(x =group, y = BDI)) +  geom_jitter(width = .1) + theme_classic2()+
  geom_hline(yintercept = 14, color = 'red', size = 1, linetype = "dotted") + geom_hline(yintercept = 20, color = 'red', size = 1, linetype = "dotted") +
  geom_hline(yintercept = 29, color = 'red', size = 1, linetype = "dotted") +
  annotate("text",
           x = c(0.65, 0.65, 0.65),
           y = c(15.5, 21.5, 30.5),
           label = c("Mild", "Mod.", "Sev."),
           family = "", fontface = 2, size=4, colour = "red") + 
  ylab("BDI Score") + xlab("Cohort") + ggpubr::theme_pubr() + labs(subtitle = 'A')

income = ggplot(data, aes(x = IncNeed_6m, fill = group)) + 
  geom_histogram(color = 'black') + 
  scale_fill_manual(values = wesanderson::wes_palette("GrandBudapest1", n = 2)) + 
  theme(legen.position = "south") + 
  labs(fill = 'Cohort') + 
  geom_vline(xintercept = 1, color = 'red', size = 1, linetype = "dotted") +
  annotate("text",
           x = c(2),
           y = c(10),
           label = c("Poverty Line"),
           family = "", fontface = 2, size=4, colour = "red") + 
  ylab("Count") + xlab("Income to Needs Ratio") + 
  ggpubr::theme_pubr() + 
  theme(legend.position = "bottom") + theme(legend.title = element_blank()) + 
  labs(subtitle = 'B')

chaos <- ggplot(data, aes(x = CHAOS_6m, fill = group)) +  geom_histogram(color = 'black') + 
  ylab("Count") + xlab("CHAOS")    + ggpubr::theme_pubr() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = wesanderson::wes_palette("GrandBudapest1", n = 2)) + 
  labs(fill = 'Cohort') + 
  ggpubr::theme_pubr() + 
  theme(legend.position = "bottom") + theme(legend.title = element_blank()) + 
  labs(subtitle = 'C')


education <- ggplot(data, aes(x = EdMadre_6m, fill = group)) +  geom_histogram(color = 'black') + 
  ylab("Count") + xlab("Mother's Education")    + ggpubr::theme_pubr() + theme(legend.title = element_blank()) +
  scale_fill_manual(values = wesanderson::wes_palette("GrandBudapest1", n = 2)) + 
  labs(fill = 'Cohort') + 
  ggpubr::theme_pubr() + 
  theme(legend.position = "none") + theme(legend.title = element_blank()) + 
  labs(subtitle = 'D')

education_father <- ggplot(data, aes(x = EdPadre_6m, fill = group)) +  geom_histogram(color = 'black') + 
  ylab("Count") + xlab("Father's Education")    + ggpubr::theme_pubr()+ 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = wesanderson::wes_palette("GrandBudapest1", n = 2)) + 
  labs(fill = 'Cohort') + 
  ggpubr::theme_pubr() + 
  theme(legend.position = "none") + theme(legend.title = element_blank()) + 
  labs(subtitle = '')

education_merge <- ggpubr::ggarrange(education, education_father)

occupation <- ggplot(data, aes(x = LabMadre_6m, fill = group)) +  geom_histogram(color = 'black') + 
  ylab("Count") + xlab("Mother's Occupation")    + ggpubr::theme_pubr() + 
  theme(legend.position = "bottom") + theme(legend.title = element_blank()) +
  scale_fill_manual(values = wesanderson::wes_palette("GrandBudapest1", n = 2)) + 
  theme(legen.position = "south") + 
  labs(fill = 'Cohort') + 
  ggpubr::theme_pubr() + 
  theme(legend.position = "none") + theme(legend.title = element_blank()) + 
  labs(subtitle = 'E')

occupation_father <- ggplot(filter(data, !is.na(LabPadre_6m)), aes(x = LabPadre_6m, fill = group)) +  geom_histogram(color = 'black') + 
  ylab("Count") + xlab("Father's Occupation")    + ggpubr::theme_pubr() + 
  theme(legend.position = "bottom") + theme(legend.title = element_blank()) +
  scale_fill_manual(values = wesanderson::wes_palette("GrandBudapest1", n = 2)) + 
  labs(fill = 'Cohort') + 
  ggpubr::theme_pubr() + 
  theme(legend.position = "none") + theme(legend.title = element_blank()) + 
  labs(subtitle = '')
  
  occupation_merge <- ggpubr::ggarrange(occupation, occupation_father)
  
  top_plot    <- ggpubr::ggarrange(bdi, income, chaos, ncol = 3, common.legend = T)
  bottom_plot <- ggpubr::ggarrange(education_merge, occupation_merge, ncol = 2)
  full_plot   <- ggpubr::ggarrange(top_plot, bottom_plot, nrow = 2)
  ggsave("Desktop/Tau_data/ses_bdi_cohorts.svg", full_plot, width = 10, height = 6)
  
