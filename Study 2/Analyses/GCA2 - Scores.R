###############################################################################
# Gaze Contingent Avoidance Project
# Sabrina Gado & Yannik Stegmann
###############################################################################

library(readxl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(cowplot)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
path = getwd()

# Read data from Excel file
df_scores <- read_excel(file.path(path, "Questionnaires", "questionnaires.xlsx"))

# Recode variables from numeric to string
df_scores$gender <- recode(df_scores$gender, `1` = "male", `2` = "female", `3` = "diverse")
df_scores$handedness <- recode(df_scores$handedness, `1` = "right", `2` = "left")
df_scores$smoking <- recode(df_scores$smoking, `1` = "smoker", `2` = "no smoker")


# Calculate Questionnaire Scores
# --> Recoding of reversed items has already been done in sosci

# SPAI
df_spai <- df_scores %>%
  select(contains('spai'))
df_spai <- df_spai - 1 # Adapt scaling (from 1-7 to 0-6)

# Calculate means of nested items
columns_spai_multi <- names(df_spai)[grepl('_', names(df_spai))]
df_spai_multi <- df_spai %>%
  select(all_of(columns_spai_multi))

df_spai <- df_spai %>% select(-(columns_spai_multi))

items <- unique(as.integer(substring(names(df_spai_multi), 5, regexpr('_', names(df_spai_multi)) - 1)))

for (item in items) {
  df_spai_multi_subset <- df_spai_multi %>%
    select(contains(sprintf('spai%d_', item)))
  df_spai <- df_spai %>%
    mutate("{sprintf('spai%d', item)}" := rowMeans(df_spai_multi_subset, na.rm = TRUE))
}

df_spai <- df_spai %>%
  mutate(SPAI = rowMeans(df_spai, na.rm = TRUE)) %>% # Calculate mean of scale
  select(SPAI)

# SIAS
df_sias <- df_scores %>%
  select(contains('sias'))

df_sias <- df_sias - 1 # Adapt scaling (from 1-5 to 0-4)
df_sias <- df_sias %>% 
  mutate(SIAS = rowSums(.)) %>% # Calculate sum of items
  select(SIAS)

# STAI-T
df_stai <- df_scores %>%
  select(contains(c('stai')))

df_stai <- df_stai %>%
  mutate(STAI_T = rowSums(.)) %>% 
  select(STAI_T)

# Intolerance of Uncertainty
df_ui <- df_scores %>%
  select(contains('ui'))

df_ui <- df_ui %>%
  mutate(UI = rowSums(.)) %>% # Calculate sum of items
  select(UI)

# VAS: State Anxiety, Nervousness, Distress, Stress
df_vas <- df_scores %>%
  select(contains('VAS')) %>%
  mutate(
    VAS_diff_anxiety = VAS_end_anxiety - VAS_start_anxiety,
    VAS_diff_nervous = VAS_end_nervous - VAS_start_nervous,
    VAS_diff_distress = VAS_end_distress - VAS_start_distress,
    VAS_diff_stress = VAS_end_stress - VAS_start_stress
  )

# Debriefing / Control Questions
rep_str = c('-'='',
            '\r\n'=',',
            ';'=',',
            'ü' = 'ue',
            'ä' = 'ae',
            'ö' = 'oe',
            'ß' = 'ss')
df_debriefing <- df_scores %>%
  select(purpose, variables) %>%
  mutate(across('purpose', str_replace_all, rep_str)) %>% 
  mutate(across('variables', str_replace_all, rep_str))

# Labbook
df_labbook <- df_scores %>%
  select(labbook) %>% 
  mutate(across('labbook', str_replace_all, rep_str))

# Demographic Information
df_scores$VP <- as.integer(df_scores$ID)
df_demo <- df_scores %>%
  select(VP, gender, age, motivation, tiredness, handedness, motivation_points)

# Create Summary
df_summary <- bind_cols(df_demo, df_debriefing, df_labbook, df_vas, df_spai, df_sias, df_stai, df_ui)
print(paste("Age = ", round(mean(df_summary$age), 1), ", SD = ", round(sd(df_summary$age), 1), ", Range = ", round(min(df_summary$age)), " - ", round(max(df_summary$age))))
table(df_summary$gender)
prop.table(table(df_summary$gender)) * 100

table(df_summary$handedness)

# Write summary to CSV
write.csv2(df_summary, file.path(path, "Questionnaires", "demo_scores.csv"), row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8")

# # Write summary to CSV for jamovi
# df_summary <- df_summary %>% 
#   select(-c(purpose, variables, labbook))
# write.csv2(df_summary, file.path(path, "Questionnaires", "demo_scores_jamovi.csv"), row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8")

# Histogram
plot_spai <- ggplot(df_summary, aes(x = SPAI)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.2, color="grey", size=0.3) +
  geom_vline(aes(xintercept = median(SPAI), color = "Median"), size=1) +
  geom_vline(aes(xintercept = 2.79, color = "Cut-Off (Remission)"), size=1) +
  geom_density(size = 0.5, color="darkgrey") +
  scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, 1)) +
  labs(x = "SPAI (Social Anxiety)", y = "Density") +
  theme_minimal() +
  scale_colour_manual(values = c("Cut-Off (Remission)" ='yellowgreen', "Median" ='deepskyblue4'), name="")
ggsave(file.path(path, "Plots", "Distribution_SPAI.png"),  width=1800, height=1000, units="px")


plot_sias <- ggplot(df_summary, aes(x = SIAS)) +
  geom_histogram(aes(y = ..density..), binwdith=2, color="grey", size=0.3) +
  geom_vline(aes(xintercept = median(SIAS), color = "Median"), size=1) +
  geom_vline(aes(xintercept = 30, color = "Clinical Cut-Off"), size=1) +
  geom_density(size = 0.5, color="darkgrey") +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 65, 5)) +
  labs(x = "SIAS (Social Anxiety)", y = "Density") +
  theme_minimal() +
  scale_colour_manual(values = c("Clinical Cut-Off" ='yellowgreen', "Median" ='deepskyblue4'), name="")
ggsave(file.path(path, "Plots", "Distribution_SIAS.png"),  width=1800, height=1000, units="px")


plot_stai <- ggplot(df_summary, aes(x = STAI_T)) +
  geom_histogram(aes(y = ..density..), color="grey", size=0.3) +
  geom_vline(aes(xintercept = median(STAI_T), color = "Median"), size=1) +
  geom_vline(aes(xintercept = 39, color = "Clinical Cut-Off"), size=1) +
  geom_density(size = 0.5, color="darkgrey") +
  scale_x_continuous(limits = c(20, 70), breaks = seq(20, 70, 5)) +
  labs(x = "STAI (Trait Anxiety)", y = "Density") +
  theme_minimal() +
  scale_colour_manual(values = c("Clinical Cut-Off" ='yellowgreen', "Median" ='deepskyblue4'), name="")
ggsave(file.path(path, "Plots", "Distribution_STAI.png"),  width=1800, height=1000, units="px")

score_plot <- plot_grid(plot_spai, plot_stai, ncol = 2, labels=c("C", "D"), align="vh")
ggsave(file.path(path, "Plots", "scores.png"),  width=3500, height=1000, units="px")