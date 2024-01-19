###############################################################################
# Gaze Contingent Avoidance Project
# Sabrina Gado & Yannik Stegmann
###############################################################################

library(readxl)
library(dplyr)
library(tidyr)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
path = getwd()

# Read data from Excel file
df_scores <- read_excel(file.path(path, "questionnaires.xlsx"))

# Recode variables from numeric to string
df_scores$gender <- recode(df_scores$gender, `1` = "male", `2` = "female", `3` = "diverse")
df_scores$handedness <- recode(df_scores$handedness, `1` = "right", `2` = "left")
df_scores$smoking <- recode(df_scores$smoking, `1` = "smoker", `2` = "no smoker")

# SPAI
df_spai <- df_scores %>%
  select(contains('spai'))

# Adapt scaling (from 1-7 to 0-6)
df_spai <- df_spai - 1

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

# Calculate mean of scale
df_spai <- df_spai %>%
  mutate(SPAI = rowMeans(df_spai, na.rm = TRUE)) %>%
  select(SPAI)

# SIAS
df_sias <- df_scores %>%
  select(contains('sias'))

# Adapt scaling (from 1-5 to 0-4)
df_sias <- df_sias - 1

# Calculate sum of items
df_sias <- df_sias %>%
  mutate(SIAS = rowSums(.)) %>%
  select(SIAS)

# STAI-T
df_stai <- df_scores %>%
  select(contains('stai'))

# Calculate sum of respective items
df_stai <- df_stai %>%
  mutate(STAI_T = rowSums(.)) %>%
  select(STAI_T)

# Intolerance of Uncertainty
df_ui <- df_scores %>%
  select(contains('ui'))

# Calculate sum of respective items
df_ui <- df_ui %>%
  mutate(UI = rowSums(.)) %>%
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
  select(labbook, digitimer, temperature, humidity) %>% 
  mutate(across('labbook', str_replace_all, rep_str))

# Demographic Information
df_scores$VP <- as.integer(df_scores$ID)
df_demo <- df_scores %>%
  select(VP, gender, age, motivation, tiredness, handedness)

# Create Summary
df_summary <- bind_cols(df_demo, df_debriefing, df_labbook, df_vas, df_spai, df_sias, df_stai, df_ui)
print(paste("Age = ", round(mean(df_summary$age), 1), ", SD = ", round(sd(df_summary$age), 1), ", Range = ", round(min(df_summary$age)), " - ", round(max(df_summary$age))))
prop.table(table(df_summary$gender)) * 100

# Write summary to CSV
write.csv2(df_summary, file.path(path, "scores_summary.csv"), row.names=FALSE, quote=FALSE, fileEncoding = "UTF-8")

