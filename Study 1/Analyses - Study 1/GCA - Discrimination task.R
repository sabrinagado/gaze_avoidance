###############################################################################
# Gaze Contingent Avoidance Project
# Sabrina Gado & Yannik Stegmann
###############################################################################

# Read Packages
library(tidyverse)
library(cowplot)
library(viridis)
library(ez)
library(apa)


# Functions
se <- function(x, na.rm = TRUE) {
  sd(x, na.rm) / sqrt(if(!na.rm) length(x) else sum(!is.na(x)))
}


# global variables
theme_set(theme_minimal(base_size = 16))
vps_summary = tibble()
code_times = c()

# read experiment files
path.behavior = file.path("Study 1", "Experiment", "gaze_discrimination_task", "data")

filemat = list.files(path.behavior, pattern = "*.csv") # read all files from data-folder into a single list


# cycle through all subject csv files
# Initialise a progress bar
cat("\n", "Read in subject csv files for discrimination performance:", "\n")
pb <- txtProgressBar(min = 1, max = length(filemat), style = 3)
i = 1
for (subject in filemat){
  # subject <- filemat[2]  #use  file for practice and coding changes
  
  vp_data <- read.csv(file.path(path.behavior, subject)) # read data from files per subject
  
  participant <- vp_data$participant %>% unique()
  
  vp_summary <- vp_data %>% 
    filter(str_detect(trialtype,"cs_")) %>% 
    mutate(condition = ifelse(trialtype == comptype, "same", "different"),
           social = ifelse(trialtype %in% c("cs_plus_s","cs_minus_s"), "social","non-social")) %>%
    select(participant, condition, social, correctness, key.rt) %>%
    summarise(correct_sum = sum(correctness == "correct"),
              correct_mean = sum(correctness == "correct")/n(),
              rt_mean = mean(key.rt,na.rm = T),
              .by = c(participant, social, condition))
  
  # Create summary data
  vps_summary <- rbind(vps_summary, vp_summary)
  
  ggplot(vp_summary, aes(x = social, fill = condition, y = correct_mean)) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "red") + 
    geom_col(position = position_dodge()) +
    scale_y_continuous("% correct", expand = c(0,0), limits = c(0,1.05)) +
    scale_x_discrete("Category") +
    scale_fill_viridis_d("Condition",end = 0.75)
    # theme_classic()
  ggsave(file.path("Study 1", "Plots", "Discrimination", "individuals", paste0(participant, ".png")), width=6, height=5)
  
  setTxtProgressBar(pb, i)
  i = i+1
}

close(pb)

vps_summary <- vps_summary %>% 
  ungroup

# Accuracy
vps_summary_grouped <- vps_summary %>%
  summarise(correct_mean = mean(correct_mean), rt_mean = mean(rt_mean) * 1000, .by=c(participant, social))

print(vps_summary_grouped %>% summarise(correct_sd = sd(correct_mean) * 100, correct_mean = mean(correct_mean) * 100, .by=c(social)))

cat("\n\n", "T-test of Discrimination task", "\n")
vps_summary_grouped.test <- vps_summary_grouped %>%
  mutate(social = str_remove_all(social, "-")) %>% 
  pivot_wider(id_cols = participant, names_from = social, values_from = correct_mean)
apa::t_test(vps_summary_grouped.test$social, vps_summary_grouped.test$nonsocial, alternative = "two.sided", paired=TRUE) %>% apa::t_apa()

vps_summary_grouped %>%   
  group_by(social) %>%
  summarise(correct = mean(correct_mean, na.rm = T), se = se(correct_mean)) %>% 
  ggplot(aes(x = social, y = correct, fill = social)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "black") + 
  geom_col(position = position_dodge(), alpha=0.5) +
  geom_errorbar(aes(ymax = correct + se, ymin = correct - se), width = 0.4) +
  geom_point(data = vps_summary_grouped, aes(x = social, y = correct_mean), size = 2, shape = 21, color = "black", alpha=0.5, position = position_jitter(width=0.2, height=0.005)) +
  # geom_beeswarm(data = vps_summary, aes(x = social, y = correct_mean), alpha = 0.5, mapping = aes(x = social)) + 
  labs(x = "Category", y = "% correct") + # title = paste("Proportion of Correct Discrimination (N = ", n_distinct(vps_summary_grouped$participant), ")", sep=""), 
  # scale_y_continuous("% correct", expand = c(0,0), limits = c(0,1.05)) +
  scale_x_discrete(NULL, labels = c("Non-Social","Social")) +
  # scale_fill_viridis_d("Condition", end = 0.25, begin = 0.25, guide = "none" ) + 
  # theme_minimal() + 
  # scale_color_viridis_d() +
  theme(legend.position = "none")
ggsave(file.path("Study 1", "Plots", "Discrimination", "ga_correct.png"), width=1600, height=1600, units="px")


vps.discrimination.score <- vps_summary %>%
  summarise(discrimination = mean(correct_mean), .by=c("participant", "social"))
names(vps.discrimination.score) = c("subject","condition_social","discrimination")
write.csv2(vps.discrimination.score, file.path("Study 1", "Behavior", "discrimination.csv"), row.names=FALSE, quote=FALSE)
