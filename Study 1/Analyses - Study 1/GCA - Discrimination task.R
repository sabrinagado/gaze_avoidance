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


# global variables
vps_summary = tibble()
code_times = c()

# read experiment files
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
path = getwd()

path.behavior = file.path(path, "Experiment", "gaze_discrimination_task", "data")

filemat = list.files(path.behavior, pattern = "*.csv") # read all files from data-folder into a single list


#cycle through all subject csv files
for (subject in filemat){
  # subject <- filemat[2]  #use  file for practice and coding changes
   
  start.time <- Sys.time()
  
  vp_data <- read.csv(file.path(path.behavior, subject)) # read data from files per subject
  
  participant <- vp_data$participant %>% unique()
  
  vp_summary <- vp_data %>% 
    filter(str_detect(trialtype,"cs_")) %>% 
    mutate(condition = ifelse(trialtype == comptype, "same", "different"),
           social = ifelse(trialtype %in% c("cs_plus_s","cs_minus_s"), "social","non-social")) %>%
    select(participant, condition, social, correctness, key.rt) %>%
    group_by(participant,social,condition) %>%
    summarise(correct_sum = sum(correctness == "correct"),
              correct_mean = sum(correctness == "correct")/n(),
              rt_mean = mean(key.rt,na.rm = T))
  
  # Create summary data
  vps_summary <- rbind(vps_summary,vp_summary)
  
  ggplot(vp_summary, aes(x = social, fill = condition, y = correct_mean)) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "red") + 
    geom_col(position = position_dodge()) +
    scale_y_continuous("% correct", expand = c(0,0), limits = c(0,1.05)) +
    scale_x_discrete("Category") +
    scale_fill_viridis_d("Condition",end = 0.75) + 
    theme_classic()
  ggsave(file.path(path, "Plots", "Discrimination", "individuals", paste0(participant, ".png")))
    
  
  # stuff for calculating computing duration
  end.time <- Sys.time()  
  time.taken <- end.time - start.time
  code_times = c(code_times,round(as.numeric(time.taken),1))
  mean_code_time = mean(code_times)
  
  
  if(match(subject,filemat) == 1) {
    print(paste0(round(match(subject,filemat)/length(filemat)*100,1), "% ..... ",
                 match(subject,filemat), " of ", length(filemat), " files processed!"))
  } else {
    print(paste0(round(match(subject,filemat)/length(filemat)*100,1), "% ..... ",
                 match(subject,filemat), " of ", length(filemat), " files processed! ..... ", 
                 round((mean_code_time*length(filemat)-mean_code_time * match(subject,filemat))/60,1), " min remaining"))
  }
} #end-loop

vps_summary <- vps_summary %>% 
  ungroup

# Accuracy
vps_summary_grouped <- vps_summary %>%
  summarise(correct_mean = mean(correct_mean), rt_mean = mean(rt_mean) * 1000, .by=c(participant, social))

vps_summary_grouped %>%
  summarise(correct_sd = sd(correct_mean) * 100, correct_mean = mean(correct_mean) * 100, .by=c(social))

vps_summary_grouped %>%
  t_test(data=., correct_mean ~ social, paired = T, alternative = "two.sided") %>% 
  apa::t_apa()

vps_summary_grouped %>% 
  group_by(social) %>%
  summarise(correct = mean(correct_mean, na.rm = T), se = sd(correct_mean)/sqrt(n())) %>% 
  ggplot(aes(x = social, y = correct, fill = social)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "red") + 
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymax = correct + se, ymin = correct - se), width = 0.4) +
  geom_point(data = vps_summary_grouped, aes(x = social, y = correct_mean), size = 2, shape = 21, color = "black", alpha=0.3, position = position_jitter(width=0.2, height=0.005)) +
  # geom_beeswarm(data = vps_summary, aes(x = social, y = correct_mean), alpha = 0.5, mapping = aes(x = social)) + 
  labs(title = paste("Proportion of Correct Discrimination (N = ", n_distinct(vps_summary_grouped$participant), ")", sep=""), x = "Category", y = "% correct") +
  # scale_y_continuous("% correct", expand = c(0,0), limits = c(0,1.05)) +
  scale_x_discrete(NULL, labels = c("non-social","social")) +
  scale_fill_viridis_d("Condition", end = 0.25, begin = 0.25, guide = "none" ) + 
  theme_minimal() + 
  scale_fill_viridis_d() + 
  scale_color_viridis_d() +
  theme(legend.position = "none")
ggsave(file.path(path, "Plots", "Discrimination", "ga_correct.png"), width=1600, height=1600, units="px")

# vps_summary %>%
#   ez::ezANOVA(dv=.(correct_mean),
#               wid=.(participant),
#               within=.(social, condition),
#               # between=.(SPAI),
#               detailed=T, type=3) %>%
#   apa::anova_apa() %>%
#   ungroup()


# # Response Times
# vps_summary_grouped %>%
#   summarise(correct_rt = sd(rt_mean) * 100, correct_rt = mean(rt_mean) * 100, .by=c(social))
# 
# vps_summary_grouped %>%
#   t_test(data=., rt_mean ~ social, paired = T, alternative = "two.sided") %>% 
#   apa::t_apa()
# 
# vps_summary_grouped %>%
#   group_by(social) %>%
#   summarise(correct = mean(rt_mean, na.rm = T), se = sd(rt_mean)/sqrt(n())) %>%
#   ggplot(aes(x = social, y = correct, fill = social)) +
#   geom_col(position = position_dodge()) +
#   geom_errorbar(aes(ymax = correct + se, ymin = correct - se), width = 0.4) +
#   geom_point(data = vps_summary_grouped, aes(x = social, y = rt_mean), size = 2, shape = 21, color = "black", alpha=0.3, position = position_jitter(width=0.2, height=0.005)) +
#   labs(title = paste("Reaction Times (N = ", n_distinct(vps_summary_grouped$participant), ")", sep=""), x = "Category", y = "RT [ms]") +
#   scale_x_discrete("Category", labels = c("non-social","social")) +
#   scale_fill_viridis_d("Condition",end = 0.25,begin = 0.25, guide = "none" ) +
#   theme_minimal() +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   theme(legend.position = "none")
# ggsave(file.path(path, "Plots", "Discrimination", "ga_rt.png"), width=1800, height=2000, units="px")
# 
# vps_summary %>%
#   ez::ezANOVA(dv=.(rt_mean),
#               wid=.(participant),
#               within=.(social, condition),
#               # between=.(SPAI),
#               detailed=T, type=3) %>%
#   apa::anova_apa() %>%
#   ungroup()


vps.discrimination.score <- vps_summary %>%
  summarise(discrimination = mean(correct_mean), .by=c("participant", "social"))
names(vps.discrimination.score) = c("subject","condition_social","discrimination")
write.csv2(vps.discrimination.score, file.path(path, "Behavior", "discrimination.csv"), row.names=FALSE, quote=FALSE)
