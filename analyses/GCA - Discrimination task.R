
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
path.behavior = file.path(path, "gaze_discrimination_task", "data")

filemat = list.files(path.behavior, pattern = "*.csv") # read all files from data-folder into a single list


#cycle through all subject csv files

for (subject in filemat){
  #subject <- filemat[2]  #use  file for practice and coding changes
   
  start.time <- Sys.time()    
  
  # 1. General Data  
  
  vp_data <- read.csv(file.path(path.behavior, subject)) %>% #read data from files per subject
    filter(str_detect(trialtype,"cs_")) %>% 
    mutate(condition = ifelse(trialtype == comptype, "same", "different"),
           social = ifelse(trialtype %in% c("cs_plus_s","cs_minus_s"), "social","nonsocial")) %>%
    select(participant, condition, social, correctness, key.rt) %>%
    group_by(participant,social,condition) %>%
    summarise(correct_sum = sum(correctness == "correct"),
              correct_mean = sum(correctness == "correct")/n(),
              rt_mean = mean(key.rt,na.rm = T))
  
  
  participant <- vp_data$participant %>% unique()
  
  vp_summary <- vp_data
  
  # Create summary data
  
  vps_summary <- rbind(vps_summary,vp_summary)
  
  ggplot(vp_data, aes(x = social, fill = condition, y = correct_mean)) +
    geom_hline(yintercept = 0.5, linetype = 2, color = "red") + 
    geom_bar(stat = "identity", position = position_dodge(), color= "black") +
    scale_y_continuous("% correct", expand = c(0,0), limits = c(0,1)) +
    scale_x_discrete("Category") +
    scale_fill_viridis_d("Condition",end = 0.75) + 
    theme_classic()
  ggsave(file.path(path, "plots", "discrimination-task", "individuals", paste0(participant, ".png")))
    
  
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


vps_summary_long <- vps_summary

vps_summary_long %>% 
  group_by(social) %>%
  summarise(correct = mean(correct_mean, na.rm = T), se = sd(correct_mean)/sqrt(n())) %>% 
  ggplot(aes(x = social, y = correct, fill = social)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "red") + 
  geom_bar(stat = "identity", position = position_dodge(), color= "black") +
  geom_errorbar(aes(ymax = correct + se, ymin = correct - se), width = 0.4) +
  geom_point(data = vps_summary_long, aes(x = social, y = correct_mean), alpha = 0.5, shape = 21, position = position_jitter(0.2)) +
  scale_y_continuous("% correct", expand = c(0,0), limits = c(0,1)) +
  scale_x_discrete("Category", labels = c("non-social","social")) +
  scale_fill_viridis_d("Condition",end = 0.25,begin = 0.25, guide = "none" ) + 
  theme_classic()
ggsave(file.path(path, "plots", "discrimination-task", "ga_correct.png"))


# vps_summary_long %>% 
#   group_by(participant,social) %>%
#   summarise(correct_mean = mean(correct_mean, na.rm = T)) %>% 
#   ggplot(aes(x = social, y = correct_mean, fill = social)) + facet_wrap(~participant) + 
#   geom_hline(yintercept = 0.5, linetype = 2, color = "red") + 
#   geom_bar(stat = "identity", position = position_dodge(), color= "black" ) +
#   scale_y_continuous("% correct", expand = c(0,0), limits = c(0,1)) +
#   scale_x_discrete("Category", labels = c("non-social","social")) +
#   scale_fill_viridis_d("Condition",end = 0.25,begin = 0.25, guide = "none" ) + 
#   theme_classic()
# ggsave(file.path(path, "plots", "discrimination-task", "ga_correct_ind.png"))


# Response Times

vps_summary_long %>% 
  group_by(social) %>%
  summarise(correct = mean(rt_mean, na.rm = T), se = sd(rt_mean)/sqrt(n())) %>% 
  ggplot(aes(x = social, y = correct, fill = social)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "red") + 
  geom_bar(stat = "identity", position = position_dodge(), color= "black") +
  geom_errorbar(aes(ymax = correct + se, ymin = correct - se), width = 0.4) +
  geom_point(data = vps_summary_long, aes(x = social, y = rt_mean), alpha = 0.5, shape = 21, position = position_jitter(0.2)) +
  scale_y_continuous("RT [s]", expand = c(0,0), limits = c(0,1)) +
  scale_x_discrete("Category", labels = c("non-social","social")) +
  scale_fill_viridis_d("Condition",end = 0.25,begin = 0.25, guide = "none" ) + 
  theme_classic()
ggsave(file.path(path, "plots", "discrimination-task", "ga_rt.png"))


# vps_summary_long %>% 
#   group_by(participant,social) %>%
#   summarise(rt_mean = mean(rt_mean, na.rm = T)) %>% 
#   ggplot(aes(x = social, y = rt_mean, fill = social)) + facet_wrap(~participant) + 
#   geom_hline(yintercept = 0.5, linetype = 2, color = "red") + 
#   geom_bar(stat = "identity", position = position_dodge(), color= "black" ) +
#   scale_y_continuous("RT [s]", expand = c(0,0), limits = c(0,1)) +
#   scale_x_discrete("Category", labels = c("non-social","social")) +
#   scale_fill_viridis_d("Condition",end = 0.25,begin = 0.25, guide = "none" ) + 
#   theme_classic()
# ggsave(file.path(path, "plots", "discrimination-task", "ga_rt_ind.png"))


# vps_summary_long %>% 
#   group_by(participant,social) %>%
#   summarise(rt_mean = mean(rt_mean, na.rm = T)) %>% 
#   ggplot(aes(x = social, y = rt_mean, fill = social)) + facet_wrap(~participant) + 
#   geom_hline(yintercept = 0.5, linetype = 2, color = "red") + 
#   geom_violin(color= "black" ) +
#   scale_y_continuous("RT [s]", expand = c(0,0), limits = c(0,1)) +
#   scale_x_discrete("Category", labels = c("non-social","social")) +
#   scale_fill_viridis_d("Condition",end = 0.25,begin = 0.25, guide = "none" ) + 
#   theme_classic()
# ggsave(file.path(path, "plots", "discrimination-task", "ga_rt_ind.png"))

