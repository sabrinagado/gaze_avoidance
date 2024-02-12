###############################################################################
# Gaze Contingent Avoidance Project
# Sabrina Gado & Yannik Stegmann
# Code adapted from Daniel Gromer
###############################################################################

library(tidyverse)
library(mc2d) # for pempiricalD()
library(afex)

# For parallelization
library(parallel)
library(doSNOW)

sample_rate = 50

# Functions --------------------------------------------------------------------

## General functions -----------------------------------------------------------

# Find clusters in the vector of F-statistics. `x` must be a call to
# pt_Ftest_statistic() and critical_value` must be a call to pt_critical_F() 
# respectively.
# 
# Returns a tibble with the following columns:
# `start`  index of the start of the cluster(s)
# `end`    index of the end of the cluster(s)
# `length` length of the cluster(s)
# (all indices with respect to the series of time points)
clusters_over_time <- function(x, critical_value)
  
{
  # Find clusters in the vector that are above/below the threshold
  x <- rle(abs(x) > critical_value)
  
  # Extract only clusters that are above the threshold
  cluster_lengths <- x$lengths[x$values]
  
  # Find the start and end of clusters via the cumsum of cluster lengths. For
  # cluster start, we'll put a 0/FALSE in front, so that we don't run into
  # indexing vec[0]
  cluster_start <- c(0L, cumsum(x$lengths))[which(c(FALSE, x$values)) - 1] + 1L
  cluster_end <- cumsum(x$lengths)[x$values]
  
  tibble(cluster = seq_along(cluster_lengths), start = cluster_start,
         end = cluster_end, length = cluster_lengths)
}

# clusters <- clusters_over_time(t, pt_critical_t_between(tmp, id = "id"))

## Functions for within-subject ANOVA (two factors) -------------------

# This is an implementation of a cluster-based permutation test for within-
# subject comparisons of 2d data (signal x time)
# 
# The data needs to be in a tibble in long format with the following columns and
# must not have any missing data (NAs or missing rows, i.e., each subject needs
# to have 2 * "number of sample" rows):
# 
# (1) A column containing the dependent variable (dv) (e.g., heart rate, skin
#     conductance)
# (2) A column indicating the time points / samples (time)
# (3) Two columns identifying the ANOVA factors (factor1, factor2)
# (4) A column containing the participant identifier (id)
#
# All functions expect the column names as character.

# Calculate the anova for each time point
perform_anova <- function(df, dv, id, factor1, factor2) {
  anova <- aov_ez(df, dv=dv, id=id, within = c(factor1, factor2))
  anova <- anova$anova_table
  F_values <- anova %>% select("F")
  return(t(F_values))
}

pt_Ftest_statistic <- function(data, dv, id, factor1, factor2, time)
{
  data_stat <- data %>%
    group_by(!!sym(time)) %>% 
    group_modify(~ as_tibble(perform_anova(.x, dv=dv, id=id, factor1=factor1, factor2=factor2)))
  colnames(data_stat) <- c("samplepoint", "F_main1", "F_main2", "F_int")
  return(data_stat)
}

# Calculate the critical F-value for a given sample size and factor levels
pt_critical_F <- function(data, id, factor1, factor2, alpha = .05)
{
  p = length(unique(data[[factor1]]))
  q = length(unique(data[[factor2]]))
  main1 = qf(1 - alpha / 2, p - 1, length(unique(data[[id]])) - 1)
  main2 = qf(1 - alpha / 2, q - 1, length(unique(data[[id]])) - 1)
  int = qf(1 - alpha / 2, (p-1) * (q-1), length(unique(data[[id]])) - 1)
  return(c(main1, main2, int))
}

# Calculation of the null distribution (H0) of cluster lengths
pt_null_distribution <- function(data, dv, within, factor1, factor2, trial, time, id, nperm = 100)
{
  # Create a vector to store the null distribution of cluster lengths
  null_distribution_F_main1 <- vector("integer", nperm)
  null_distribution_F_main2 <- vector("integer", nperm)
  null_distribution_F_int <- vector("integer", nperm)
  
  critical_F_values <- pt_critical_F(data, id = id, factor1 = factor1, factor2 = factor2)
  critical_F_main1 <- critical_F_values[[1]] 
  critical_F_main2 <- critical_F_values[[2]]
  critical_F_int <- critical_F_values[[3]]
  
  # Make sure that `data` is arranged correctly
  data <- data |> arrange(!!!syms(c(id, within, time)))
  
  for (i in seq_len(nperm))
  {
    cat(i, " - ")
    # Create a random permutation of the conditions (= permute the `within` column)
    shuffled_conds <- data %>%
      # Ensure there's a unique row for each combination of id and trial for shuffling
      distinct(!!sym(id), !!sym(trial), .keep_all = TRUE) %>%
      # Shuffle the condition labels within each subject
      group_by(!!sym(id)) %>%
      mutate(shuffled_condition = sample(condition, size = n(), replace = FALSE)) %>%
      ungroup() %>% 
      select(c(!!sym(id), !!sym(trial), shuffled_condition))
    
    data_shuffled <- data %>% 
      left_join(shuffled_conds, by = c(id, trial))
    
    data_shuffled <- data_shuffled %>% 
      mutate(shuffled_condition_social = if_else(str_detect(shuffled_condition, "non-social"), "non-social", "social")) %>% 
      mutate(shuffled_condition_threat = if_else(str_detect(shuffled_condition, "pos"), "pos", "neg"))
    
    # Find clusters in the random permutation
    F_dist <- pt_Ftest_statistic(data_shuffled, dv=dv, id=id, factor1=paste0("shuffled_", factor1), factor2=paste0("shuffled_", factor2), time=time)
    clusters_F_main1 <- clusters_over_time(F_dist$F_main1, critical_F_main1)
    clusters_F_main2 <- clusters_over_time(F_dist$F_main2, critical_F_main2)
    clusters_F_int <- clusters_over_time(F_dist$F_int, critical_F_int)
    
    # If clusters were found, take the one with the maximum length and save it
    null_distribution_F_main1[i] <- if (nrow(clusters_F_main1) > 0) max(clusters_F_main1$length) else 0
    null_distribution_F_main2[i] <- if (nrow(clusters_F_main2) > 0) max(clusters_F_main2$length) else 0
    null_distribution_F_int[i] <- if (nrow(clusters_F_int) > 0) max(clusters_F_int$length) else 0
  }
  
  null_distributions <- data.frame(null_distribution_F_main1, null_distribution_F_main2, null_distribution_F_int)
  colnames(null_distributions) <- c("F_main1", "F_main2", "F_int")
  return(null_distributions)
}

# Prepare data -----------------------------------------------------------------

ga_unified <- readRDS("ET_pupil_ga_unified.RData")
pupil_df <- readRDS("ET_pupil_df.Rdata")

responder <- pupil_df %>% filter(condition %in% c(1,2,3)) %>% group_by(ID) %>% 
  summarise(dilators = mean(as.numeric(interpolated))) %>% filter(dilators < 0.5) %>% .$ID

rep_str = c('10' = "Shock",
            '11' = "Reward", 
            '12' = "No Feedback",
            '2' = 'CSneg, social, Acq',
            '3' = 'CSpos, social, Acq',
            '4' = 'CSneg, non-social, Acq',
            '5' = 'CSpos, non-social, Acq',
            '6' = 'CSneg, social, Test',
            '7' = 'CSpos, social, Test',
            '8' = 'CSneg, non-social, Test',
            '9' = 'CSpos, non-social, Test')

pupil_long <- ga_unified %>% 
  filter(valid == TRUE) %>%
  filter(ID %in% responder) %>%
  .$diameter %>% bind_rows() %>%
  mutate(condition = as.factor(condition)) %>%
  filter(condition %in% c(6,7,8,9)) %>% 
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

# Calculate null distributions -------------------------------------------------
iterations <- 1000

# with parallelization:
# Use parallelization to run the permutations. As this was written for Windows, the doSNOW backend is used in the current implementation.
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, list("pt_null_distribution", "pt_critical_F", "pt_Ftest_statistic", "perform_anova"))
registerDoSNOW(cl)

pupil_null_dist <- data.frame()

tmp <- foreach(i = 1:length(cl), .combine = "rbind", .packages = c("tidyverse", "afex")) %dopar% {
  y <- pt_null_distribution(pupil_long, dv = "diameter", within = "condition", factor1 = "condition_threat", factor2 = "condition_social",
                                   time = "samplepoint", trial = "trial", id = "subject", nperm = as.integer(ceiling(iterations/length(cl))))
  return(y)
}
stopCluster(cl)

pupil_null_dist <- tmp[1:iterations,]

# # without parallelization:
# pupil_null_dist <- pt_null_distribution(pupil_long, dv = "diameter", within = "condition", factor1 = "condition_threat", factor2 = "condition_social",
#                                                time = "time", trial = "trial", id = "subject", nperm = iterations)

# Create a tibble with the 95% quantile of cluster lengths under H0
pupil_critical_length_null_F_main1 <- quantile(pupil_null_dist$F_main1, .95)
pupil_critical_length_null_F_main2 <- quantile(pupil_null_dist$F_main2, .95)
pupil_critical_length_null_F_int <- quantile(pupil_null_dist$F_int, .95)

pupil_critical_lengths_null <- data.frame(pupil_critical_length_null_F_main1, pupil_critical_length_null_F_main2, pupil_critical_length_null_F_int)
colnames(pupil_critical_lengths_null) <- c("F_main1", "F_main2", "F_int")

# Save distributions and stop cluster ----------------------------------------
save(pupil_null_dist, pupil_critical_lengths_null,
     file = "pupil_null_distributions.RData")

# Cluster-based permutation tests ---------------------------------------------
load("pupil_null_distributions.RData")

# Find clusters
F_tests <- pt_Ftest_statistic(pupil_long, dv = "diameter", id = "subject", factor1 = "condition_threat", factor2 = "condition_social", time = "samplepoint")
criticalFs <- pt_critical_F(pupil_long, id = "subject", factor1 = "condition_threat", factor2 = "condition_social")
critical_F_main1 <- criticalFs[[1]] 
critical_F_main2 <- criticalFs[[2]]
critical_F_int <- criticalFs[[3]]

clusters_F_main1 <- clusters_over_time(F_tests$F_main1, critical_F_main1)
clusters_F_main1 <- clusters_F_main1 %>% 
  mutate(critical_length = pupil_critical_lengths_null$F_main1) %>% 
  mutate(p = 1 - pempiricalD(length, critical_length)) %>% 
  filter(p < .05) %>% 
  mutate(times_start = start / sample_rate - 0.5,
         times_end = end / sample_rate - 0.5)

clusters_F_main2 <- clusters_over_time(F_tests$F_main2, critical_F_main2)
clusters_F_main2 <- clusters_F_main2 %>% 
  mutate(critical_length = pupil_critical_lengths_null$F_main2) %>% 
  mutate(p = 1 - pempiricalD(length, critical_length)) %>% 
  filter(p < .05) %>% 
  mutate(times_start = start / sample_rate - 0.5,
         times_end = end / sample_rate - 0.5)

clusters_F_int <- clusters_over_time(F_tests$F_int, critical_F_int)
clusters_F_int <- clusters_F_int %>% 
  mutate(critical_length = pupil_critical_lengths_null$F_int) %>% 
  mutate(p = 1 - pempiricalD(length, critical_length)) %>% 
  filter(p < .05) %>% 
  mutate(times_start = start / sample_rate - 0.5,
         times_end = end / sample_rate - 0.5)

save(clusters_F_main1, clusters_F_main2, clusters_F_int,
     file = "pupil_cluster.RData")

load("pupil_cluster.RData")

ga_unified %>% filter(valid == TRUE) %>%
  filter(ID %in% responder) %>%
  .$diameter %>% bind_rows() %>%
  mutate(condition = as.factor(condition)) %>%
  filter(condition %in% c(6,7,8,9)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  group_by(condition,samplepoint) %>% 
  summarise(diameter.mean = mean(diameter), diameter.se = sd(diameter)/sqrt(length(diameter)), time = mean(time)) %>%
  {ggplot(., aes(x=time, y=diameter.mean)) +
      geom_vline(xintercept=0, color="black",linetype="solid") + #zero = picture onset
      geom_vline(xintercept=10, color="black",linetype="solid") + #picture offset
      geom_line(aes(colour=condition)) +
      geom_ribbon(aes(ymin=diameter.mean-diameter.se, ymax=diameter.mean+diameter.se, colour=condition, fill=condition), color = NA, alpha=.2) +
      geom_segment(data = clusters_F_main2, aes(x=times_start, xend = times_end, y=-90, yend=-90, size="Main Effect Social"), colour = "#ff8383", linewidth = 1, inherit.aes=FALSE) +
      # geom_segment(data = clusters_F_main1, aes(x=times_start, xend = times_end, y=-93, yend=-93, size="Main Effect Threat"), colour = "#e874ff", linewidth = 1, inherit.aes=FALSE) +
      # geom_segment(data = clusters_F_int, aes(x=times_start, xend = times_end, y=-96, yend=-96, size="Social x Threat Interaction "), colour = "#ffdd74", linewidth = 1, inherit.aes=FALSE) +
      scale_x_continuous("Time [s]",limits=c(-0.5, 11), minor_breaks=c(0,1,2,3,4,5,6,7,8,9,10), breaks=c(0, 2, 4, 6, 8, 10)) +
      scale_y_continuous("Pupil Diameter", breaks=c(-80,-40, 0, 40), minor_breaks=c(-80, -60, -40, -20, 0, 20, 40, 60)) +
      scale_color_viridis_d(aesthetics = c("colour", "fill")) +
      theme_bw() +
      scale_size_manual("effects", values=rep(1,4), guide=guide_legend(override.aes = list(colour=c("#ff8383")))) # , "#e874ff", "#ffdd74"
  }
ggsave("../plots/Pupil/cs_test_cluster.png",type="cairo-png", width=2500/400, height=1080/300, dpi=300)


# Bins
pupil_df_long_test <- pupil_df %>%
  filter(ID %in% responder) %>%
  pivot_longer(dilation_0:dilation_19, names_to = "timebin", values_to ="diameter") %>%
  separate(timebin,c("quark","timebin")) %>% mutate(timebin = as.numeric(timebin)) %>%
  # filter(timebin > 4 & timebin < 21) %>% #filter(valid == T) %>%
  filter(condition %in% c(6,7,8,9)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  select(ID, trial, condition, timebin, diameter) %>% 
  mutate(diameter = as.numeric(diameter[,1])) %>% 
  mutate(time = timebin * 0.5) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

main_effect_threat = list()
main_effect_social = list()
interaction_effect = list()
alpha = .05 / length(unique(pupil_df_long_test$time))
for (timepoint in unique(pupil_df_long_test$time)) {
  # timepoint = 0
  data = pupil_df_long_test %>% 
    filter(time == timepoint) %>% 
    mutate(ID = as.factor(ID), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
    group_by(ID, condition_social, condition_threat)
  
  # data %>% aov_ez(dv="diameter", id="ID", within = c("condition_social","condition_threat"))
  
  model = lmer(diameter ~ condition_social + condition_threat + condition_social:condition_threat + (1|ID), data)
  anova = anova(model, type=2)
  
  p_threat = anova["condition_threat", "Pr(>F)"]
  p_social = anova["condition_social", "Pr(>F)"]
  p_interaction = anova["condition_social:condition_threat", "Pr(>F)"]
  
  if (p_threat < alpha) {
    main_effect_threat = c(main_effect_threat, timepoint)
  }
  if (p_social < alpha) {
    main_effect_social = c(main_effect_social, timepoint)
  }
  if (p_interaction < alpha) {
    interaction_effect = c(interaction_effect, timepoint)
  }
}
main_effect_threat = as.numeric(main_effect_threat)
main_effect_social = as.numeric(main_effect_social)
interaction_effect = as.numeric(interaction_effect)

main_effect_threat = tibble(times_start = as.numeric(main_effect_threat), times_end = as.numeric(main_effect_threat) + 0.5)
main_effect_social = tibble(times_start = as.numeric(main_effect_social), times_end = as.numeric(main_effect_social) + 0.5)
interaction_effect = tibble(times_start = as.numeric(interaction_effect), times_end = as.numeric(interaction_effect) + 0.5)

ga_unified %>% filter(valid == TRUE) %>%
  filter(ID %in% responder) %>%
  .$diameter %>% bind_rows() %>%
  mutate(condition = as.factor(condition)) %>% 
  # filter(trial < 60) %>%
  filter(condition %in% c(6,7,8,9)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  group_by(condition,samplepoint) %>% 
  summarise(diameter.mean = mean(diameter), diameter.se = sd(diameter)/sqrt(length(diameter)), time = mean(time)) %>%
  {ggplot(., aes(x=time, y=diameter.mean)) +
      geom_vline(xintercept=0, color="black",linetype="solid") + #zero = picture onset
      geom_vline(xintercept=10, color="black",linetype="solid") + #picture offset
      geom_line(aes(colour=condition)) +
      geom_ribbon(aes(ymin=diameter.mean-diameter.se, ymax=diameter.mean+diameter.se, colour=condition, fill=condition), color = NA, alpha=.2) +
      geom_segment(data = main_effect_social, aes(x=times_start, xend = times_end, y=-90, yend=-90, size="Main Effect Social"), colour = "#ff8383", linewidth = 1, inherit.aes=FALSE) +
      geom_segment(data = main_effect_threat, aes(x=times_start, xend = times_end, y=-93, yend=-93, size="Main Effect Threat"), colour = "#e874ff", linewidth = 1, inherit.aes=FALSE) +
      # geom_segment(data = interaction_effect, aes(x=times_start, xend = times_end, y=-96, yend=-96, size="Social x Threat Interaction "), colour = "#ffdd74", linewidth = 1, inherit.aes=FALSE) +
      scale_x_continuous("Time [s]",limits=c(-0.5, 11), minor_breaks=c(0,1,2,3,4,5,6,7,8,9,10), breaks=c(0, 2, 4, 6, 8, 10)) +
      scale_y_continuous("Pupil Diameter", breaks=c(-80,-40, 0, 40), minor_breaks=c(-80, -60, -40, -20, 0, 20, 40, 60)) +
      scale_color_viridis_d(aesthetics = c("colour", "fill")) +
      theme_bw() +
      scale_size_manual("effects", values=rep(1,4), guide=guide_legend(override.aes = list(colour=c("#ff8383", "#e874ff")))) # , "#ffdd74"
  }
ggsave("../plots/Pupil/cs_test_bins.png",type="cairo-png", width=2500/400, height=1080/300, dpi=300)

# # Acquisition
# # Long format for statistical testing
# pupil_df_long_acq <- pupil_df %>%
#   # filter(outcome == "no outcome") %>%
#   #filter(ID %in% responder) %>%
#   pivot_longer(dilation_0:dilation_15, names_to = "timebin", values_to ="diameter") %>%
#   separate(timebin,c("quark","timebin")) %>% mutate(timebin = as.numeric(timebin)) %>%
#   # filter(timebin > 4 & timebin < 21) %>% #filter(valid == T) %>%
#   filter(condition %in% c(2,3,4,5)) %>%
#   mutate(across('condition', str_replace_all, rep_str)) %>% 
#   select(ID, trial, condition, outcome, timebin, diameter) %>% 
#   mutate(diameter = as.numeric(diameter[,1])) %>% 
#   mutate(time = timebin * 0.5) %>% 
#   mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
#   mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))
# 
# main_effect_threat = list()
# main_effect_social = list()
# interaction_effect = list()
# alpha = .05 / length(unique(pupil_df_long_acq$time))
# for (timepoint in unique(pupil_df_long_acq$time)) {
#   # timepoint = 3
#   data = pupil_df_long_acq %>% 
#     filter(time == timepoint) %>% 
#     mutate(ID = as.factor(ID), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
#     group_by(ID, condition_social, condition_threat)
#   model = lmer(diameter ~ condition_social + condition_threat + condition_social:condition_threat + (1|ID), data)
#   anova = anova(model, type=2)
#   
#   p_threat = anova["condition_threat", "Pr(>F)"]
#   p_social = anova["condition_social", "Pr(>F)"]
#   p_interaction = anova["condition_social:condition_threat", "Pr(>F)"]
#   if (p_threat < alpha) {
#     main_effect_threat = c(main_effect_threat, timepoint)
#   }
#   if (p_social < alpha) {
#     main_effect_social = c(main_effect_social, timepoint)
#   }
#   if (p_interaction < alpha) {
#     interaction_effect = c(interaction_effect, timepoint)
#   }
# }
# main_effect_threat = as.numeric(main_effect_threat)
# main_effect_social = as.numeric(main_effect_social)
# interaction_effect = as.numeric(interaction_effect)
# 
# main_effect_threat = tibble(times_start = as.numeric(main_effect_threat), times_end = as.numeric(main_effect_threat) + 0.5)
# main_effect_social = tibble(times_start = as.numeric(main_effect_social), times_end = as.numeric(main_effect_social) + 0.5)
# interaction_effect = tibble(times_start = as.numeric(interaction_effect), times_end = as.numeric(interaction_effect) + 0.5)
# 
# ga_unified %>% filter(valid == TRUE) %>%
#   # filter(outcome == "no outcome") %>%
#   #filter(ID %in% responder) %>%
#   .$diameter %>% bind_rows() %>%
#   mutate(condition = as.factor(condition)) %>%
#   filter(condition %in% c(2,3,4,5)) %>%
#   mutate(across('condition', str_replace_all, rep_str)) %>% 
#   group_by(condition,samplepoint) %>% 
#   summarise(diameter.mean = mean(diameter), diameter.se = sd(diameter)/sqrt(length(diameter)), time = mean(time)) %>%
#   {ggplot(., aes(x=time, y=diameter.mean)) +
#       geom_vline(xintercept=0, colour="black",linetype="solid") + #zero
#       geom_line(aes(colour=condition)) +
#       geom_ribbon(aes(ymin=diameter.mean-diameter.se, ymax=diameter.mean+diameter.se, colour=condition, fill=condition), color = NA, alpha=.2) +
#       geom_segment(data = main_effect_social, aes(x=times_start, xend = times_end, y=-70, yend=-70, size="Main Effect Social"), colour = "#ff8383", linewidth = 1, inherit.aes=FALSE) +
#       geom_segment(data = main_effect_threat, aes(x=times_start, xend = times_end, y=-75, yend=-75, size="Main Effect Threat"), colour = "#e874ff", linewidth = 1, inherit.aes=FALSE) +
#       geom_segment(data = interaction_effect, aes(x=times_start, xend = times_end, y=-80, yend=-80, size="Social x Threat Interaction "), colour = "#ffdd74", linewidth = 1, inherit.aes=FALSE) +
#       scale_x_continuous("Time [s]",limits=c(-0.5, 8), minor_breaks=c(0,1,2,3,4,5,6,7,8), breaks=c(0, 2, 4, 6, 8)) +
#       scale_y_continuous("Pupil Diameter",limits=c(-80, 150)) +
#       scale_color_viridis_d(aesthetics = c("colour", "fill")) +
#       theme_bw() + 
#       # labs(title = paste("Pupil (N = ", n_distinct(pupil_df_long_acq$ID), ")", sep="")) +
#       scale_size_manual("effects", values=rep(1,3), guide=guide_legend(override.aes = list(colour=c("#ff8383", "#e874ff", "#ffdd74")))) # 
#   }
# ggsave("../plots/Pupil/cs_acq_bins.png",type="cairo-png", width=2500/400, height=1080/300, dpi=300)
