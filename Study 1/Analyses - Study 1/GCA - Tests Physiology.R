###############################################################################
# Gaze Contingent Avoidance Project
# Sabrina Gado & Yannik Stegmann
# Code adapted from Daniel Gromer
###############################################################################

library(tidyverse)
library(mc2d) # for pempiricalD()
library(afex)
library(lme4)
library(lmerTest)
library(apa)
library(cowplot)

# For parallelization
library(parallel)
library(doSNOW)

theme_set(theme_bw(base_size = 16))

iterations <- 1000

rep_str = c('10' = "Shock",
            '11' = "Reward", 
            '12' = "No Feedback",
            '2' = 'CSneg, social, Acq',
            '3' = 'CSpos, social, Acq',
            '4' = 'CSneg, non-social, Acq',
            '5' = 'CSpos, non-social, Acq',
            '6' = 'CSneg, social',
            '7' = 'CSpos, social',
            '8' = 'CSneg, non-social',
            '9' = 'CSpos, non-social')

{ # Functions ---
  
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
    anova <- aov_ez(df, dv=dv, id=id, within = c(factor1, factor2), na.rm = TRUE)
    anova <- anova$anova_table
    F_values <- anova %>% select("F")
    return(t(F_values))
  }
  
  perform_lmm <- function(df, dv, id, factor1, factor2) {
    formula_lmm <- as.formula(paste0(dv, " ~ ", factor1, " + ", factor2," + ", factor1, ":", factor2, "+ (1|", id, ")"))
    model = lmer(formula_lmm, df)
    anova = anova(model, type=2)
    F_values <- anova %>% select("F value")
    return(t(F_values))
  }
  
  pt_Ftest_statistic <- function(data, dv, id, factor1, factor2, time)
  {
    data_stat <- data %>%
      group_by(!!sym(time)) %>% 
      group_modify(~ as_tibble(perform_lmm(.x, dv=dv, id=id, factor1=factor1, factor2=factor2)))
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
    
    pb <- txtProgressBar(min = 1, max = nperm, style = 3)
    
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
      
      setTxtProgressBar(pb, i)
    }
    
    close(pb)
    
    null_distributions <- data.frame(null_distribution_F_main1, null_distribution_F_main2, null_distribution_F_int)
    colnames(null_distributions) <- c("F_main1", "F_main2", "F_int")
    return(null_distributions)
  }
}

################################################################################
# PUPIL ------------------------------------------------------------------------
################################################################################
cat("\n", "Cluster-based permutation tests of pupil data", "\n")
sample_rate = 50

# Prepare data -----------------------------------------------------------------
ga_unified <- readRDS(file.path("Study 1", "Physio", "ET_pupil_ga_unified.RData"))
pupil_df <- readRDS(file.path("Study 1", "Physio", "ET_pupil_df.Rdata"))

responder <- pupil_df %>% filter(condition %in% c(1,2,3)) %>% group_by(ID) %>% 
  summarise(dilators = mean(as.numeric(interpolated))) %>% filter(dilators < 0.5) %>% .$ID

# CLUSTER
pupil_long <- ga_unified %>% 
  filter(valid == TRUE) %>%
  filter(ID %in% responder) %>%
  .$diameter %>% bind_rows() %>%
  mutate(condition = as.factor(condition)) %>%
  filter(condition %in% c(6,7,8,9)) %>% 
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>% 
  filter(time < 10)

# Calculate null distributions -------------------------------------------------
# with parallelization:
# Use parallelization to run the permutations. As this was written for Windows, the doSNOW backend is used in the current implementation.
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, list("pt_null_distribution", "pt_critical_F", "pt_Ftest_statistic", "perform_anova", "perform_lmm"))
registerDoSNOW(cl)

pupil_null_dist <- data.frame()

tmp <- foreach(i = 1:length(cl), .combine = "rbind", .packages = c("tidyverse", "afex", "lme4")) %dopar% {
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
save(pupil_null_dist, pupil_critical_lengths_null, file = file.path("Study 1", "Physio", "pupil_null_distributions.RData"))

# Cluster-based permutation tests ---------------------------------------------
load(file.path("Study 1", "Physio", "pupil_null_distributions.RData"))

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

save(clusters_F_main1, clusters_F_main2, clusters_F_int, file = file.path("Study 1", "Physio", "pupil_cluster.RData"))

load(file.path("Study 1", "Physio", "pupil_cluster.RData"))

plot_pupil_cluster <- ga_unified %>% filter(valid == TRUE) %>%
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
      geom_segment(data = clusters_F_main2, aes(x=times_start, xend = times_end, y=-.11, yend=-.11, size="Main Effect Stimulus Type"), colour = "#ff8383", linewidth = 1, inherit.aes=FALSE) +
      # geom_segment(data = clusters_F_main1, aes(x=times_start, xend = times_end, y=-93, yend=-93, size="Main Effect Conditioning"), colour = "#e874ff", linewidth = 1, inherit.aes=FALSE) +
      # geom_segment(data = clusters_F_int, aes(x=times_start, xend = times_end, y=-96, yend=-96, size="Stimulus Type x Conditioning "), colour = "#ffdd74", linewidth = 1, inherit.aes=FALSE) +
      scale_x_continuous("Time [s]",limits=c(-0.5, 11), minor_breaks=c(0,1,2,3,4,5,6,7,8,9,10), breaks=c(0, 2, 4, 6, 8, 10)) +
      scale_y_continuous("∆ Pupil Diameter [mm]", limits=c(-0.115, 0.1), breaks=c(-0.1, -0.05, 0, 0.05, 0.1), minor_breaks=c(-0.1, -0.075, -0.05, -0.025, 0, 0.025, 0.05, 0.075, 0.1)) +
      scale_color_viridis_d(NULL, aesthetics = c("colour", "fill")) +
      # theme_bw() +
      scale_size_manual(NULL, values=rep(1,4), guide=guide_legend(override.aes = list(colour=c("#ff8383")))) # , "#e874ff", "#ffdd74"
  }

# ggsave(file.path("Study 1", "Plots", "Pupil", "cs_test_cluster.png") ,type="cairo-png", width=2500/400, height=1080/300, dpi=300)

################################################################################
# EDA --------------------------------------------------------------------------
################################################################################
cat("\n", "Cluster-based permutation tests of EDA data", "\n")
sample_rate = 50

# Prepare data -----------------------------------------------------------------

eda_unified <- readRDS(file.path("Study 1", "Physio", "EDA_unified.RData"))
eda_df <- readRDS(file.path("Study 1", "Physio", "EDA_df.Rdata"))

responder <- eda_df %>%
  mutate(scl_bl = abs(SCL-Baseline)) %>%
  group_by(ID) %>%
  summarise(scl_avg = mean(scl_bl)) %>% 
  filter(scl_avg >= .02) %>% .$ID

# CLUSTER
eda_long <- eda_unified %>%
  filter(ID %in% responder) %>%
  .$EDA %>% bind_rows() %>%
  mutate(condition = as.factor(condition)) %>%
  filter(condition %in% c(6,7,8,9)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>%
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

# Calculate null distributions -------------------------------------------------
# with parallelization:
# Use parallelization to run the permutations. As this was written for Windows, the doSNOW backend is used in the current implementation.
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, list("pt_null_distribution", "pt_critical_F", "pt_Ftest_statistic", "perform_anova", "perform_lmm"))
registerDoSNOW(cl)

eda_null_dist <- data.frame()

tmp <- foreach(i = 1:length(cl), .combine = "rbind", .packages = c("tidyverse", "afex", "lme4")) %dopar% {
  y <- pt_null_distribution(eda_long, dv = "EDA", within = "condition", factor1 = "condition_threat", factor2 = "condition_social",
                            time = "samplepoint", trial = "trial", id = "subject", nperm = as.integer(ceiling(iterations/length(cl))))
  return(y)
}
stopCluster(cl)

eda_null_dist <- tmp[1:iterations,]

# # without parallelization:
# eda_null_dist <- pt_null_distribution(eda_long, dv = "diameter", within = "condition", factor1 = "condition_threat", factor2 = "condition_social",
#                                       time = "time", trial = "trial", id = "subject", nperm = iterations)

# Create a tibble with the 95% quantile of cluster lengths under H0
eda_critical_length_null_F_main1 <- quantile(eda_null_dist$F_main1, .95)
eda_critical_length_null_F_main2 <- quantile(eda_null_dist$F_main2, .95)
eda_critical_length_null_F_int <- quantile(eda_null_dist$F_int, .95)

eda_critical_lengths_null <- data.frame(eda_critical_length_null_F_main1, eda_critical_length_null_F_main2, eda_critical_length_null_F_int)
colnames(eda_critical_lengths_null) <- c("F_main1", "F_main2", "F_int")

# Save distributions and stop cluster ----------------------------------------
save(eda_null_dist, eda_critical_lengths_null, file = file.path("Study 1", "Physio", "eda_null_distributions.RData"))

# Cluster-based permutation tests ---------------------------------------------
load(file.path("Study 1", "Physio", "eda_null_distributions.RData"))

# Find clusters
F_tests <- pt_Ftest_statistic(eda_long, dv = "EDA", id = "subject", factor1 = "condition_threat", factor2 = "condition_social", time = "samplepoint")
criticalFs <- pt_critical_F(eda_long, id = "subject", factor1 = "condition_threat", factor2 = "condition_social")
critical_F_main1 <- criticalFs[[1]] 
critical_F_main2 <- criticalFs[[2]]
critical_F_int <- criticalFs[[3]]

clusters_F_main1 <- clusters_over_time(F_tests$F_main1, critical_F_main1)
clusters_F_main1 <- clusters_F_main1 %>% 
  mutate(critical_length = eda_critical_lengths_null$F_main1) %>% 
  mutate(p = 1 - pempiricalD(length, critical_length)) %>% 
  filter(p < .05) %>% 
  mutate(times_start = start / sample_rate - 0.5,
         times_end = end / sample_rate - 0.5)

clusters_F_main2 <- clusters_over_time(F_tests$F_main2, critical_F_main2)
clusters_F_main2 <- clusters_F_main2 %>% 
  mutate(critical_length = eda_critical_lengths_null$F_main2) %>% 
  mutate(p = 1 - pempiricalD(length, critical_length)) %>% 
  filter(p < .05) %>% 
  mutate(times_start = start / sample_rate - 0.5,
         times_end = end / sample_rate - 0.5)

clusters_F_int <- clusters_over_time(F_tests$F_int, critical_F_int)
clusters_F_int <- clusters_F_int %>% 
  mutate(critical_length = eda_critical_lengths_null$F_int) %>% 
  mutate(p = 1 - pempiricalD(length, critical_length)) %>% 
  filter(p < .05) %>% 
  mutate(times_start = start / sample_rate - 0.5,
         times_end = end / sample_rate - 0.5)

save(clusters_F_main1, clusters_F_main2, clusters_F_int, file = file.path("Study 1", "Physio", "eda_cluster.RData"))

load(file.path("Study 1", "Physio", "eda_cluster.RData"))

plot_eda_cluster <- eda_unified %>%
  filter(ID %in% responder) %>%
  .$EDA %>% bind_rows() %>%
  mutate(condition = as.factor(condition)) %>% 
  filter(condition %in% c(6,7,8,9)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  group_by(condition,sample) %>% 
  summarise(EDA.mean = mean(EDA), EDA.se = sd(EDA)/sqrt(length(EDA)), time = mean(time)) %>%
  {ggplot(., aes(x=time, y=EDA.mean)) +
      geom_vline(xintercept=0, color="black",linetype="solid") + #zero = picture onset
      geom_vline(xintercept=10, color="black",linetype="solid") + #picture offset
      geom_line(aes(colour=condition)) +
      geom_ribbon(aes(ymin=EDA.mean-EDA.se, ymax=EDA.mean+EDA.se, colour=condition, fill=condition), color = NA, alpha=.2) +
      # geom_segment(data = clusters_F_main2, aes(x=times_start, xend = times_end, y=-0.2, yend=-0.2, size="Main Effect Stimulus Type"), colour = "#ff8383", linewidth = 1, inherit.aes=FALSE) +
      # geom_segment(data = clusters_F_main1, aes(x=times_start, xend = times_end, y=-0.23, yend=-0.23, size="Main Effect Conditioning"), colour = "#e874ff", linewidth = 1, inherit.aes=FALSE) +
      # geom_segment(data = clusters_F_int, aes(x=times_start, xend = times_end, y=-0.26, yend=-0.26, size="Stimulus Type x Conditioning Interaction "), colour = "#ffdd74", linewidth = 1, inherit.aes=FALSE) +
      scale_x_continuous("Time [s]",limits=c(-0.5, 11), minor_breaks=c(0,1,2,3,4,5,6,7,8,9,10), breaks=c(0, 2, 4, 6, 8, 10)) +
      scale_y_continuous("∆ Skin Conductance [μS]") + #, breaks=c(-80,-40, 0, 40), minor_breaks=c(-80, -60, -40, -20, 0, 20, 40, 60)) +
      scale_color_viridis_d(aesthetics = c("colour", "fill"))
      # theme_bw() +
      # scale_size_manual("effects", values=rep(1,4), guide=guide_legend(override.aes = list(colour=c("#ff8383", "#e874ff", "#ffdd74"))))
  }

# ggsave(file.path("Study 1", "Plots", "EDA", "cs_test_cluster.png"),type="cairo-png", width=2500/400, height=1080/300, dpi=300)

################################################################################
# HR ---------------------------------------------------------------------------
################################################################################
cat("\n", "Cluster-based permutation tests of HR data", "\n")
sample_rate = 2

# Prepare data -----------------------------------------------------------------

heart.wide <- readRDS(file.path("Study 1", "Physio", "HR.RData"))

step_plotting = 0.5
scaling.window = scaling.window = c(seq(-4, 12, by=step_plotting))
baselineWindow = c(-.5, 0) #correct for Baseline in this time window

# CLUSTER
hr_long <- heart.wide %>%
  pivot_longer(hr.1:hr.32, names_to = "samplepoint", values_to ="HR", names_prefix="hr.") %>%
  mutate(samplepoint = as.numeric(samplepoint)) %>%
  filter(condition %in% c(6,7,8,9)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>%
  # mutate(ID = subject) %>% 
  select(subject, trial, condition, outcome, samplepoint, HR) %>%
  mutate(time = ((samplepoint - 1) * step_plotting) + min(scaling.window)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>% 
  filter(time >= 0, time < 10)

heart = heart.wide %>% gather(key="time", value="HRchange", matches("hr\\.\\d+")) %>% tibble() %>% 
  mutate(time = time %>% gsub("hr.", "", .) %>% as.integer() %>% {. * step_plotting + min(scaling.window)} %>% round(2), # %>% as.factor(),
         condition = as.factor(condition)) %>% 
  select(-contains("hr.bin."))

# Calculate null distributions -------------------------------------------------

# with parallelization:
# Use parallelization to run the permutations. As this was written for Windows, the doSNOW backend is used in the current implementation.
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, list("pt_null_distribution", "pt_critical_F", "pt_Ftest_statistic", "perform_anova", "perform_lmm"))
registerDoSNOW(cl)

hr_null_dist <- data.frame()

tmp <- foreach(i = 1:length(cl), .combine = "rbind", .packages = c("tidyverse", "afex", "lme4")) %dopar% {
  y <- pt_null_distribution(hr_long, dv = "HR", within = "condition", factor1 = "condition_threat", factor2 = "condition_social",
                            time = "samplepoint", trial = "trial", id = "subject", nperm = as.integer(ceiling(iterations/length(cl))))
  return(y)
}

stopCluster(cl)
hr_null_dist <- tmp[1:iterations,]

# # without parallelization:
# hr_null_dist <- pt_null_distribution(hr_long, dv = "HR", within = "condition", factor1 = "condition_threat", factor2 = "condition_social",
#                                                time = "time", trial = "trial", id = "subject", nperm = iterations)

# Create a tibble with the 95% quantile of cluster lengths under H0
hr_critical_length_null_F_main1 <- quantile(hr_null_dist$F_main1, .95)
hr_critical_length_null_F_main2 <- quantile(hr_null_dist$F_main2, .95)
hr_critical_length_null_F_int <- quantile(hr_null_dist$F_int, .95)

hr_critical_lengths_null <- data.frame(hr_critical_length_null_F_main1, hr_critical_length_null_F_main2, hr_critical_length_null_F_int)
colnames(hr_critical_lengths_null) <- c("F_main1", "F_main2", "F_int")

# Save distributions and stop cluster ----------------------------------------
save(hr_null_dist, hr_critical_lengths_null, file = file.path("Study 1", "Physio", "hr_null_distributions.RData"))

# Cluster-based permutation tests ---------------------------------------------
load(file.path("Study 1", "Physio", "hr_null_distributions.RData"))

# Find clusters
F_tests <- pt_Ftest_statistic(hr_long, dv = "HR", id = "subject", factor1 = "condition_threat", factor2 = "condition_social", time = "samplepoint")
criticalFs <- pt_critical_F(hr_long, id = "subject", factor1 = "condition_threat", factor2 = "condition_social")
critical_F_main1 <- criticalFs[[1]] 
critical_F_main2 <- criticalFs[[2]]
critical_F_int <- criticalFs[[3]]

clusters_F_main1 <- clusters_over_time(F_tests$F_main1, critical_F_main1)
clusters_F_main1 <- clusters_F_main1 %>% 
  mutate(critical_length = hr_critical_lengths_null$F_main1) %>% 
  mutate(p = 1 - pempiricalD(length, critical_length)) %>% 
  filter(p < .05) %>% 
  mutate(times_start = start / sample_rate - 0.5,
         times_end = end / sample_rate - 0.5)

clusters_F_main2 <- clusters_over_time(F_tests$F_main2, critical_F_main2)
clusters_F_main2 <- clusters_F_main2 %>% 
  mutate(critical_length = hr_critical_lengths_null$F_main2) %>% 
  mutate(p = 1 - pempiricalD(length, critical_length)) %>% 
  filter(p < .05) %>% 
  mutate(times_start = start / sample_rate - 0.5,
         times_end = end / sample_rate - 0.5)

clusters_F_int <- clusters_over_time(F_tests$F_int, critical_F_int)
clusters_F_int <- clusters_F_int %>% 
  mutate(critical_length = hr_critical_lengths_null$F_int) %>% 
  mutate(p = 1 - pempiricalD(length, critical_length)) %>% 
  filter(p < .05) %>% 
  mutate(times_start = start / sample_rate - 0.5,
         times_end = end / sample_rate - 0.5)

save(clusters_F_main1, clusters_F_main2, clusters_F_int, file = file.path("Study 1", "Physio", "hr_cluster.RData"))

load(file.path("Study 1", "Physio", "hr_cluster.RData"))

plot_hr_cluster <- heart %>% 
  filter(condition %in% c(6,7,8,9)) %>%
  filter(time <= 10) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>%
  group_by(condition, time) %>% 
  summarise(HR.se = sd(HRchange, na.rm=T)/sqrt(n()), HR.mean = mean(HRchange, na.rm=T)) %>% 
  mutate(time = time %>% as.character() %>% as.numeric()) %>% 
  {ggplot(., aes(x=time, y=HR.mean)) +
      geom_vline(xintercept=0, color="black",linetype="solid") + #zero = picture onset
      geom_vline(xintercept=10, color="black",linetype="solid") + #picture offset
      geom_line(aes(colour=condition)) +
      geom_ribbon(aes(ymin=HR.mean-HR.se, ymax=HR.mean+HR.se, colour=condition, fill=condition), color = NA, alpha=.2) +
      # geom_segment(data = clusters_F_main2, aes(x=times_start, xend = times_end, y=-2.7, yend=-2.7, size="Main Effect Stimulus Type"), colour = "#ff8383", linewidth = 1, inherit.aes=FALSE) +
      # geom_segment(data = clusters_F_main1, aes(x=times_start, xend = times_end, y=-3, yend=-3, size="Main Effect Conditioning"), colour = "#e874ff", linewidth = 1, inherit.aes=FALSE) +
      # geom_segment(data = clusters_F_int, aes(x=times_start, xend = times_end, y=-3.3, yend=-3.3, size="Stimulus Type x Conditioning Interaction "), colour = "#ffdd74", linewidth = 1, inherit.aes=FALSE) +
      scale_x_continuous("Time [s]",limits=c(-1, 11), minor_breaks=c(0,1,2,3,4,5,6,7,8,9,10), breaks=c(0, 2, 4, 6, 8, 10)) +
      scale_y_continuous("∆ Heart Rate [BPM]") + #, breaks=c(-80,-40, 0, 40), minor_breaks=c(-80, -60, -40, -20, 0, 20, 40, 60)) +
      scale_color_viridis_d(aesthetics = c("colour", "fill")) +
      # theme_bw() +
      scale_size_manual("effects", values=rep(1,4), guide=guide_legend(override.aes = list(colour=c("#ff8383", "#e874ff", "#ffdd74"))))
  }

# ggsave(file.path("Study 1", "Plots", "HR", "cs_test_cluster.png"), type="cairo-png", width=2500/400, height=1080/300, dpi=300)
