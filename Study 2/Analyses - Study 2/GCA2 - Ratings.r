###############################################################################
# Gaze Contingent Avoidance Project
# Sabrina Gado & Yannik Stegmann
# Code adapted from Mario Reutter & Janna Teigeler
###############################################################################

if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)

requirePackage = function(name, load=T) {
  package = as.character(name)
  if (package %in% rownames(installed.packages()) == FALSE) install.packages(package)
  if (load) library(package, character.only=TRUE)
}

{ #install packages needed ---------------------------------------------------------------
  requirePackage("scales", load=F)
  requirePackage("cowplot", load=F)
  requirePackage("readxl", load=F)
  requirePackage("apaTables", load=F)
  requirePackage("schoRsch", load=F)
  requirePackage("apa", load=F)
  requirePackage("ez", load=F)
  requirePackage("TOSTER", load=F)
  requirePackage("psych", load=F)
  requirePackage("DescTools")
  requirePackage("stringr", load=T)
  requirePackage("emmeans", load=T)
  requirePackage("cowplot")
}

{ # Variables ---------------------------------------------------------------
  # a priori exclusions for all variables
  exclusions = c(100) %>% unique() %>% sort()
  
  acq1End = 32
  acq2End = acq1End + 32
  test1End = acq2End + 48
  test2End = test1End + 48
  # trials.n = 112 # number of trials that shall be analyzed (if more trials, last ones will be taken)
}

{ # Paths -------------------------------------------------------------------
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd('..')
  path = getwd()
  
  #load behavioral data (logs)
  path.logs = file.path(path, "Experiment", "attentional_competition_task", "data")
  files.log.prefix = "attentional_competition_task_"
  files.log.extension = ".log"
  
  path.conds = file.path(path, "Experiment", "attentional_competition_task", "data")
  files.cond.prefix = "attentional_competition_task_"
  files.cond.extension = ".csv"
  
  files.rating.prefix = "attentional_competition_task_"
  files.rating.extension = ".csv"
  
  
  # Get Scores
  path.scores = file.path(path, "Questionnaires", "demo_scores.csv")
  scores = read_delim(path.scores, delim=";", locale=locale(decimal_mark=","), na=".", show_col_types=F)
  scores$subject <- scores$VP
  scores <- scores %>%
    select(subject, gender, age, motivation_points, SPAI, SIAS, STAI_T, UI, motivation, tiredness) %>% 
    mutate(motivation_points = as.numeric(motivation_points))
  
  # Files
  files.rating = list.files(path.logs, pattern=paste0("^", files.rating.prefix, ".*", files.rating.extension, "$"), full.names=TRUE)
  files.rating <- files.rating[-c(exclusions)]
  
  files.log = list.files(path.logs, pattern=paste0("^", files.log.prefix, ".*", files.log.extension, "$"), full.names=TRUE)
  files.log <- files.log[-c(exclusions)]
  
  files.cond = list.files(path.conds, pattern=paste0("^", files.cond.prefix, ".*", files.cond.extension, "$"), full.names=TRUE)
  files.cond <- files.cond[-c(exclusions)]
}

{ # Functions ---------------------------------------------------------------
  loadConditions = function(files.cond.path) {
    conds <- data.frame(
      subject = integer(0),
      trial = integer(0),
      phase = character(0),
      condition = character(0),
      score = integer(0),
      outcome = character()
    )
    for (file.cond.path in files.cond.path) {
      # file.cond.path = files.cond.path[1]
      cond = read_csv(file.cond.path)
      cond <- cond %>% drop_na("score")
      cond <- cond %>%
        select(stim, stim1, stim2, stim3, stim4, score, outcome)
      cond <- cond %>% 
        mutate(condition = ifelse(is.na(stim), "all", stim),
               subject = file.cond.path %>% sub(".*attentional_competition_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer())
      cond <- cond %>% 
        mutate(new_trial = ifelse((stim1 != lag(stim1)) | (stim2 != lag(stim2)) | (stim3 != lag(stim3)) | (stim4 != lag(stim4)), TRUE, FALSE)) %>% 
        mutate(new_trial = ifelse(is.na(new_trial), TRUE, new_trial))
      
      if (sum(cond$new_trial) < 160) {
        cond$new_trial = TRUE
      }
      cond <- cond %>% 
        filter(new_trial)
      cond <- cond %>% 
        mutate(trial = seq(1, nrow(cond)))
      cond <- cond %>% 
        mutate(outcome = ifelse(trial <= acq2End, outcome, NA),
               phase = ifelse(trial <= acq2End, "acquisition", "test"))
      cond <- cond %>%
        select(subject, trial, phase, condition, score, outcome)
      
      conds <- rbind(conds, cond)
    }
    
    conds <- conds %>% arrange(subject, trial)
    
    return(conds)
  }
  
  loadRatings = function(files.rating.path, files.log.path) {
    ratings <- data.frame(
      subject = integer(0),
      phase = character(0),
      condition = character(0),
      rating = numeric(0)
    )
    for (file.rating.path in files.rating.path) {
      # file.rating.path = files.rating[1]
      rating = read_csv(file.rating.path)
      rating <- rating %>%
        select(stim, sliderStim.response)
      rating <- na.omit(rating)
      rating <- rating %>%
        group_by(stim) %>%
        mutate(timepoint = row_number()) %>%
        ungroup()
      rating <- rating %>%
        mutate(timepoint = ifelse(grepl("new", stim), 3, timepoint))
      rating <- rating %>%
        mutate(phase = "baseline") %>%
        mutate(phase = ifelse(timepoint==2, "acquisition", phase)) %>% 
        mutate(phase = ifelse(timepoint==3, "test", phase))
      
      rating <- rating %>% 
        mutate(subject = file.rating.path %>% sub(".*attentional_competition_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer())
      rating <- rating %>%
        select(subject, phase, stim, sliderStim.response)
      names(rating) = c("subject","phase","condition", "rating")
      
      ratings <- rbind(rating, ratings)
    }
    
    ratings <- ratings %>% arrange(subject)
    
    return(ratings)
  }
  
  partial_eta_squared_ci <- function(x)
  {
    # Calculate upper and lower limit of Lambda for 90%-CI
    limits <- MBESS::conf.limits.ncf(x$F, .95, x$DFn, x$DFd)
    
    # Convert Lambda to partial eta-squared
    ci <- c(limits$Lower.Limit, limits$Upper.Limit) /
      (c(limits$Lower.Limit, limits$Upper.Limit) + x$DFn + x$DFd + 1)
    
    # Replace NA with zero (`conf.limits.ncf` returns NA instead of 0)
    ci[is.na(ci)] <- 0
    
    ci <- sprintf("%.2f", round(ci, 2))
    
    paste0("[", ci[1], "; ", ci[2], "]")
  }
}

rep_str = c('cs_pos_ns_new'='CSpos,\nnon-social,\nnew','cs_pos_s_new'='CSpos,\nsocial,\nnew',
            'cs_neg_ns_new'='CSneg,\nnon-social,\nnew', 'cs_neg_s_new'='CSneg,\nsocial,\nnew',
            'cs_pos_ns'='CSpos,\nnon-social','cs_pos_s'='CSpos,\nsocial',
            'cs_neg_ns'='CSneg,\nnon-social', 'cs_neg_s'='CSneg,\nsocial')

###############################################################################
# Ratings
###############################################################################
ratings = loadRatings(files.rating, files.log)  %>% # ratings
  mutate(across('condition', str_replace_all, rep_str))

# ratings.wide <- ratings %>%
#   pivot_wider(names_from = condition, values_from = rating)
# 
# ratings.wide <- ratings.wide %>%
#   left_join(scores, by="subject") %>%
#   rename_with(~ gsub(",\n", "_", .x, fixed = TRUE))
# 
# write.csv2(ratings.wide, file.path(path, "ratings_wide.csv"), row.names=FALSE, quote=FALSE)

ratings <- ratings %>%
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>%
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>%
  mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>% 
  mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))

ratings <- ratings %>%
  left_join(scores, by="subject")

write.csv2(ratings, file.path(path, "Ratings", "ratings.csv"), row.names=FALSE, quote=FALSE)

# # Filter for Motivation
# ratings <- ratings %>% filter(motivation_points > 4)


plots <- list()
i = 1

for (p in c("Baseline", "Acquisition", "Test")) {
  # p = "Test"
  ratings.phase <- ratings %>%
    filter(phase == tolower(p))
  
  if (p == "Baseline") {
    y_label = "Rating"
  } else {
    y_label = NULL
  }

  ratings.phase.summary <- ratings.phase %>%
    summarise(Mean = mean(rating), SD = sd(rating), .by=c(condition, condition_social, condition_threat, condition_novelty)) %>% 
    mutate(condition = str_remove(condition, ",\nnew"))
  ratings.phase <- ratings.phase %>% 
    mutate(condition = str_remove(condition, ",\nnew"))

  plot <- ggplot(ratings.phase.summary, aes(x = condition, y = Mean, fill = condition)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_errorbar(
      aes(ymin = Mean - SD, ymax = Mean + SD),
      position = position_dodge(width = 0.7),
      width = 0.25) +
    #geom_line(data=ratings.phase, aes(y = rating, group = subject), alpha=0.1) +
    geom_point(data=ratings.phase, aes(y = rating), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
    # labs(title = paste("Rating ", p, " Phase (N = ", n_distinct(ratings.phase$subject), ")", sep=""), x = "Conditions", y = "Rating") +
    labs(title = paste(p, "Phase", sep=" "), x = NULL, y = y_label) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_viridis_d() +
    scale_color_viridis_d() + 
    coord_cartesian(ylim = c(1, 10)) + 
    scale_y_continuous(limits=c(0, 10), breaks=c(2, 4, 6, 8, 10), minor_breaks=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
  
  if(p=="Test") {
    plot <- plot + facet_grid(cols = vars(condition_novelty))
    ggsave(file.path(path, "Plots", "Ratings", paste0("ratings_", tolower(p), ".png")), width=2800, height=2000, units="px")
    anova <- ratings.phase %>%
      group_by(subject, condition_social, condition_threat, condition_novelty) %>%
      mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat), condition_novelty = as.factor(condition_novelty)) %>%
      ez::ezANOVA(dv=.(rating),
                  wid=.(subject),
                  within=.(condition_threat, condition_social, condition_novelty),
                  # between=.(SPAI),
                  detailed=T, type=3)
    anova %>% apa::anova_apa()
    print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
    print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
    print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # novelty
    print(anova$ANOVA[5,] %>% partial_eta_squared_ci()) # threat * social
    print(anova$ANOVA[6,] %>% partial_eta_squared_ci()) # threat * novelty
    print(anova$ANOVA[7,] %>% partial_eta_squared_ci()) # social * novelty
    print(anova$ANOVA[8,] %>% partial_eta_squared_ci()) # threat * social * novelty
    print(cor(ratings.phase$rating, ratings.phase$SPAI))
    
  } else {
    ggsave(file.path(path, "Plots", "Ratings", paste0("ratings_", tolower(p), ".png")), width=1500, height=2000, units="px")
    anova <- ratings.phase %>%
      group_by(subject, condition_social, condition_threat) %>%
      mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
      ez::ezANOVA(dv=.(rating),
                  wid=.(subject),
                  within=.(condition_threat, condition_social),
                  # between=.(SPAI),
                  detailed=T, type=3)
    anova %>% apa::anova_apa()
    print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
    print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
    print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # interaction
    print(cor(ratings.phase$rating, ratings.phase$SPAI))
  }
  
  plots[[i]] <- plot
  i = i+1
}

rating_plot <- plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 3, rel_widths = c(1, 1, 2))
ggsave(file.path(path, "Plots", "Ratings", paste0("Ratings.png")), width=3500, height=1200, units="px")
