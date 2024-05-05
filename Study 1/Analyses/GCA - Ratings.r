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
  requirePackage("cowplot", load=T)
}

{ # Variables ---------------------------------------------------------------
  # a priori exclusions for all variables
  exclusions = c() %>% unique() %>% sort()
  
  vps = vps <- seq(1, 14)
  exclusions.eye.num = c() %>% c(exclusions) #a priori exclusions, e.g. calibration not successful
  # numToEye = function(nums) return(nums %>% formatC(width=2, format="d", flag="0") %>% paste0("gca_", .)) #add leading zeros and prefix "vp"
  # exclusions.eye = exclusions.eye.num %>% numToEye()
  vps <- vps[!(vps %in% exclusions.eye.num)]
}

{ # Paths -------------------------------------------------------------------
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd('..')
  path = getwd()
  
  #load behavioral data (logs)
  path.logs = file.path(path, "Experiment", "gaze_avoidance_task", "data")
  files.log.prefix = "gca_avoidance_task_"
  files.log.extension = ".log"
  
  files.rating.prefix = "gca_avoidance_task_"
  files.rating.extension = ".csv"
  
  # Get Scores
  path.scores = file.path(path, "Questionnaires", "demo_scores.csv")
  scores = read_delim(path.scores, delim=";", locale=locale(decimal_mark=","), na=".", show_col_types=F)
  scores$subject <- scores$VP
  scores <- scores %>%
    select(subject, gender, age, digitimer, temperature, humidity, SPAI, SIAS, STAI_T, UI, motivation, tiredness)
  
  # Files
  files.rating = list.files(path.logs, pattern=paste0("^", files.rating.prefix, ".*", files.rating.extension, "$"), full.names=TRUE)
}

{ # Functions ---------------------------------------------------------------
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
        select(stimtype, sliderStim.response)
      rating <- na.omit(rating)
      rating <- rating %>%
        group_by(stimtype) %>%
        mutate(timepoint = row_number()) %>%
        ungroup()
      rating <- rating %>%
        mutate(phase = "baseline") %>%
        mutate(phase = ifelse(timepoint==2, "acquisition", phase)) %>% 
        mutate(phase = ifelse(timepoint==3, "test", phase))
      
      rating <- rating %>% 
        mutate(subject = file.rating.path %>% sub(".*gca_avoidance_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer())
      rating <- rating %>%
        select(subject, phase, stimtype, sliderStim.response)
      names(rating) = c("subject","phase","condition", "rating")
      
      ratings <- rbind(rating, ratings)
    }
    
    ratings <- ratings %>% arrange(subject)
    
    return(ratings)
  }
}


rep_str = c('cs_minus_ns'='CSpos,\nnon-social','cs_minus_s'='CSpos,\nsocial',
            'cs_plus_ns'='CSneg,\nnon-social', 'cs_plus_s'='CSneg,\nsocial')

###############################################################################
# Ratings
###############################################################################
ratings = loadRatings(files.rating, files.log)  %>% # ratings
  mutate(across('condition', str_replace_all, rep_str))

# ratings.wide <- ratings %>% 
#   pivot_wider(names_from = condition, values_from = rating)
# 
# ratings.wide <- ratings.wide %>% 
#   left_join(scores, by="subject")
# 
# write.csv2(ratings.wide, file.path(path, "ratings_wide.csv"), row.names=FALSE, quote=FALSE)

ratings <- ratings %>%
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>%
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

ratings <- ratings %>%
  left_join(scores, by="subject")

write.csv2(ratings, file.path(path, "Ratings", "ratings.csv"), row.names=FALSE, quote=FALSE)

plots <- list()
i = 1
for (p in c("Baseline", "Acquisition", "Test")) {
  # p = "Baseline"
  ratings.phase <- ratings %>%
    filter(phase == tolower(p))

  ratings.phase.summary <- ratings.phase %>%
    summarise(Mean = mean(rating), SD = sd(rating), .by=condition)
  if (p == "Baseline") {
    y_label = "Rating"
  } else {
    y_label = NULL
  }

  plot <- ggplot(ratings.phase.summary, aes(x = condition, y = Mean, fill = condition)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_errorbar(
      aes(ymin = Mean - SD, ymax = Mean + SD),
      position = position_dodge(width = 0.7),
      width = 0.25) +
    geom_line(data=ratings.phase, aes(y = rating, group = subject), alpha=0.1) +
    geom_point(data=ratings.phase, aes(y = rating), size = 2, shape = 21, color = "black", alpha=0.3) + #, position=position_jitter(width=0.05)) +
    # labs(title = paste("Rating ", p, " Phase (N = ", n_distinct(ratings.phase$subject), ")", sep=""), x = "Conditions", y = "Rating") +
    labs(title = paste(p, "Phase", sep=" "), x = NULL, y = y_label) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_viridis_d() +
    scale_color_viridis_d()

  # ggsave(file.path(path, "Plots", "Ratings", paste0("ratings_", tolower(p), ".png")), width=1500, height=2000, units="px")
  plots[[i]] <- plot

  print(p)
  ratings.phase %>%
    group_by(subject, condition_social, condition_threat) %>%
    mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
    ez::ezANOVA(dv=.(rating),
                wid=.(subject),
                within=.(condition_social, condition_threat),
                between=.(SPAI),
                detailed=T, type=3) %>%
    apa::anova_apa()
  
  print("t-test threat")
  ratings.phase %>%
    summarize(rating = mean(rating), .by=c(subject, condition_threat)) %>% 
    t_test(data=., rating ~ condition_threat, paired = T, alternative = "two.sided") %>% 
    apa::t_apa()
  
  print("t-test stimulus")
  ratings.phase %>%
    summarize(rating = mean(rating), .by=c(subject, condition_social)) %>% 
    t_test(data=., rating ~ condition_social, paired = T, alternative = "two.sided") %>% 
    apa::t_apa()
  
  print(cor(ratings.phase$rating, ratings.phase$SPAI))
  i = i + 1
}

rating_plot <- plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol = 3)
ggsave(file.path(path, "Plots", "Ratings", paste0("Ratings.png")), width=3500, height=1200, units="px")
