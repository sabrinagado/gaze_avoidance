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

{ #load behavioral data (logs)
  path.logs = file.path("Study 1", "Experiment", "gaze_avoidance_task", "data")
  files.log.prefix = "gca_avoidance_task_"
  files.log.extension = ".log"
  
  files.rating.prefix = "gca_avoidance_task_"
  files.rating.extension = ".csv"
  
  # Get Scores
  path.scores = file.path("Study 1", "Questionnaires", "demo_scores.csv")
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
    
    cat("\n", "Read in subject csv files for ratings:", "\n")
    pb <- txtProgressBar(min = 1, max = length(files.rating.path), style = 3)
    i = 1
    for (file.rating.path in files.rating.path) {
      # file.rating.path = files.rating[1]
      rating = suppressMessages(read_csv(file.rating.path, show_col_types=F))
      
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
      
      setTxtProgressBar(pb, i)
      i = i+1
    }
    
    close(pb)
    
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

write.csv2(ratings, file.path("Study 1", "Ratings", "ratings.csv"), row.names=FALSE, quote=FALSE)

for (p in c("Baseline", "Acquisition", "Test")) {
  # p = "Baseline"
  cat("\n\n", "Analyses of the ratings in the", tolower(p), "phase\n")
  
  ratings.phase <- ratings %>%
    filter(phase == tolower(p))

  ratings.phase.summary <- ratings.phase %>%
    summarise(Mean = mean(rating), SD = sd(rating), .by=condition)
  if (p == "Baseline") {y_label = "Rating"} else {y_label = NULL}
  
  anova <- ratings.phase %>%
    # group_by(subject, condition_social, condition_threat) %>%
    mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
    ez::ezANOVA(dv=.(rating),
                wid=.(subject),
                within=.(condition_threat, condition_social),
                # between=.(SPAI),
                detailed=T,
                type=3)
  anova %>% apa::anova_apa()
  print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
  print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
  print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # interaction
  
  cat("\n\n", "Correlation between ratings and SPAI\n")
  print(cor.test(ratings.phase$rating, ratings.phase$SPAI) %>% apa::cor_apa())
}

ratings <- ratings %>% mutate(phase = str_to_title(phase))
ratings <- ratings %>% mutate(phase = factor(phase, levels=c("Baseline", "Acquisition", "Test"), ordered=TRUE))

ratings <- ratings %>% mutate(condition = str_replace_all(condition, "\n", " "))

ratings.summary <- ratings %>%
  summarise(mean_rating = mean(rating), sd_rating = sd(rating), .by=c(phase, condition))

plot_ratings_exp1 <- ggplot(ratings.summary, aes(x = phase, y = mean_rating, group = condition, color=condition)) +
  geom_point(data = ratings, aes(x = phase, y = rating, group = condition, color = condition), 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7), alpha = 0.2, size=1) +
  geom_line(position = position_dodge(width = 0.7), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.7), shape = "square", size=3) +
  geom_errorbar(
    aes(ymin = mean_rating - sd_rating, ymax = mean_rating + sd_rating),
    position = position_dodge(width = 0.7),
    width = 0.25, linewidth = 0.5) +
  labs(y = "Rating", x = "") +
  guides(color = guide_legend(title = ""), size="none") +
  theme_minimal() +
  # theme(legend.position="bottom") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()