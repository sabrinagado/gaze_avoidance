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
  requirePackage("lme4", load=T)
  # requirePackage("lmerTest", load=T)
}

{ # Variables ---------------------------------------------------------------
  # a priori exclusions for all variables
  exclusions = c(100) %>% unique() %>% sort()
  
  acq1End = 32
  acq2End = acq1End + 32
  test1End = acq2End + 48
  test2End = test1End + 48
  # trials.n = 112 # number of trials that shall be analyzed (if more trials, last ones will be taken)
  
  sample.rate = 1000 #samples/second
  
  screen.height = 1080 #height of screen in pix
  screen.width  = 1920 # width of screen in pix
  
  image.acq.height = 450 # height of picture in acq phase in px
  image.acq.width = 426 # width of picture in acq phase in px
  
  exclusions.eye.num = c() %>% c(exclusions) #a priori exclusions, e.g. calibration not successful
  
  #baseline validation
  baseline = c(700, 1000) #Baseline in ms relative to trial onset; min(baseline) = start; max(baseline) = end
  # cross_offset = 1000
  # useAllBaselines = list("83" = 2, "84" = 3) #manually allow all baselines of vp83 phase 2 (Gen1) & vp84 phase 3 (Gen2)
  saveBaselinePlots = TRUE
  driftPlots = T #c("vp30", "vp33")
  maxDeviation_rel = 3 #max abs value of z-score
  outlierLimit.eye = .5 #maximum percentage of invalid baselines per subject
  maxSpread = 150 #maximum spread of valid baselines (diameter / edge length)
  usePointDistance = T #use point (vector) distance instead of coordinates independently?
  
  # #ROI analysis
  validFixTime.trial = .5 #percentage of valid fixation time within trial in order to analyze trial
  validFixTime.subj = .5 #percentage of trials with sufficient valid fixation time in order to analyze subject
  # diagnosticDwell = .5 #percentage of trials per subject that need at least one fixation towards the diagnostic ROI
  # showTrialPlots = F
  # unifyRois = T
  # binResolution = 500 #resolution of temporal bins in ms
  # bins = seq(0, trialEnd, binResolution) #bins for temporal dynamic analysis (in ms)
  
  # screen: 24" ASUS VG248QE
  screen.width <- 1920
  screen.height <- 1080
  screen.width.cm <- 53.136
  screen.height.cm <- 29.889
  pixsize_in_cm <- screen.width.cm / screen.width  # Pixel size in cm
  distance <- 500 # Distance Camera - Eye 500 mm (chin rest with remote tracking)
  
  # #statistical analysis
  # z.max = 2 #winsorize dependent variables to a z value of 2
  # q.max = pnorm(z.max * c(-1, 1)) #needed for DescTools::Winsorize function
  # q.max = c(0, 1) #switch off Winsorizing (comment out to switch on)
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
  
  path.eye = file.path(path, "Experiment", "attentional_competition_task", "data") #eye tracking data
  path.pupil = file.path(path, "Experiment", "attentional_competition_task", "data") #pupil data
  path.plots = file.path(path, "Plots", "Gaze")
  
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
  outlier_remove <- function(x, z=3) {
    # Remove outlier iteratively
    insample = 1 - is.na(x) #insample <- rep(1,length(x)); insample[is.na(x)] <- 0
    ok <- FALSE
    while (!ok && sum(insample) > 3) {
      xminpos <- (1:length(x))[x==min(x[insample==1], na.rm=T)] #can contain NAs
      xminpos <- xminpos[xminpos %in% (1:length(insample))[insample==1]][1] #eliminates NAs (and takes only first minimum if several are present)
      xmaxpos <- (1:length(x))[x==max(x[insample==1], na.rm=T)] #can contain NAs
      xmaxpos <- xmaxpos[xmaxpos %in% (1:length(insample))[insample==1]][1] #eliminates NAs (and takes only first maximum if several are present)
      tempinsample <- insample; tempinsample[c(xminpos,xmaxpos)] <- 0
      subx <- x[tempinsample==1]
      
      if (x[xminpos] < (mean(subx) - z*sd(subx))) {
        insample[xminpos] <- 0
        out1 <- TRUE
      } else {
        out1 <- FALSE
      }
      
      if (x[xmaxpos] > (mean(subx) + z*sd(subx))) {
        insample[xmaxpos] <- 0
        out2 <- TRUE
      } else {
        out2 <- FALSE
      }
      
      if (!out1 & !out2) { ok <- TRUE }
    }
    return(insample==1)
  }
  
  #removes values until (absolute) deviation is satisfied
  outlier_remove_spread = function(x, deviation, insample=NULL) {
    if (is.null(insample)) insample = T %>% rep(length(x))
    insample[is.na(x)] = F
    
    spread = x[insample] %>% range() %>% diff()
    
    while (spread > deviation) {
      m = x[insample] %>% mean()
      x[!insample] = m #set outliers to m to ignore them in next step
      out.next = {x - m} %>% abs() %>% which.max()
      
      insample[out.next] = F #crucial: use x, not x[insample] to match vector length of insample
      spread = x[insample] %>% range() %>% diff()
    }
    return(insample)
  }
    
  
  #TODO after outlier removal is finished, look once (?) more if "outliers" can be INCLUDED (can happen for biased distributions)
  outlier_remove_abs <- function(x, deviation) {
    # Remove outlier iteratively
    insample = 1 - is.na(x) #insample <- rep(1,length(x)); insample[is.na(x)] <- 0
    ok <- FALSE
    while (!ok && sum(insample) > 2) {
      xminpos <- (1:length(x))[x==min(x[insample==1], na.rm=T)]
      xminpos <- xminpos[xminpos %in% (1:length(insample))[insample==1]][1]
      xmaxpos <- (1:length(x))[x==max(x[insample==1], na.rm=T)]
      xmaxpos <- xmaxpos[xmaxpos %in% (1:length(insample))[insample==1]][1]
      tempinsample <- insample; tempinsample[c(xminpos,xmaxpos)] <- 0
      subx <- x[tempinsample==1]
      
      if (x[xminpos] < (mean(subx) - deviation)) {
        insample[xminpos] <- 0
        out1 <- TRUE
      } else {
        out1 <- FALSE
      }
      
      if (x[xmaxpos] > (mean(subx) + deviation)) {
        insample[xmaxpos] <- 0
        out2 <- TRUE
      } else {
        out2 <- FALSE
      }
      
      if (!out1 & !out2) { ok <- TRUE }
    }
    return(insample==1)
  }
  
  outlier_remove_point = function(x, y, z=3, insample=NULL) {
    if (length(x) != length(y)) warning("outlier_remove_point: x- and y-coordinates of unequal length.")
    if (is.null(insample)) insample = T %>% rep(min(c(length(x), length(y))))
    insample[is.na(x) | is.na(y)] = F
    
    ok = F
    while (!ok && sum(insample) > 3) {
      xmean = mean(x[insample]); ymean = mean(y[insample]) #calculate centroid
      distances = pointDistance(x, y, xmean, ymean) #have to be updated every iteration because centroid shifts
      distances[!insample] = 0 #ignore previous outliers for this iteration
      out.next = distances %>% which.max() #find most extreme point (that has not yet been removed)
      insample[out.next] = F #temporarily remove point
      
      if (distances[out.next] <= mean(distances[insample]) + z*sd(distances[insample])) {
        insample[out.next] = T
        ok = T
      }
    }
    
    return(insample)
  }
  
  outlier_remove_point_spread = function(x, y, deviation, insample=NULL) {
    if (length(x) != length(y)) warning("outlier_remove_point_spread: x- and y-coordinates of unequal length.")
    if (is.null(insample)) insample = T %>% rep(min(c(length(x), length(y))))
    insample[is.na(x) | is.na(y)] = F
    
    ok = F
    while (!ok && sum(insample) > 2) {
      xmean = mean(x[insample]); ymean = mean(y[insample]) #calculate centroid
      distances = pointDistance(x, y, xmean, ymean) #have to be updated every iteration because centroid shifts
      distances[!insample] = 0 #ignore previous outliers for this iteration
      out.next = distances %>% which.max() #find most extreme point (that has not yet been removed)
      insample[out.next] = F #temporarily remove point
      
      if (distances[out.next] <= deviation/2) { #deviation = diameter but with distances only radius can be checked => deviation/2
        insample[out.next] = T
        ok = T
      }
    }
    
    return(insample)
  }
  
  pointInBounds = function(point.x, point.y, x.bounds, y.bounds) {
    if (is.na(point.x) || is.na(point.y) || length(point.x)==0 || length(point.y)==0) return(FALSE)
    return(point.x >= min(x.bounds) & point.x <= max(x.bounds) & 
             point.y >= min(y.bounds) & point.y <= max(y.bounds))
  }
  
  degToCm = function(angle, distance) {
    return(distance*tan(angle*pi/180))
  }
  
  cmToDeg = function(cm, distance) {
    return(tan(cm/distance)*180/pi)
  }
  
  degToPix = function(angle, distance, resolution, screenSize) {
    return(degToCm(angle, distance) * resolution / screenSize)
  }
  
  cmToPix <- function(point_in_cm, screenSize_px, screenSize_cm) {
    return(point_in_cm * screenSize_px / screenSize_cm)
  }
  
  pixToCm <- function(point_in_pixels, screenSize_px, screenSize_cm) {
    return(point_in_pixels * screenSize_cm / screenSize_px)
  }
  
  # Calculate distance to monitor for deflected fixations
  distfix <- function(x,y,distance) {
    distcalc <- sqrt(x^2+y^2+distance^2)
    return(distcalc)
  }
  
  # Calculate angle between 2 points
  # https://www.sr-research.com/eye-tracking-blog/background/visual-angle/
  pointDistance = function(x1, y1, x2, y2) {
    return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
  }
  
  visangle <- function(x1, y1, x2, y2, screen.distance){
    # x1 = 27.79677; x2 = 9.287730; y1 = 14.28584; y2 = 23.695335
    size <- pointDistance(x1, y1, x2, y2)
    theta <- atan(size/screen.distance)
    theta.deg <- theta/pi*180
    as.numeric(theta.deg)
  }
  
  # # Calculate angle between 2 vectors (to compute scanpath)
  # angle <- function(x,y){
  #   dot.prod <- x%*%y 
  #   norm.x <- norm(x,type="2")
  #   norm.y <- norm(y,type="2")
  #   theta <- acos(dot.prod / (norm.x * norm.y))
  #   theta.deg <- theta/pi*180
  #   as.numeric(theta.deg)
  # }
  
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
        select(stim, stim1, stim2, stim3, stim4, score, outcome, learning_trials.thisN, trials.thisN, testtrials.thisN, testtrials_novelty.thisN) %>% 
        mutate(trial = ifelse(!is.na(learning_trials.thisN), learning_trials.thisN + 1,
                              ifelse(!is.na(trials.thisN), trials.thisN + 1 + acq1End,
                                     ifelse(!is.na(testtrials.thisN), testtrials.thisN + 1 + acq2End, testtrials_novelty.thisN + 1 + test1End)))) %>% 
        select(-c(learning_trials.thisN, trials.thisN, testtrials.thisN, testtrials_novelty.thisN))
      
      cond <- cond %>% 
        mutate(condition = ifelse(is.na(stim), "all", stim),
               subject = file.cond.path %>% sub(".*attentional_competition_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer())
      cond <- cond %>% 
        mutate(new_trial = ifelse((trial != lag(trial)), TRUE, FALSE)) %>% 
        mutate(new_trial = ifelse(is.na(new_trial), TRUE, new_trial))
      
      cond <- cond %>% 
        filter(new_trial) %>% 
        select(-new_trial)

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
  
  loadFixations = function(filePath, screen.height) {
    rename_list <- list(
      "subject" = "RECORDING_SESSION_LABEL", 
      "trial" = "TRIAL_INDEX",
      "start" = "CURRENT_FIX_START",
      "end" = "CURRENT_FIX_END",
      "x" = "CURRENT_FIX_X",
      "y" = "CURRENT_FIX_Y",
      "eye" = "EYE_USED"
    )
    
    fixations = read_delim(filePath, delim="\t", col_names=T, locale=locale(decimal_mark=","), na=".", show_col_types=F)
    fixations <- fixations %>% rename(!!!rename_list)
    
    fixations = fixations %>%
      mutate(subject = subject %>% sub("attentional_competition_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer(),
             y = screen.height - y)
    fixations = fixations[c("subject","trial","start","end","x","y", "eye")]
    fixations <- fixations %>% arrange(subject, trial)
    return(fixations)
  }
  
  loadSaccades = function(filePath, screen.height) {
    rename_list <- list(
      "subject" = "RECORDING_SESSION_LABEL", 
      "trial" = "TRIAL_INDEX",
      "contains_blink" = "CURRENT_SAC_CONTAINS_BLINK",
      "start_time" = "CURRENT_SAC_START_TIME",
      "end_time" = "CURRENT_SAC_END_TIME",
      "start_x" = "CURRENT_SAC_START_X",
      "end_x" = "CURRENT_SAC_END_X",
      "start_y" = "CURRENT_SAC_START_Y",
      "end_y" = "CURRENT_SAC_END_Y"
    )
    
    saccades = read_delim(filePath, delim="\t", col_names=T, locale=locale(decimal_mark=","), na=".", show_col_types=F)
    saccades <- saccades %>% rename(!!!rename_list)
    saccades = saccades %>%
      mutate(subject = subject %>% sub("attentional_competition_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer(),
             start_y = screen.height - start_y, end_y = screen.height - end_y)
    saccades = saccades[c("subject","trial","contains_blink", "start_time", "end_time", "start_x", "start_y","end_x","end_y")]
    saccades <- saccades %>% arrange(subject, trial)
    return(saccades)
  }
  
  loadMessages = function(filePath) {
    rename_list <- list(
      "subject" = "RECORDING_SESSION_LABEL", 
      "trial" = "TRIAL_INDEX",
      "time" = "CURRENT_MSG_TIME",
      "event" = "CURRENT_MSG_TEXT"
    )
    
    messages = read_delim(filePath, delim="\t", col_names=T, locale=locale(decimal_mark=","), na=".", show_col_types=F)
    messages <- messages %>% rename(!!!rename_list)
    messages = messages %>%
      mutate(subject = subject %>% sub("attentional_competition_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer())
    messages = messages[c("subject", "trial", "time","event")]
    messages <- messages %>% arrange(subject, trial)
    return(messages)
  }
  
  validateBaselines = function(fixs, mess, exclusions, maxDeviation_rel, maxSpread, saveBaselinePlots=TRUE, postfix="") {
    if (saveBaselinePlots) dir.create(file.path(path.plots, "baseline"), showWarnings=FALSE)
    
    vpn = fixs$subject %>% unique() %>% sort() %>% as.character() #all subjects in fixations
    vpn = vpn[!(vpn %in% exclusions)] #minus a priori exclusions
    vpn.n = length(vpn)
    
    trials.n = fixs %>% group_by(subject) %>% summarise(trials = max(trial) - min(trial) + 1) %>% .$trials %>% max()
    
    baselines.trial <- data.frame(
      subject = integer(0),
      trial = integer(0),
      x = numeric(0),
      y = numeric(0),
      x_divergence = numeric(0),
      y_divergence = numeric(0),
      condition = character(0),
      blok = logical(0)
    )
    
    # baselines.trial = vector("list", length(vpn.n)) #list of evaluations of baselines per trial for every subject
    baselines.summary = data.frame(subject=character(vpn.n), ntrials=numeric(vpn.n), nValid=numeric(vpn.n), 
                                   invalid=numeric(vpn.n), na=numeric(vpn.n), 
                                   sd_x=numeric(vpn.n), sd_y=numeric(vpn.n),
                                   range_x=numeric(vpn.n), range_y=numeric(vpn.n),
                                   stringsAsFactors=FALSE)
    baselines.summary[,] = NA
    
    print("Validating baselines")
    for (vpIndex in 1:vpn.n) { #TODO make inner of loop enclosed function and call it with apply?
      # vpIndex = 1
      code = vpn[vpIndex]
      cat(code, " ... ")
      
      # Determine trial number
      fix.vp = fixs %>% filter(subject==code)
      msg.vp = mess %>% filter(subject==code)
      trials.min = min(fix.vp$trial)
      trials.max = max(fix.vp$trial)
      
      # baseline
      baseline.vp = data.frame(x=numeric(trials.n), y=numeric(trials.n), condition=character(trials.n), stringsAsFactors=FALSE)
      baseline.vp[,] = NA
      
      # Loop over trials to determine trial-by-trial baselines
      for (trial in fix.vp$trial %>% unique()) {
        # trial=1
        # Select trial data
        fix.vp.trial = fix.vp[fix.vp$trial==trial,] #fix.vp %>% filter(trial==trial) #doesn't work
        msg.vp.trial = msg.vp[msg.vp$trial==trial,] #msg.vp %>% filter(trial==trial)
        
        condition = fix.vp.trial$condition %>% unique()
        
        # # Determine onset (in ms)
        # onset = cross_offset # Stop Fixation Cross
        # # Subtract onset from startamps
        # fix.vp.trial$start  <- fix.vp.trial$start - onset
        # fix.vp.trial$end <- fix.vp.trial$end - onset
        
        # Caluculate baseline as weighted average of fixations
        fix.vp.trial.bl <- fix.vp.trial[fix.vp.trial$start<max(baseline) & fix.vp.trial$end>min(baseline),]
        if (nrow(fix.vp.trial.bl) > 0) {
          # Restrict fixation data to baseline
          fix.vp.trial.bl$start[1] <- max(c(min(baseline), first(fix.vp.trial.bl$start)))
          fix.vp.trial.bl$end[nrow(fix.vp.trial.bl)] <- min(c(max(baseline), last(fix.vp.trial.bl$end))) #ifelse(tail(fix.vp.trial.bl$end,1)>blen,blen,tail(fix.vp.trial.bl$end,1)) # = min(blen, tail...) ?
          fix.vp.trial.bl$dur <- fix.vp.trial.bl$end - fix.vp.trial.bl$start
          
          # Calculate baseline coordinates
          xbl <- sum(fix.vp.trial.bl$x*fix.vp.trial.bl$dur)/sum(fix.vp.trial.bl$dur)
          ybl <- sum(fix.vp.trial.bl$y*fix.vp.trial.bl$dur)/sum(fix.vp.trial.bl$dur)
          
          # Store values
          baseline.vp[trial-trials.min+1,] = c(xbl,ybl, condition)
        } else {
          # When no valid fixations are available store NA as baseline for current trial
          baseline.vp[trial-trials.min+1, "condition"] = condition
        }
      }
      
      baseline.vp <- baseline.vp %>% mutate(x = as.numeric(x), y = as.numeric(y))
      baseline.vp <- baseline.vp %>% mutate(x_divergence = screen.width/2 - x, y_divergence = screen.height/2 - y)
      
      # Determine outlier
      if (usePointDistance) {
        baseline.vp = baseline.vp %>% mutate(
          blok = outlier_remove_point(x, y, maxDeviation_rel) %>%
            outlier_remove_point_spread(x, y, maxSpread, insample=.))
      } else {
        blxok = outlier_remove(baseline.vp$x, maxDeviation_rel) %>% outlier_remove_spread(baseline.vp$x, maxSpread, .)
        blyok = outlier_remove(baseline.vp$y, maxDeviation_rel) %>% outlier_remove_spread(baseline.vp$y, maxSpread, .)
        baseline.vp$blok = blxok & blyok #Baseline is valid when x and y coordinates are ok (i.e. no outlier)
      }
      
      # Store number of valid baselines per subject
      nValid = sum(baseline.vp$blok)
      invalid = 1 - nValid / trials.n
      nas = sum(is.na(baseline.vp$x)) / trials.n
      x.coords = baseline.vp %>% filter(blok) %>% .$x; mean_x = mean(x.coords); sd_x = sd(x.coords); range_x = x.coords %>% range() %>% diff()
      y.coords = baseline.vp %>% filter(blok) %>% .$y; mean_y = mean(y.coords); sd_y = sd(y.coords); range_y = y.coords %>% range() %>% diff()
      baselines.summary[vpIndex, 2:9] = c(trials.n, nValid, invalid, nas, sd_x, sd_y, range_x, range_y); baselines.summary[vpIndex , 1] = code
      
      baselines.trial <- rbind(baselines.trial, baseline.vp %>% mutate(subject = as.numeric(code), trial = row_number()))
      
      #plot subject
      if (saveBaselinePlots==TRUE || code %in% saveBaselinePlots || as.character(code) %in% saveBaselinePlots) {
        filename = file.path(path.plots, "baseline", postfix, sprintf("gca_%s_%s.png", code, postfix))
          
        borders.rel.x = c(mean_x - maxDeviation_rel*sd_x, mean_x + maxDeviation_rel*sd_x)
        borders.rel.y = c(mean_y - maxDeviation_rel*sd_y, mean_y + maxDeviation_rel*sd_y)
        borders.abs.x = c(x.coords %>% range() %>% mean() - maxSpread/2, 
                          x.coords %>% range() %>% mean() + maxSpread/2)
        borders.abs.y = c(y.coords %>% range() %>% mean() - maxSpread/2, 
                          y.coords %>% range() %>% mean() + maxSpread/2)
        
        if (driftPlots==T || code %in% driftPlots || as.character(code) %in% driftPlots) {
          blplot = baseline.vp %>% mutate(trial=1:n()) %>% 
            ggplot(aes(x=x, y=y, color=trial)) + 
            xlim(0, screen.width) + ylim(0, screen.height) + #restrict area to screen
            geom_rect(xmin=0, ymin=0, xmax=screen.width, ymax=screen.height, color="black", fill=NA) +
            geom_hline(yintercept=mean_y, linetype="longdash") + geom_vline(xintercept = mean_x, linetype="longdash") #centroid of valid baselines
          
          #add borders depending on whether point distance was used (circles) or not (rectangles)
          if (usePointDistance) {
            distances = pointDistance(x.coords, y.coords, mean_x, mean_y)
            r_abs = mean(distances) + maxDeviation_rel*sd(distances)
            
            blplot = blplot + #circles
              ggforce::geom_circle(aes(x0=mean_x, y0=mean_y, r=r_abs), color="red", fill=NA, inherit.aes=F) + #borders for validity (relative)
              ggforce::geom_circle(aes(x0=mean_x, y0=mean_y, r=maxSpread/2), color="orange", fill=NA, inherit.aes=F) #borders for validity (absolute)
          } else {
            borders.rel.x = c(mean_x - maxDeviation_rel*sd_x, mean_x + maxDeviation_rel*sd_x)
            borders.rel.y = c(mean_y - maxDeviation_rel*sd_y, mean_y + maxDeviation_rel*sd_y)
            borders.abs.x = c(x.coords %>% range() %>% mean() - maxSpread/2, 
                              x.coords %>% range() %>% mean() + maxSpread/2)
            borders.abs.y = c(y.coords %>% range() %>% mean() - maxSpread/2, 
                              y.coords %>% range() %>% mean() + maxSpread/2)
            
            blplot = blplot + #rectangles
              geom_rect(xmin=min(borders.rel.x), ymin=min(borders.rel.y), xmax=max(borders.rel.x), ymax=max(borders.rel.y), color="red", fill=NA) + #borders for validity (relative)
              geom_rect(xmin=min(borders.abs.x), ymin=min(borders.abs.y), xmax=max(borders.abs.x), ymax=max(borders.abs.y), color="orange", fill=NA) #borders for validity (absolute)
          }
          
          #add points
          blplot = blplot + 
            geom_point() + scale_color_continuous(low="blue", high="green") + #all points (color coded by trial)
            geom_point(data=baseline.vp %>% filter(blok==F), mapping=aes(x=x, y=y), color="red") + #invalid baselines red
            geom_point(x=screen.width/2, y=screen.height/2, shape="+", size=5, color="black") + #fixation cross
            #theme(panel.border = element_rect(color = "black", fill=NA, size=5)) +
            coord_fixed() +
            ggtitle(paste0(code, "_", postfix, " (", round(invalid, digits=2)*100, "% out, ", round(nas, digits=2)*100, "% NAs)"))
          blplot %>% ggsave(filename=filename, plot=., device="png", dpi=300, units="in", width=1920/300, height = 1080/300)
        } else {
          borders.rel.x = c(mean_x - maxDeviation_rel*sd_x, mean_x + maxDeviation_rel*sd_x)
          borders.rel.y = c(mean_y - maxDeviation_rel*sd_y, mean_y + maxDeviation_rel*sd_y)
          borders.abs.x = c(x.coords %>% range() %>% mean() - maxSpread/2, 
                            x.coords %>% range() %>% mean() + maxSpread/2)
          borders.abs.y = c(y.coords %>% range() %>% mean() - maxSpread/2, 
                            y.coords %>% range() %>% mean() + maxSpread/2)
          
          #x11() #plot in new windows (max of 63)
          png(filename, 
              width=1920, height=1080)
          plot(baseline.vp$x, baseline.vp$y, pch=16, col=ifelse(baseline.vp$blok==0, "red", "black"), xlab="x (px)", ylab="y (px)", xlim=c(0, screen.width), ylim=c(0, screen.height), asp=1)
          title(paste0(code, "_", postfix, " (", round(invalid, digits=2)*100, "% out, ", round(nas, digits=2)*100, "% NAs)"))
          abline(h=c(mean_y, screen.height), v=mean_x); abline(v=borders.rel.x, h=borders.rel.y, col="red"); points(x=mean(c(0, screen.width)), y=screen.height/2, pch=3, col="blue")
          
          dev.off()
        }
      }
    }
    print("Baseline validation finished")
    result_list <- list(baseline.summary = baselines.summary, baseline.trials = baselines.trial %>% select(subject, trial, condition, x, y, x_divergence, y_divergence, blok))
    return(result_list)
  }
  
  loadRoisAcq = function(files.roi.path) {
    rois_acq <- data.frame(
      subject = integer(0),
      trial = integer(0),
      stim = character(0),
      roi.xleft = integer(0),
      roi.xright = integer(0),
      roi.ybottom = integer(0),
      roi.ytop = integer(0),
      quadrant.xleft = integer(0),
      quadrant.xright = integer(0),
      quadrant.ybottom = integer(0),
      quadrant.ytop = integer(0)
    )
    for (file.roi.path in files.roi.path) {
      # file.roi.path = files.roi.path[1]
      roi_acq = read_csv(file.roi.path)
      roi_acq <- roi_acq %>%
        select(position, stim)
      roi_acq <- na.omit(roi_acq)
      # In normalised (‘norm’) units the window ranges in both x and y from -1 to +1. That is, the top right of the window has coordinates (1,1), the bottom left is (-1,-1).
      roi_acq <- roi_acq %>%
        mutate(trial = row_number(),
               x_position_norm = position %>% sub("\\(", "", .) %>% sub(", .*", "", .) %>% as.numeric(),
               y_position_norm = position %>% sub(".*, ", "", .) %>% sub("\\)", "", .) %>% as.numeric(),
               x_position = screen.width / 2 + (x_position_norm * screen.width / 2),
               y_position = screen.height / 2 + (y_position_norm * screen.height / 2),
               roi.xleft = x_position - image.acq.width / 2,
               roi.xright = x_position + image.acq.width / 2,
               roi.ytop = y_position - image.acq.height / 2,
               roi.ybottom = y_position + image.acq.height / 2,
               quadrant.xleft = ifelse(x_position_norm < 0, 0, screen.width / 2),
               quadrant.xright = ifelse(x_position_norm < 0, screen.width / 2, screen.width),
               quadrant.ytop = ifelse(y_position_norm < 0, 0, screen.height / 2),
               quadrant.ybottom = ifelse(y_position_norm < 0, screen.height / 2, screen.height))
      roi_acq <- roi_acq %>% 
        mutate(subject = file.roi.path %>% sub(".*attentional_competition_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer())
      roi_acq <- roi_acq %>%
        select(subject, trial, stim, roi.xleft, roi.xright, roi.ybottom, roi.ytop, quadrant.xleft, quadrant.xright, quadrant.ybottom, quadrant.ytop)
      
      rois_acq <- rbind(rois_acq, roi_acq)
    }
    
    return(rois_acq)
  }
  
  loadRoisTest = function(files.roi.path) {
    rois_test <- data.frame(
      subject = integer(0),
      trial = integer(0),
      stim1 = character(0),
      roi1.xleft = integer(0),
      roi1.xright = integer(0),
      roi1.ybottom = integer(0),
      roi1.ytop = integer(0),
      stim2 = character(0),
      roi2.xleft = integer(0),
      roi2.xright = integer(0),
      roi2.ybottom = integer(0),
      roi2.ytop = integer(0),
      stim3 = character(0),
      roi3.xleft = integer(0),
      roi3.xright = integer(0),
      roi3.ybottom = integer(0),
      roi3.ytop = integer(0),
      stim4 = character(0),
      roi4.xleft = integer(0),
      roi4.xright = integer(0),
      roi4.ybottom = integer(0),
      roi4.ytop = integer(0)
    )
    for (file.roi.path in files.roi.path) {
      # file.roi.path = files.roi.path[4]
      roi_test = read_csv(file.roi.path)
      roi_test <- roi_test %>%
        select(position1, stim1, position2, stim2, position3, stim3, position4, stim4, testtrials.thisN, testtrials_novelty.thisN) %>% 
        mutate(trial = ifelse(is.na(testtrials.thisN), testtrials_novelty.thisN + test1End + 1, testtrials.thisN + acq2End + 1)) %>% 
        select(-c(testtrials.thisN, testtrials_novelty.thisN))
      
      roi_test <- na.omit(roi_test)
      
      roi_test <- roi_test %>% 
        mutate(new_trial = ifelse((trial != lag(trial)), TRUE, FALSE)) %>% 
        mutate(new_trial = ifelse(is.na(new_trial), TRUE, new_trial))
      
      # if (sum(roi_test$new_trial) < 96) {
      #   roi_test$new_trial = TRUE
      # }
      roi_test <- roi_test %>% 
        filter(new_trial) %>% 
        select(-new_trial)
      
      # In normalised (‘norm’) units the window ranges in both x and y from -1 to +1. That is, the top right of the window has coordinates (1,1), the bottom left is (-1,-1).
      roi_test <- roi_test %>%
        mutate(x_position1_norm = position1 %>% sub("\\(", "", .) %>% sub(", .*", "", .) %>% as.numeric(),
               y_position1_norm = position1 %>% sub(".*, ", "", .) %>% sub("\\)", "", .) %>% as.numeric(),
               x_position1 = screen.width / 2 + (x_position1_norm * screen.width / 2),
               y_position1 = screen.height / 2 + (y_position1_norm * screen.height / 2),
               roi1.xleft = x_position1 - image.acq.width / 2,
               roi1.xright = x_position1 + image.acq.width / 2,
               roi1.ytop = y_position1 - image.acq.height / 2,
               roi1.ybottom = y_position1 + image.acq.height / 2,
               
               x_position2_norm = position2 %>% sub("\\(", "", .) %>% sub(", .*", "", .) %>% as.numeric(),
               y_position2_norm = position2 %>% sub(".*, ", "", .) %>% sub("\\)", "", .) %>% as.numeric(),
               x_position2 = screen.width / 2 + (x_position2_norm * screen.width / 2),
               y_position2 = screen.height / 2 + (y_position2_norm * screen.height / 2),
               roi2.xleft = x_position2 - image.acq.width / 2,
               roi2.xright = x_position2 + image.acq.width / 2,
               roi2.ytop = y_position2 - image.acq.height / 2,
               roi2.ybottom = y_position2 + image.acq.height / 2,
               
               x_position3_norm = position3 %>% sub("\\(", "", .) %>% sub(", .*", "", .) %>% as.numeric(),
               y_position3_norm = position3 %>% sub(".*, ", "", .) %>% sub("\\)", "", .) %>% as.numeric(),
               x_position3 = screen.width / 2 + (x_position3_norm * screen.width / 2),
               y_position3 = screen.height / 2 + (y_position3_norm * screen.height / 2),
               roi3.xleft = x_position3 - image.acq.width / 2,
               roi3.xright = x_position3 + image.acq.width / 2,
               roi3.ytop = y_position3 - image.acq.height / 2,
               roi3.ybottom = y_position3 + image.acq.height / 2,
               
               x_position4_norm = position4 %>% sub("\\(", "", .) %>% sub(", .*", "", .) %>% as.numeric(),
               y_position4_norm = position4 %>% sub(".*, ", "", .) %>% sub("\\)", "", .) %>% as.numeric(),
               x_position4 = screen.width / 2 + (x_position4_norm * screen.width / 2),
               y_position4 = screen.height / 2 + (y_position4_norm * screen.height / 2),
               roi4.xleft = x_position4 - image.acq.width / 2,
               roi4.xright = x_position4 + image.acq.width / 2,
               roi4.ytop = y_position4 - image.acq.height / 2,
               roi4.ybottom = y_position4 + image.acq.height / 2)
      
      roi_test <- roi_test %>% 
        mutate(subject = file.roi.path %>% sub(".*attentional_competition_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer())
      roi_test <- roi_test %>%
        select(subject, trial,
               stim1, roi1.xleft, roi1.xright, roi1.ybottom, roi1.ytop,
               stim2, roi2.xleft, roi2.xright, roi2.ybottom, roi2.ytop,
               stim3, roi3.xleft, roi3.xright, roi3.ybottom, roi3.ytop,
               stim4, roi4.xleft, roi4.xright, roi4.ybottom, roi4.ytop)
      rois_test <- rbind(rois_test, roi_test)
    }
    
    return(rois_test)
  }
  
  checkRoi <- function(x, y, rect_xleft, rect_xright, rect_ybottom, rect_ytop) {
    return(x >= rect_xleft & x <= rect_xright & y <= rect_ybottom & y >= rect_ytop)
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
# Eye-Tracking
###############################################################################
conditions = loadConditions(files.cond) %>%  # conditions
  mutate(across('condition', str_replace_all, rep_str))
write.csv2(conditions, file.path(path, "Physio", "Trigger", "conditions.csv"), row.names=FALSE, quote=FALSE)

roisAcq = loadRoisAcq(files.cond) # rois
roisTest = loadRoisTest(files.cond) # rois

path.eye = file.path(path, "Experiment", "attentional_competition_task", "data") # eye tracking data
messages = loadMessages(file.path(path.eye, "messages.txt")) %>% # messages for timepoints
  mutate(trial = trial - 2) %>% 
  filter(trial > 0) %>% 
  mutate(across('event', str_replace_all, rep_str)) %>% 
  filter(time != 0)
write.csv2(messages, file.path(path, "Gaze", "messages.csv"), row.names=FALSE, quote=FALSE)

fixations = loadFixations(file.path(path.eye, "fixations.txt"), screen.height) %>% # fixations
  mutate(trial = trial - 2) %>% 
  filter(trial > 0) %>% 
  left_join(conditions, by=c("subject", "trial")) %>% # get conditions & other variables
  left_join(messages %>% filter(grepl("ImageOnset", event)) %>% # get image onset
              rename(picOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  left_join(messages %>% filter(grepl("FeedbackOnset", event)) %>% # get feedback onset
              rename(feedbackOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  left_join(messages %>% filter(grepl("ImageTestOnset", event)) %>% # get test onset
              rename(picTestOnset = time) %>% distinct(subject, trial, .keep_all=T) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  mutate(picOnset = ifelse(is.na(picOnset), picTestOnset, picOnset)) %>% 
  select(-picTestOnset) %>% 
  mutate(feedbackOnset = ifelse(is.na(feedbackOnset), picOnset + 2000, feedbackOnset))
vpn.eye = fixations$subject %>% unique() %>% setdiff(exclusions.eye.num) %>% sort()

###############################################################################
# Acquisition
###############################################################################
# determine valid acquisition trials by fixation time (invalid trials need not be evaluated by their baseline)
eye.valid.trial = fixations %>% filter(grepl("acquisition", phase)) %>% # filter for acq trials
  mutate(feedbackOnset = feedbackOnset - picOnset, # realign such that 0 = picture onset
         start = start - picOnset, end = end - picOnset, # realign such that 0 = picture onset
         start = ifelse(start < 0, 0, start), # discard fraction of fixation before stimulus
         end = ifelse(end > feedbackOnset + 200, feedbackOnset + 200, end), # discard fraction of fixation after feedback onset
         dur = end - start) %>% filter(dur > 0) %>% 
  group_by(subject, trial) %>% summarise(valid = sum(dur) / mean(feedbackOnset + 200))
with(eye.valid.trial, hist(valid, breaks=seq(0, 1, length.out=20+1), main=paste0("Valid Fixation Time"))); abline(v=validFixTime.trial, col="red", lwd=3, lty=2)

#exclude trials with insufficient valid fixations (need not be validated for their baseline)
fixations.acq.valid = eye.valid.trial %>% filter(valid > validFixTime.trial) %>% 
  left_join(fixations %>% mutate(trial = as.numeric(trial)), by=c("subject", "trial")) %>% select(-valid)

# Baseline Validation
baseline.acq.validation = validateBaselines(fixations.acq.valid, messages, exclusions, maxDeviation_rel, maxSpread, saveBaselinePlots, postfix="acq")
baseline.acq.summary <- baseline.acq.validation$baseline.summary
baseline.acq.trials <- baseline.acq.validation$baseline.trials
baseline.acq.trials <- baseline.acq.trials %>% 
  left_join(baseline.acq.summary %>% select(subject, invalid) %>% mutate(subject = subject %>% as.numeric())) %>% 
  mutate(blok = ifelse(invalid > outlierLimit.eye, FALSE, blok))

baseline.acq.summary = baseline.acq.summary %>% 
  mutate(included = invalid <= outlierLimit.eye & range_x <= maxSpread & range_y <= maxSpread)
with(baseline.acq.summary, hist(invalid, breaks=max(ntrials), main=paste0("Valid Baselines"))); abline(v = outlierLimit.eye, col="red", lwd=2, lty=2)

baseline.acq.summary %>% summarise(totalN = n(), includedN = sum(included), includedP = mean(included))
baseline.acq.summary %>% group_by(subject) %>% summarise(invalid = mean(invalid)) %>% arrange(desc(invalid))
baseline.acq.summary %>% summarise(mean_bl = mean(invalid)*100, sd_bl = sd(invalid)*100)
excluded_subjects = baseline.acq.summary %>% filter(included == FALSE)
excluded_subjects = as.numeric(excluded_subjects$subject)
print(paste0("Number of excluded subjects (Acquisition): ", length(excluded_subjects)))


### SACCADES
saccades = loadSaccades(file.path(path.eye, "saccades.txt"), screen.height) %>% 
  mutate(trial = trial - 2) %>% 
  filter(trial > 0) %>%
  left_join(conditions, by=c("subject", "trial")) %>% # get conditions & other variables
  left_join(messages %>% filter(grepl("ImageOnset", event)) %>% # get image onset
              rename(picOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  left_join(messages %>% filter(grepl("FeedbackOnset", event)) %>% # get feedback onset
              rename(feedbackOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  left_join(messages %>% filter(grepl("ImageTestOnset", event)) %>% # get test onset
              rename(picTestOnset = time) %>% distinct(subject, trial, .keep_all=T) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  mutate(picOnset = ifelse(is.na(picOnset), picTestOnset, picOnset)) %>% 
  select(-picTestOnset) %>% 
  mutate(feedbackOnset = ifelse(is.na(feedbackOnset), picOnset + 2000, feedbackOnset))

# ggplot(saccades, aes(x = end_x, y = end_y, color = subject)) +
#   geom_point(size=0.1) +
#   xlim(0, screen.width) +
#   ylim(0, screen.height)
# ggsave(file.path(path, "Plots", "Gaze", "fixations.png"), width=screen.width, height=screen.height, units="px")

# Baseline OK?
saccades.acq.valid <- saccades %>% 
  left_join(baseline.acq.trials %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial")) %>% 
  filter(trial <= acq2End) %>% 
  filter(!subject %in% excluded_subjects)

# Correction for baseline deviations 
saccades.acq.valid <- saccades.acq.valid %>% 
  mutate(start_x_corr = start_x + x_divergence, end_x_corr = end_x + x_divergence,
         start_y_corr = start_y + y_divergence, end_y_corr = end_y + y_divergence)

# Correct for picture onset
saccades.acq.valid = saccades.acq.valid %>%
  mutate(feedbackOnset = feedbackOnset - picOnset, # realign such that 0 = picture onset
         start_time = start_time - picOnset, end_time = end_time - picOnset, # realign such that 0 = picture onset
         start_time = ifelse(start_time < 0, 0, start_time), # discard fraction of fixation before stimulus
         picOnset = picOnset - picOnset,
         dur = end_time - start_time)

saccades.acq.valid = saccades.acq.valid %>% filter(dur > 0) # drop saccades before picture onset

saccades.acq.valid = saccades.acq.valid %>% filter(start_time < feedbackOnset)

# Outcome OK?
saccades.acq.valid <- saccades.acq.valid %>% 
  left_join(roisAcq, by=c("subject", "trial"))

saccades.acq.valid <- saccades.acq.valid %>% 
  mutate(ROI = checkRoi(end_x_corr, end_y_corr, roi.xleft, roi.xright, roi.ybottom, roi.ytop),
         Quadrant = checkRoi(end_x_corr, end_y_corr, quadrant.xleft, quadrant.xright, quadrant.ybottom, quadrant.ytop))

saccades.acq.valid <- saccades.acq.valid %>% 
  mutate(outcome_corr = ifelse(ROI & str_detect(condition, "neg"), "loss", ifelse(ROI & str_detect(condition, "pos"), "reward", NA)))

saccades.acq.valid <- saccades.acq.valid %>% 
  group_by(subject) %>%
  arrange(trial, match(ROI, c(TRUE, FALSE))) %>% 
  distinct(trial, .keep_all = TRUE) %>% 
  arrange(subject, trial) %>% 
  ungroup

saccades.acq.valid <- saccades.acq.valid %>% 
  mutate(outcome = coalesce(outcome, "none")) %>% 
  mutate(outcome_corr = coalesce(outcome_corr, "none")) %>%
  mutate(outcomeok = ifelse(outcome == outcome_corr, TRUE, FALSE))

# Calculate lengths of saccades
saccades.acq.analysis = saccades.acq.valid %>% 
  mutate(start_x_corr_cm = pixToCm(start_x_corr, screen.width, screen.width.cm),
         end_x_corr_cm = pixToCm(end_x_corr, screen.width, screen.width.cm),
         start_y_corr_cm = pixToCm(start_y_corr, screen.height, screen.height.cm),
         end_y_corr_cm = pixToCm(end_y_corr, screen.height, screen.height.cm)) %>% 
  mutate(angle = visangle(start_x_corr_cm, start_y_corr_cm, end_x_corr_cm, end_y_corr_cm, distance))

saccades.acq.analysis <- saccades.acq.analysis %>% 
  mutate(block = ifelse(trial <= acq1End, 0, 1)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

saccades.acq.analysis <- saccades.acq.analysis %>% 
  mutate(correct_reaction = ifelse(str_detect(condition_threat, "neg") & str_detect(outcome_corr, "none"), TRUE, ifelse(str_detect(condition_threat, "pos") & str_detect(outcome_corr, "reward"), TRUE, FALSE)))

# Proportion of wrong outcome and invalid baselines
saccades.acq.analysis %>% 
  summarise(mean_outcome_trial = mean(outcomeok), .by=c(subject, trial)) %>% 
  summarise(mean_outcome_subj = (1-mean(mean_outcome_trial)) * 100, .by=c(subject)) %>% 
  summarise(mean_outcome = mean(mean_outcome_subj), sd_outcome = sd(mean_outcome_subj))

# Add scores and write saccades to CSV
saccades.acq.analysis <- saccades.acq.analysis %>%
  left_join(scores, by="subject")

write.csv2(saccades.acq.analysis, file.path(path, "Gaze", "saccades_acq.csv"), row.names=FALSE, quote=FALSE)

# Percentage of saccades going towards the stimuli
saccades.acq.prop <- saccades.acq.analysis %>%
  # filter(motivation_points > 4) %>%  # Filter for Motivation
  # filter (trial >= 32) %>%
  filter(blok) %>%
  filter(outcomeok) %>%
  filter(!contains_blink) %>%
  summarise(absolute_frequency_ROI = sum(angle >= 1 & ROI & start_time < feedbackOnset), relative_frequency_ROI = mean(angle >= 1 & ROI & start_time < feedbackOnset), relative_frequency_Quadrant = mean(Quadrant & start_time < feedbackOnset), .by=c(subject, SPAI, condition, condition_social, condition_threat))

# saccades.acq.prop.wide <- saccades.acq.prop %>%
#   select(c(subject, condition, relative_frequency_ROI)) %>%
#   pivot_wider(names_from = condition, values_from = relative_frequency_ROI)
# 
# saccades.acq.prop.wide <- saccades.acq.prop.wide %>%
#   left_join(scores, by="subject") %>% 
#   rename_with(~ gsub(",\n", "_", .x, fixed = TRUE))
# 
# write.csv2(saccades.acq.prop.wide, file.path(path, "sacc_prop_acq_wide.csv"), row.names=FALSE, quote=FALSE)

saccades.acq.roi.summary <- saccades.acq.prop %>%
  summarise(Mean = mean(relative_frequency_ROI), SD = sd(relative_frequency_ROI), .by=condition)

plot_saccades_acq <- ggplot(saccades.acq.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.acq.prop, aes(y = relative_frequency_ROI), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
  # labs(title = paste("Proportion of Trials with a Saccade to the Stimulus (N = ", n_distinct(saccades.acq.prop$subject), ")", sep=""), x = "Conditions", y = "Proportion") +
  labs(title = paste("Proportion of Approach Trials", sep=""), x = NULL, y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

ggsave(file.path(path, "Plots", "Gaze", "acquisition", "saccades_proportion_roi.png"), width=1800, height=2000, units="px")

anova <- saccades.acq.prop %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
  ez::ezANOVA(dv=.(relative_frequency_ROI),
              wid=.(subject),
              within=.(condition_threat, condition_social),
              # between=.(SPAI),
              detailed=T, type=3)
anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # threat x social

saccades.acq.prop %>%
  filter(condition_threat == "neg") %>%
  apa::t_test(data=., relative_frequency_ROI ~ condition_social, paired = T, alternative = "two.sided") %>% 
  apa::t_apa()

# saccades.acq.prop %>%
#   filter(condition_threat == "neg") %>%
#   filter(condition_social == "social") %>%
#   cor.test(x=.$relative_frequency_ROI, y=.$SPAI, method="pearson", alternative = "two.sided") %>% 
#   apa::cor_apa()

# # Saccades to ROI incl. Microsaccades
# saccades.acq.prop <- saccades.acq.analysis %>%
#   filter(blok) %>% 
#   filter(outcomeok) %>%
#   filter(!contains_blink) %>%
#   summarise(relative_frequency_ROI = mean(ROI), .by=c(subject, condition)) 
# 
# saccades.acq.prop <- saccades.acq.analysis %>% 
#   # filter (trial >= 32) %>%
#   filter(blok) %>% 
#   filter(outcomeok) %>%
#   filter(!contains_blink) %>%
#   summarise(absolute_frequency_ROI = sum(ROI & start_time < feedbackOnset), relative_frequency_ROI = mean(ROI & start_time < feedbackOnset), relative_frequency_Quadrant = mean(Quadrant & start_time < feedbackOnset), .by=c(subject, SPAI, condition, condition_social, condition_threat)) 
# 
# saccades.acq.roi.summary <- saccades.acq.prop %>% 
#   summarise(Mean = mean(relative_frequency_ROI), SD = sd(relative_frequency_ROI), .by=condition)
# 
# ggplot(saccades.acq.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
#   geom_col(position = "dodge", width = 0.7) +
#   geom_errorbar(
#     aes(ymin = Mean - SD, ymax = Mean + SD),
#     position = position_dodge(width = 0.7),
#     width = 0.25) +
#   geom_point(data=saccades.acq.prop, aes(y = relative_frequency_ROI), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
#   labs(title = paste("Proportion of Trials with a Saccade to the Stimulus (N = ", n_distinct(saccades.acq.prop$subject), ")", sep=""), x = "Conditions", y = "Proportion") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_fill_viridis_d() + 
#   scale_color_viridis_d()
# 
# ggsave(file.path(path, "Plots", "Gaze", "acquisition", "microsaccades_proportion_roi.png"), width=1800, height=2000, units="px")
# 
# saccades.acq.prop %>% 
#   group_by(subject, condition_social, condition_threat) %>% 
#   mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
#   ez::ezANOVA(dv=.(relative_frequency_ROI),
#               wid=.(subject), 
#               within=.(condition_social, condition_threat), 
#               # between=.(SPAI),
#               detailed=T, type=3) %>% 
#   apa::anova_apa()
# 
# 
# # Line Plot
# saccades.acq.prop <- saccades.acq.analysis %>%
#   filter(blok) %>%
#   filter(outcomeok) %>%
#   filter(!contains_blink) %>%
#   summarise(absolute_frequency_ROI = sum(ROI & start_time < feedbackOnset), relative_frequency_ROI = mean(ROI & start_time < feedbackOnset), relative_frequency_Quadrant = mean(Quadrant & start_time < feedbackOnset), .by=c(subject, block, condition))
# 
# saccades.acq.roi.summary <- saccades.acq.prop %>%
#   summarise(Mean = mean(relative_frequency_ROI), SD = sd(relative_frequency_ROI), .by=c(condition, block))
# 
# ggplot(saccades.acq.roi.summary, aes(x = block, y = Mean, group = condition, color = condition)) +
#   geom_point(position = position_dodge(width = 0.5), shape = "square", size = 5) +
#   geom_errorbar(
#     aes(ymin = Mean - SD, ymax = Mean + SD),
#     position = position_dodge(width = 0.5),
#     width = 0.25, linewidth = 1) +
#   geom_line(position = position_dodge(width = 0.5), linewidth = 1) +
#   labs(title = paste("Proportion of Trials with a Saccade to the Stimulus (N = ", n_distinct(avoidance.acq.prop$subject), ")", sep=""), x = "Blocks", y = "Proportion") +
#   guides(color = guide_legend(title = "Conditions")) +
#   theme_minimal() +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   scale_x_continuous(breaks=c(0,1,2), labels=c("Block 1\n1-32", "Block 2\n33-72", "Block 3\n73-112"))
# 
# ggsave(file.path(path, "Plots", "Gaze", "acquisition", "saccades_proportion_roi_block.png"), width=2800, height=2000, units="px")

# Latency to first saccade going towards the stimuli
saccades.acq.lat.roi <- saccades.acq.analysis %>%
  # filter(motivation_points > 4) %>%  # Filter for Motivation
  # filter (trial >= 32) %>%
  filter(blok) %>%
  filter(outcomeok) %>%
  filter(!contains_blink) %>%
  filter(angle >= 1) %>%
  filter(ROI & start_time < feedbackOnset) %>%
  reframe(latency = mean(start_time), .by=c(subject, SPAI, condition, condition_social, condition_threat))

saccades.acq.lat.roi.filter <- saccades.acq.lat.roi%>%
  reframe(n_cond = length(unique(condition)), .by=c(subject))

saccades.acq.lat.roi.summary <- saccades.acq.lat.roi %>%
  summarise(Mean = mean(latency), SD = sd(latency), .by=c(condition, condition_social, condition_threat))

plot_latencies_acq <- ggplot(saccades.acq.lat.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.acq.lat.roi, aes(y = latency), size = 2, shape = 21, color = "black", alpha=0.3, position = position_jitter(width=0.2, height=0.005)) +
  # labs(title = paste("Latency to First Saccade to the Stimulus (N = ", n_distinct(saccades.acq.lat.roi$subject), ")", sep=""), x = "Conditions", y = "Latency [ms]") +
  labs(title = paste("Latency of First Saccade", sep=""), x = NULL, y = "Latency [ms]") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()


ggsave(file.path(path, "Plots", "Gaze", "acquisition", "saccades_latency_roi.png"), width=1800, height=2000, units="px")

anova <- saccades.acq.lat.roi %>%
  left_join(saccades.acq.lat.roi.filter, by=c("subject")) %>%
  filter(n_cond == 4) %>%
  group_by(subject, condition_social, condition_threat) %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
  ez::ezANOVA(dv=.(latency),
              wid=.(subject),
              within=.(condition_threat, condition_social),
              # between=.(SPAI),
              detailed=T, type=3)
anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # threat x social

# saccades.acq.lat.roi %>%
#   filter(condition_social == "social") %>%
#   t_test(data=., latency ~ condition_threat, paired = T, alternative = "two.sided") %>% 
#   apa::t_apa()
# 
# saccades.acq.lat.roi %>%
#   filter(condition_threat == "neg") %>%
#   t_test(data=., latency ~ condition_social, paired = T, alternative = "two.sided") %>% 
#   apa::t_apa()

# # Length of saccade going towards the stimuli
# saccades.acq.len.roi <- saccades.acq.analysis %>%
#   # filter (trial >= 32) %>%
#   filter(blok) %>%
#   filter(outcomeok) %>%
#   filter(!contains_blink) %>%
#   # filter(angle >= 1) %>%
#   filter(ROI & start_time < feedbackOnset) %>%
#   summarise(length = mean(angle), .by=c(subject, SPAI, condition, condition_social, condition_threat))
# 
# saccades.acq.len.roi.filter <- saccades.acq.len.roi%>%
#   reframe(n_cond = length(unique(condition)), .by=c(subject))
# 
# saccades.acq.len.roi.summary <- saccades.acq.len.roi %>%
#   summarise(Mean = mean(length), SD = sd(length), .by=c(condition, condition_social, condition_threat))
# 
# ggplot(saccades.acq.len.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
#   geom_col(position = "dodge", width = 0.7) +
#   geom_errorbar(
#     aes(ymin = Mean - SD, ymax = Mean + SD),
#     position = position_dodge(width = 0.7),
#     width = 0.25) +
#   geom_point(data=saccades.acq.len.roi, aes(y = length), size = 2, shape = 21, color = "black", alpha=0.3, position = position_jitter(width=0.2, height=0.005)) +
#   labs(title = paste("Length of Saccade to the Stimulus (N = ", n_distinct(saccades.acq.len.roi$subject), ")", sep=""), x = "Conditions", y = "Length [degree visual angle]") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d()
# 
# ggsave(file.path(path, "Plots", "Gaze", "acquisition", "saccades_length_roi.png"), width=1800, height=2000, units="px")
# 
# 
# saccades.acq.len.roi %>%
#   left_join(saccades.acq.len.roi.filter, by=c("subject")) %>%
#   filter(n_cond == 4) %>%
#   group_by(subject, condition_social, condition_threat) %>%
#   mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
#   ez::ezANOVA(dv=.(length),
#               wid=.(subject),
#               within=.(condition_social, condition_threat),
#               # between=.(SPAI),
#               detailed=T, type=3) %>%
#   apa::anova_apa()


# Percentage of avoidance trials
# Save dataframe
avoidance.acq.prop <- saccades.acq.analysis %>%
  filter(blok) %>%
  filter(outcomeok) %>%
  summarise(relative_frequency_av = mean(outcome_corr=="none"), .by=c(subject, condition))

# avoidance.acq.prop.wide <- avoidance.acq.prop %>%
#   pivot_wider(names_from = condition, values_from = relative_frequency_av)
# 
# avoidance.acq.prop.wide <- avoidance.acq.prop.wide %>%
#   left_join(scores, by="subject")
# 
# write.csv2(avoidance.acq.prop.wide, file.path(path, "avoidance_wide.csv"), row.names=FALSE, quote=FALSE)

# Bar Plot
avoidance.acq.prop <- saccades.acq.analysis %>%
  # filter (trial >= 32) %>%
  filter(blok) %>%
  filter(outcomeok) %>%
  summarise(absolute_frequency_av = sum(outcome_corr=="none"), relative_frequency_av = mean(outcome_corr=="none"), .by=c(subject, SPAI, condition, condition_social, condition_threat))

avoidance.acq.prop.summary <- avoidance.acq.prop %>%
  summarise(Mean = mean(relative_frequency_av), SD = sd(relative_frequency_av), .by=condition)

# ggplot(avoidance.acq.prop.summary, aes(x = condition, y = Mean, fill = condition)) +
#   geom_col(position = "dodge", width = 0.7) +
#   geom_errorbar(
#     aes(ymin = Mean - SD, ymax = Mean + SD), #ymin = Mean - SD,
#     position = position_dodge(width = 0.7),
#     width = 0.25) +
#   geom_line(data=avoidance.acq.prop, aes(y = relative_frequency_av, group = subject), alpha=0.1) +
#   geom_point(data=avoidance.acq.prop, aes(y = relative_frequency_av), size = 2, shape = 21, color = "black", alpha=0.3) + # , position=position_jitter(width=0.05)) +
#   labs(title = paste("Proportion of Avoidance-Trials (N = ", n_distinct(avoidance.acq.prop$subject), ")", sep=""), x = "Conditions", y = "Proportion") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d()
# 
# ggsave(file.path(path, "Plots", "Gaze", "acquisition", "avoidance_proportion_roi.png"), width=1800, height=2000, units="px")
# 
# avoidance.acq.prop %>%
#   group_by(subject, condition_social, condition_threat) %>%
#   mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
#   ez::ezANOVA(dv=.(relative_frequency_av), wid=.(subject),
#               within=.(condition_social, condition_threat),
#               # between=.(SPAI),
#               detailed=T, type=3) %>%
#   apa::anova_apa()
# 
# # Line Plot
# avoidance.acq.prop <- saccades.acq.analysis %>%
#   filter(blok) %>%
#   filter(outcomeok) %>%
#   summarise(absolute_frequency_av = sum(outcome_corr=="none"), relative_frequency_av = mean(outcome_corr=="none"), .by=c(subject, block, condition))
# 
# avoidance.acq.prop.summary <- avoidance.acq.prop %>%
#   summarise(Mean = mean(relative_frequency_av), SD = sd(relative_frequency_av), .by=c(condition, block))
# 
# ggplot(avoidance.acq.prop.summary, aes(x = block, y = Mean, group = condition, color = condition)) +
#   geom_point(position = position_dodge(width = 0.5), shape = "square", size = 5) +
#   geom_errorbar(
#     aes(ymin = Mean - SD, ymax = Mean + SD),
#     position = position_dodge(width = 0.5),
#     width = 0.25, linewidth = 1) +
#   geom_line(position = position_dodge(width = 0.5), linewidth = 1) +
#   labs(title = paste("Proportion of Avoidance-Trials (N = ", n_distinct(avoidance.acq.prop$subject), ")", sep=""), x = "Blocks", y = "Proportion") +
#   guides(color = guide_legend(title = "Conditions")) +
#   theme_minimal() +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d() +
#   scale_x_continuous(breaks=c(0,1), labels=c("Block 1\n1-32", "Block 2\n33-64"))
# 
# ggsave(file.path(path, "Plots", "Gaze", "acquisition", "avoidance_proportion_roi_block.png"), width=2800, height=2000, units="px")


# Correlation Correct Avoidance and Discrimination
discrimination <- read.csv2(file.path(path, "Behavior", "discrimination.csv"), sep=";", dec=",")

avoidance_corr.acq.prop <- saccades.acq.analysis %>%
  filter(blok) %>%
  filter(outcomeok) %>%
  summarise(absolute_frequency_corr = sum(correct_reaction), relative_frequency_corr = mean(correct_reaction), .by=c(subject, condition_social))

avoidance_corr.acq.prop <- avoidance_corr.acq.prop %>%
  left_join(discrimination, by=c("subject", "condition_social")) %>%
  rename(Category = condition_social)

cor.test(avoidance_corr.acq.prop$relative_frequency_corr, avoidance_corr.acq.prop$discrimination) %>% apa::cor_apa()

ggplot(avoidance_corr.acq.prop, aes(x = discrimination, y = relative_frequency_corr, color = Category, fill = Category)) +
  geom_point(shape = 21, size = 2, color = "black", alpha=0.1) +
  geom_smooth(method = "lm", aes(fill=Category), alpha = 0.1) +
  # labs(title = paste("Correlation of Discrimination and Avoidance Performance (N = ", n_distinct(avoidance_corr.acq.prop$subject), ")", sep=""), x = "Discrimination Score", y = "Performance Score") +
  labs(title = "Correlation of Discrimination and Performance", x = "Discrimination Score", y = "Performance Score") +
  theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme(legend.position="top", legend.box = "horizontal", legend.title=element_blank())
ggsave(file.path(path, "Plots", "Gaze", "acquisition", "corr_discrimination_avoidance.png"), width=1600, height=1600, units="px")


### FIXATIONS
# Baseline OK?
fixations.acq.valid <- fixations.acq.valid %>% 
  left_join(baseline.acq.trials %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial")) %>% 
  filter(trial <= acq2End) %>% 
  filter(!subject %in% excluded_subjects)

# Correction for baseline deviations 
fixations.acq.valid <- fixations.acq.valid %>% 
  mutate(x_corr = x + x_divergence,
         y_corr = y + y_divergence)

# Correct for picture onset
fixations.acq.valid = fixations.acq.valid %>%
  mutate(feedbackOnset = feedbackOnset - picOnset, # realign such that 0 = picture onset
         start = start - picOnset, end = end - picOnset, # realign such that 0 = picture onset
         start = ifelse(start < 0, 0, start), # discard fraction of fixation before stimulus
         picOnset = picOnset - picOnset,
         dur = end - start)

fixations.acq.analysis = fixations.acq.valid %>% filter(dur > 0) # drop saccades before picture onset

# ROIs
fixations.acq.analysis <- fixations.acq.analysis %>% 
  left_join(roisAcq, by=c("subject", "trial"))

fixations.acq.analysis <- fixations.acq.analysis %>% 
  mutate(ROI = checkRoi(x_corr,y_corr, roi.xleft, roi.xright, roi.ybottom, roi.ytop),
         Quadrant = checkRoi(x_corr, y_corr, quadrant.xleft, quadrant.xright, quadrant.ybottom, quadrant.ytop))

# Add scores and write saccades to CSV
fixations.acq.analysis <- fixations.acq.analysis %>%
  left_join(scores, by="subject") %>% 
  ungroup

write.csv2(fixations.acq.analysis, file.path(path, "Gaze", "fixations_acq.csv"), row.names=FALSE, quote=FALSE)

fixations.acq.analysis <- fixations.acq.analysis %>% 
  mutate(block = ifelse(trial <= acq1End, 0, 1)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

# Mean Dwell Time on stimuli
fixations.acq.dwell <- fixations.acq.analysis %>%
  # filter(motivation_points > 4) %>%  # Filter for Motivation
  # filter (trial >= 32) %>%
  filter(blok) %>%
  filter(ROI) %>% 
  summarise(mean_dwell_time = mean(dur), .by=c(subject, SPAI, condition, condition_social, condition_threat)) %>% 
  complete(subject, condition) %>% 
  mutate(mean_dwell_time = replace_na(mean_dwell_time, 0)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

fixations.acq.roi.summary <- fixations.acq.dwell %>%
  summarise(Mean = mean(mean_dwell_time), SD = sd(mean_dwell_time), .by=condition)

plot_fixations_acq <- ggplot(fixations.acq.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=fixations.acq.dwell, aes(y = mean_dwell_time), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
  # labs(title = paste("Dwell Times on the Stimulus (N = ", n_distinct(fixations.acq.dwell$subject), ")", sep=""), x = "Conditions", y = "Dwell time [ms]") +
  labs(title = paste("Dwell Time", sep=""), x = NULL, y = "Dwell time [ms]") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

ggsave(file.path(path, "Plots", "Gaze", "acquisition", "dwell_time_roi.png"), width=1800, height=2000, units="px")

anova <- fixations.acq.dwell %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
  ez::ezANOVA(dv=.(mean_dwell_time),
              wid=.(subject),
              within=.(condition_threat, condition_social),
              # between=.(SPAI),
              detailed=T, type=3)
anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # threat x social


attention_acq_plot <- plot_grid(plot_saccades_acq, plot_latencies_acq, plot_fixations_acq, ncol = 3, labels=c("A", "B", "C"))
ggsave(file.path(path, "Plots", "Gaze", "acquisition", paste0("gaze_acq.png")), width=3500, height=1200, units="px")


###############################################################################
# Competition
###############################################################################
# Determine valid test trials by fixation time (invalid trials need not be evaluated by their baseline)
end.testtrial = 2000
eye.valid.trial = fixations  %>% filter(grepl("test", phase)) %>% # filter for test trials
  mutate(end.testtrial = end.testtrial, # realign such that 0 = picture onset
         start = start - picOnset, end = end - picOnset, # realign such that 0 = picture onset
         start = ifelse(start < 0, 0, start), # discard fraction of fixation before stimulus
         end = ifelse(end > end.testtrial, end.testtrial, end), # discard fraction of fixation after end
         dur = end - start) %>% filter(dur > 0) %>% 
  group_by(subject, trial) %>% summarise(valid = sum(dur) / mean(end.testtrial))
with(eye.valid.trial, hist(valid, breaks=seq(0, 1, length.out=20+1), main=paste0("Valid Fixation Time"))); abline(v=validFixTime.trial, col="red", lwd=3, lty=2)

#exclude trials with insufficient valid fixations (need not be validated for their baseline)
fixations.test.valid = eye.valid.trial %>% filter(valid > validFixTime.trial) %>% 
  left_join(fixations %>% mutate(trial = as.numeric(trial)), by=c("subject", "trial")) %>% select(-valid)

# Baseline Validation
baseline.test.validation = validateBaselines(fixations.test.valid, messages, exclusions, maxDeviation_rel, maxSpread, saveBaselinePlots, postfix="test")
baseline.test.summary <- baseline.test.validation$baseline.summary
baseline.test.trials <- baseline.test.validation$baseline.trials
baseline.test.trials <- baseline.test.trials %>% 
  left_join(baseline.test.summary %>% select(subject, invalid) %>% mutate(subject = subject %>% as.numeric())) %>% 
  mutate(blok = ifelse(invalid > outlierLimit.eye, FALSE, blok))

baseline.test.summary = baseline.test.summary %>% 
  mutate(included = invalid <= outlierLimit.eye & range_x <= maxSpread & range_y <= maxSpread)
with(baseline.test.summary, hist(invalid, breaks=max(ntrials), main=paste0("Valid Baselines"))); abline(v = outlierLimit.eye, col="red", lwd=2, lty=2)


### FIXATIONS
baseline.test.summary %>% summarise(totalN = n(), includedN = sum(included), includedP = mean(included))
baseline.test.summary %>% group_by(subject) %>% summarise(invalid = mean(invalid)) %>% arrange(desc(invalid))
baseline.test.summary %>% summarise(mean_bl = mean(invalid)*100, sd_bl = sd(invalid)*100)
excluded_subjects = baseline.test.summary %>% filter(included == FALSE)
excluded_subjects = as.numeric(excluded_subjects$subject)
print(paste0("Number of excluded subjects (Test): ", length(excluded_subjects)))


### SACCADES
# Baseline OK?
saccades.test.valid <- saccades %>% 
  left_join(baseline.test.trials %>% mutate(trial = trial + acq2End) %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial")) %>% 
  filter(trial > acq2End) %>% 
  filter(!subject %in% excluded_subjects) %>% 
  ungroup

# Correction for baseline deviations 
saccades.test.valid <- saccades.test.valid %>% 
  mutate(start_x_corr = start_x + x_divergence, end_x_corr = end_x + x_divergence,
         start_y_corr = start_y + y_divergence, end_y_corr = end_y + y_divergence)
         
# Correct for picture onset
end.testtrial = 2000
saccades.test.valid = saccades.test.valid %>%
  mutate(end.testtrial = end.testtrial, # realign such that 0 = picture onset
         start_time = start_time - picOnset, end_time = end_time - picOnset, # realign such that 0 = picture onset
         start_time = ifelse(start_time < 0, 0, start_time), # discard fraction of fixation before stimulus
         end_time = ifelse(end_time > end.testtrial, end.testtrial, end_time), # discard fraction of fixation after end
         picOnset = picOnset - picOnset,
         dur = end_time - start_time) %>% 
  select(-feedbackOnset)

saccades.test.valid = saccades.test.valid %>% filter(dur > 0)

# Add ROIs
saccades.test.valid <- saccades.test.valid %>% 
  left_join(roisTest, by=c("subject", "trial"))

saccades.test.valid <- saccades.test.valid %>% 
  mutate(ROI1 = checkRoi(end_x_corr, end_y_corr, roi1.xleft, roi1.xright, roi1.ybottom, roi1.ytop),
         ROI2 = checkRoi(end_x_corr, end_y_corr, roi2.xleft, roi2.xright, roi2.ybottom, roi2.ytop),
         ROI3 = checkRoi(end_x_corr, end_y_corr, roi3.xleft, roi3.xright, roi3.ybottom, roi3.ytop),
         ROI4 = checkRoi(end_x_corr, end_y_corr, roi4.xleft, roi4.xright, roi4.ybottom, roi4.ytop)) %>% 
  ungroup

saccades.test.valid <- saccades.test.valid %>% 
  mutate(condition = ifelse(ROI1, stim1, ifelse(ROI2, stim2, ifelse(ROI3, stim3, ifelse(ROI4, stim4, NA))))) %>%
  mutate(across('condition', str_replace_all, rep_str))


# Calculate lengths of saccades
saccades.test.analysis = saccades.test.valid %>% 
  mutate(start_x_corr_cm = pixToCm(start_x_corr, screen.width, screen.width.cm),
         end_x_corr_cm = pixToCm(end_x_corr, screen.width, screen.width.cm),
         start_y_corr_cm = pixToCm(start_y_corr, screen.height, screen.height.cm),
         end_y_corr_cm = pixToCm(end_y_corr, screen.height, screen.height.cm)) %>% 
  mutate(angle = visangle(start_x_corr_cm, start_y_corr_cm, end_x_corr_cm, end_y_corr_cm, distance))

saccades.test.analysis <- saccades.test.analysis %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>% 
  mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>% 
  mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))

# Add scores and write saccades to CSV
saccades.test.analysis <- saccades.test.analysis %>%
  left_join(scores, by="subject")

write.csv2(saccades.test.analysis, file.path(path, "Gaze", "saccades_test.csv"), row.names=FALSE, quote=FALSE)

# Percentage of saccades going to the stimuli

## Trials without novel stimuli
saccades.test1.prop <- saccades.test.analysis %>%
  # filter(motivation_points > 4) %>%  # Filter for Motivation
  # filter ((trial > acq2End) & (trial <= test1End)) %>%
  mutate(condition = ifelse(angle < 1, NA, condition)) %>% 
  filter(!(str_detect(stim1, "new") | str_detect(stim2, "new") | str_detect(stim3, "new") | str_detect(stim4, "new"))) %>%
  filter(blok) %>%
  filter(!contains_blink) %>%
  mutate(condition = as.factor(condition)) %>%
  group_by(subject) %>%
  count(condition) %>%
  summarise(subject=subject, condition=condition, n=n, sum = sum(n)) %>%
  mutate(rel_freq = n/sum) %>%
  drop_na(condition) %>%
  ungroup %>%
  summarise(rel_freq = mean(rel_freq), .by=c(subject, condition))

saccades.test1.prop <- bind_rows(saccades.test1.prop, saccades.test1.prop %>% tidyr::expand(subject, condition)) %>%
  distinct(subject, condition, .keep_all=T) %>%
  mutate(rel_freq = replace_na(rel_freq, 0)) %>%
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>%
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>%
  mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>%
  mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))

saccades.test1.roi.summary <- saccades.test1.prop %>%
  summarise(Mean = mean(rel_freq), SD = sd(rel_freq), .by=c(condition, condition_social, condition_threat)) %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

saccades.test1.prop <- saccades.test1.prop %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

saccades.test1.roi.filter <- saccades.test1.prop %>%
  reframe(n_cond = length(unique(condition)), .by=c(subject))


plot_prop_sacc_familiar <- ggplot(saccades.test1.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.test1.prop, aes(y = rel_freq), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
  # labs(title = paste("Proportion of Saccades to the Stimulus (N = ", n_distinct(saccades.test.prop$subject), ")", sep=""), x = "Conditions", y = "Proportion") +
  labs(title = "Proportion of Saccades", x = NULL, y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()


anova <- saccades.test1.prop %>%
  left_join(saccades.test1.roi.filter, by=c("subject")) %>%
  filter(n_cond == 4) %>%
  group_by(subject, condition_social, condition_threat) %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
  ez::ezANOVA(dv=.(rel_freq),
              wid=.(subject),
              within=.(condition_threat, condition_social),
              # between=.(SPAI),
              detailed=T, type=3)
anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # threat * social


## Trials with novel stimuli
saccades.test2.prop <- saccades.test.analysis %>%
  # filter(motivation_points > 4) %>%  # Filter for Motivation
  # filter ((trial > test1End) & (trial <= test2End)) %>%
  mutate(condition = ifelse(angle < 1, NA, condition)) %>% 
  filter(str_detect(stim1, "new") | str_detect(stim2, "new") | str_detect(stim3, "new") | str_detect(stim4, "new")) %>% 
  filter(blok) %>%
  filter(!contains_blink) %>%
  mutate(condition = as.factor(condition)) %>%
  group_by(subject) %>%
  count(condition) %>%
  summarise(subject=subject, condition=condition, n=n, sum = sum(n)) %>%
  mutate(rel_freq = n/sum) %>%
  drop_na(condition) %>%
  ungroup %>%
  summarise(rel_freq = mean(rel_freq), .by=c(subject, condition))

saccades.test2.prop <- bind_rows(saccades.test2.prop, saccades.test2.prop %>% tidyr::expand(subject, condition)) %>%
  distinct(subject, condition, .keep_all=T) %>%
  mutate(rel_freq = replace_na(rel_freq, 0)) %>%
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>%
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>%
  mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>%
  mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))

saccades.test2.roi.summary <- saccades.test2.prop %>%
  summarise(Mean = mean(rel_freq), SD = sd(rel_freq), .by=c(condition, condition_social, condition_threat, condition_novelty)) %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

saccades.test2.prop <- saccades.test2.prop %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

saccades.test2.roi.filter <- saccades.test2.prop %>%
  reframe(n_cond = length(unique(condition)), .by=c(subject))

plot_prop_sacc_novel <- ggplot(saccades.test2.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.test2.prop, aes(y = rel_freq), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
  # labs(title = paste("Proportion of Saccades to the Stimulus (N = ", n_distinct(saccades.test.prop$subject), ")", sep=""), x = "Conditions", y = "Proportion") +
  labs(title = "Proportion of Saccades", x = NULL, y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

plot_prop_sacc_novel <- plot_prop_sacc_novel + facet_grid(cols = vars(condition_novelty))
ggsave(file.path(path, "Plots", "Gaze", "test", "all_saccades_proportion_roi.png"), width=2800, height=2000, units="px")


anova <- saccades.test2.prop %>%
  left_join(saccades.test2.roi.filter, by=c("subject")) %>%
  filter(n_cond == 4) %>%
  group_by(subject, condition_social, condition_threat, condition_novelty) %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat), condition_novelty = as.factor(condition_novelty)) %>%
  ez::ezANOVA(dv=.(rel_freq),
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

# # With Microsaccades
# saccades.test.prop <- saccades.test.analysis %>%
#   filter(blok) %>%
#   filter(!contains_blink) %>%
#   mutate(condition = as.factor(condition)) %>% 
#   group_by(subject, trial) %>% 
#   count(condition) %>% 
#   summarise(subject=subject, trial=trial, condition=condition, n=n, sum = sum(n)) %>% 
#   mutate(rel_freq = n/sum) %>% 
#   drop_na(condition) %>% 
#   ungroup %>%
#   summarise(rel_freq = mean(rel_freq), .by=c(subject, condition))
# 
# saccades.test.prop <- bind_rows(saccades.test.prop, saccades.test.prop %>% expand(subject, condition)) %>% 
#   distinct(subject, condition, .keep_all=T) %>% 
#   mutate(rel_freq = replace_na(rel_freq, 0)) %>% 
#   mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
#   mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>% 
#   mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>% 
#   mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))
# 
# 
# saccades.test.roi.summary <- saccades.test.prop %>%
#   summarise(Mean = mean(rel_freq), SD = sd(rel_freq), .by=c(condition, condition_social, condition_threat, condition_novelty)) %>% 
#   mutate(condition = str_remove(condition, ",\nnew"))
# 
# saccades.test.prop <- saccades.test.prop %>% 
#   mutate(condition = str_remove(condition, ",\nnew"))
# 
# saccades.test.roi.filter <- saccades.test.prop %>%
#   reframe(n_cond = length(unique(condition)), .by=c(subject))
# 
# plot_prop_sacc <- ggplot(saccades.test.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
#   geom_col(position = "dodge", width = 0.7) +
#   geom_errorbar(
#     aes(ymin = Mean - SD, ymax = Mean + SD),
#     position = position_dodge(width = 0.7),
#     width = 0.25) +
#   geom_point(data=saccades.test.prop, aes(y = rel_freq), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
#   # labs(title = paste("Proportion of Saccades to the Stimulus (N = ", n_distinct(saccades.test.prop$subject), ")", sep=""), x = "Conditions", y = "Proportion") +
#   labs(title = "Proportion of Saccades to the Stimuli", x = NULL, y = "Proportion") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d()
# 
# plot_prop_sacc <- plot_prop_sacc + facet_grid(cols = vars(condition_novelty))
# ggsave(file.path(path, "Plots", "Gaze", "test", "all_saccades_proportion_roi.png"), width=2800, height=2000, units="px")
# 
# 
# saccades.test.prop %>%
#   left_join(saccades.test.roi.filter, by=c("subject")) %>%
#   filter(n_cond == 4) %>%
#   group_by(subject, condition_social, condition_threat, condition_novelty) %>%
#   mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat), condition_novelty = as.factor(condition_novelty)) %>%
#   ez::ezANOVA(dv=.(rel_freq),
#               wid=.(subject),
#               within=.(condition_social, condition_threat, condition_novelty),
#               # between=.(SPAI),
#               detailed=T, type=3) %>%
#   apa::anova_apa()

# # Percentage of first saccades
# saccades.test.prop <- saccades.test.analysis %>%
#   mutate(condition = ifelse(angle < 1, NA, condition)) %>%
#   filter(blok) %>%
#   filter(!contains_blink) %>%
#   mutate(condition = as.factor(condition)) %>%
#   group_by(subject, trial) %>%
#   drop_na(condition) %>%
#   mutate(new_trial = ifelse((trial != lag(trial)), TRUE, FALSE)) %>%
#   mutate(new_trial = ifelse(is.na(new_trial), TRUE, new_trial)) %>%
#   filter(new_trial) %>%
#   select(-new_trial) %>%
#   ungroup %>%
#   group_by(subject) %>%
#   count(condition) %>%
#   summarise(subject=subject, condition=condition, n=n, sum = test2End-acq2End) %>%
#   mutate(rel_freq = n/(sum)) %>%
#   ungroup %>%
#   summarise(rel_freq = mean(rel_freq), .by=c(subject, condition))
# 
# saccades.test.prop <- bind_rows(saccades.test.prop, saccades.test.prop %>% tidyr::expand(subject, condition)) %>%
#   distinct(subject, condition, .keep_all=T) %>%
#   mutate(rel_freq = replace_na(rel_freq, 0)) %>%
#   mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>%
#   mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>%
#   mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>%
#   mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))
# 
# # saccades.test.prop.wide <- saccades.test.prop %>%
# #   select(c(subject, condition, rel_freq)) %>% 
# #   pivot_wider(names_from = condition, values_from = rel_freq)
# # 
# # saccades.test.prop.wide <- saccades.test.prop.wide %>%
# #   left_join(scores, by="subject") %>% 
# #   rename_with(~ gsub(",\n", "_", .x, fixed = TRUE))
# # 
# # write.csv2(saccades.test.prop.wide, file.path(path, "first_sacc_prop_test_wide.csv"), row.names=FALSE, quote=FALSE)
# 
# saccades.test.roi.summary <- saccades.test.prop %>%
#   summarise(Mean = mean(rel_freq), SD = sd(rel_freq), .by=c(condition, condition_social, condition_threat, condition_novelty)) %>% 
#   mutate(condition = str_remove(condition, ",\nnew"))
# 
# saccades.test.prop <- saccades.test.prop %>% 
#   mutate(condition = str_remove(condition, ",\nnew"))
# 
# saccades.test.roi.filter <- saccades.test.prop %>%
#   reframe(n_cond = length(unique(condition)), .by=c(subject))
# 
# plot_prop_first_sacc <- ggplot(saccades.test.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
#   geom_col(position = "dodge", width = 0.7) +
#   geom_errorbar(
#     aes(ymin = Mean - SD, ymax = Mean + SD),
#     position = position_dodge(width = 0.7),
#     width = 0.25) +
#   geom_point(data=saccades.test.prop, aes(y = rel_freq), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
#   # labs(title = paste("Proportion of Saccades to the Stimulus (N = ", n_distinct(saccades.test.prop$subject), ")", sep=""), x = "Conditions", y = "Proportion") +
#   labs(title = "Proportion of First Saccades", x = NULL, y = "Proportion") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_fill_viridis_d() +
#   scale_color_viridis_d()
# 
# plot_prop_first_sacc <- plot_prop_first_sacc + facet_grid(cols = vars(condition_novelty))
# ggsave(file.path(path, "Plots", "Gaze", "test", "first_saccades_proportion_roi.png"), width=2800, height=2000, units="px")
# 
# saccades.test.prop %>%
#   left_join(saccades.test.roi.filter, by=c("subject")) %>%
#   filter(n_cond == 4) %>%
#   group_by(subject, condition_social, condition_threat, condition_novelty) %>%
#   mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat), condition_novelty = as.factor(condition_novelty)) %>%
#   ez::ezANOVA(dv=.(rel_freq),
#               wid=.(subject),
#               within=.(condition_social, condition_threat, condition_novelty),
#               # between=.(SPAI),
#               detailed=T, type=3) %>%
#   apa::anova_apa()

# First Saccades Novelty

## Trials without novel stimuli
saccades.test1.prop <- saccades.test.analysis %>% 
  # filter(motivation_points > 4) %>%  # Filter for Motivation
  # filter ((trial > acq2End) & (trial <= test1End)) %>%
  filter(!(str_detect(stim1, "new") | str_detect(stim2, "new") | str_detect(stim3, "new") | str_detect(stim4, "new"))) %>%
  mutate(condition = ifelse(angle < 1, NA, condition)) %>% 
  filter(blok) %>%
  filter(!contains_blink) %>%
  mutate(condition = as.factor(condition)) %>% 
  group_by(subject, trial) %>% 
  drop_na(condition) %>% 
  mutate(new_trial = ifelse((trial != lag(trial)), TRUE, FALSE)) %>% 
  mutate(new_trial = ifelse(is.na(new_trial), TRUE, new_trial)) %>% 
  filter(new_trial) %>% 
  select(-new_trial) %>%
  ungroup %>% 
  group_by(subject) %>% 
  count(condition) %>% 
  summarise(subject=subject, condition=condition, n=n, sum = test1End - acq2End) %>% 
  mutate(rel_freq = n/(sum)) %>% 
  ungroup %>%
  summarise(rel_freq = mean(rel_freq), .by=c(subject, condition))

saccades.test1.prop <- bind_rows(saccades.test1.prop, saccades.test1.prop %>% tidyr::expand(subject, condition)) %>% 
  distinct(subject, condition, .keep_all=T) %>% 
  mutate(rel_freq = replace_na(rel_freq, 0)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>% 
  mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>% 
  mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))

saccades.test1.roi.filter <- saccades.test1.prop %>%
  reframe(n_cond = length(unique(condition)), .by=c(subject))

saccades.test1.roi.summary <- saccades.test1.prop %>%
  summarise(Mean = mean(rel_freq), SD = sd(rel_freq), .by=c(condition, condition_social, condition_threat)) %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

saccades.test1.prop <- saccades.test1.prop %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

plot_prop_first_sacc_familiar <- ggplot(saccades.test1.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.test1.prop, aes(y = rel_freq), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
  # labs(title = paste("Proportion of Saccades to the Stimulus (N = ", n_distinct(saccades.test.prop$subject), ")", sep=""), x = "Conditions", y = "Proportion") +
  labs(title = "Proportion of First Saccades", x = NULL, y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

anova <- saccades.test1.prop %>%
  left_join(saccades.test1.roi.filter, by=c("subject")) %>%
  # filter(n_cond == 8) %>%
  group_by(subject, condition_social, condition_threat) %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
  ez::ezANOVA(dv=.(rel_freq),
              wid=.(subject),
              within=.(condition_threat, condition_social),
              # between=.(SPAI),
              detailed=T, type=3)
anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # threat * social

## Trials with novel stimuli
saccades.test2.prop <- saccades.test.analysis %>% 
  # filter(motivation_points > 4) %>%  # Filter for Motivation
  # filter ((trial > test1End) & (trial <= test2End)) %>%
  filter(str_detect(stim1, "new") | str_detect(stim2, "new") | str_detect(stim3, "new") | str_detect(stim4, "new")) %>% 
  mutate(condition = ifelse(angle < 1, NA, condition)) %>% 
  filter(blok) %>%
  filter(!contains_blink) %>%
  mutate(condition = as.factor(condition)) %>% 
  group_by(subject, trial) %>% 
  drop_na(condition) %>% 
  mutate(new_trial = ifelse((trial != lag(trial)), TRUE, FALSE)) %>% 
  mutate(new_trial = ifelse(is.na(new_trial), TRUE, new_trial)) %>% 
  filter(new_trial) %>% 
  select(-new_trial) %>%
  ungroup %>% 
  group_by(subject) %>% 
  count(condition) %>% 
  summarise(subject=subject, condition=condition, n=n, sum = 40) %>% 
  mutate(rel_freq = n/(sum)) %>% 
  ungroup %>%
  summarise(rel_freq = mean(rel_freq), .by=c(subject, condition))

saccades.test2.prop <- bind_rows(saccades.test2.prop, saccades.test2.prop %>% tidyr::expand(subject, condition)) %>% 
  distinct(subject, condition, .keep_all=T) %>% 
  mutate(rel_freq = replace_na(rel_freq, 0)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>% 
  mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>% 
  mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))

# saccades.test.prop.wide <- saccades.test.prop %>%
#   select(c(subject, condition, rel_freq)) %>%
#   pivot_wider(names_from = condition, values_from = rel_freq)
# 
# saccades.test.prop.wide <- saccades.test.prop.wide %>%
#   left_join(scores, by="subject") %>%
#   rename_with(~ gsub(",\n", "_", .x, fixed = TRUE))
# 
# write.csv2(saccades.test.prop.wide, file.path(path, "first_sacc_prop_test_wide_novelty.csv"), row.names=FALSE, quote=FALSE)

saccades.test2.roi.filter <- saccades.test2.prop %>%
  reframe(n_cond = length(unique(condition)), .by=c(subject))

saccades.test2.roi.summary <- saccades.test2.prop %>%
  summarise(Mean = mean(rel_freq), SD = sd(rel_freq), .by=c(condition, condition_social, condition_threat, condition_novelty)) %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

saccades.test2.prop <- saccades.test2.prop %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

plot_prop_first_sacc_nov <- ggplot(saccades.test2.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.test2.prop, aes(y = rel_freq), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
  # labs(title = paste("Proportion of Saccades to the Stimulus (N = ", n_distinct(saccades.test.prop$subject), ")", sep=""), x = "Conditions", y = "Proportion") +
  labs(title = "Proportion of First Saccades", x = NULL, y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

plot_prop_first_sacc_nov <- plot_prop_first_sacc_nov + facet_grid(cols = vars(condition_novelty))
ggsave(file.path(path, "Plots", "Gaze", "test", "first_saccades_proportion_roi_nov.png"), width=2800, height=2000, units="px")

anova <- saccades.test2.prop %>%
  left_join(saccades.test2.roi.filter, by=c("subject")) %>%
  # filter(n_cond == 8) %>%
  group_by(subject, condition_social, condition_threat, condition_novelty) %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat), condition_novelty = as.factor(condition_novelty)) %>%
  ez::ezANOVA(dv=.(rel_freq),
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


# Latency to first saccade
saccades.test.lat.roi <- saccades.test.analysis %>%
  # filter(motivation_points > 4) %>%  # Filter for Motivation
  filter(blok) %>%
  filter(!contains_blink) %>%
  distinct(subject, trial, condition, .keep_all=T) %>% 
  filter(angle >= 1) %>%
  summarise(latency = mean(start_time), .by=c(subject, SPAI, condition, condition_social, condition_threat, condition_novelty))%>% 
  drop_na(condition)

saccades.test.lat.roi.filter <- saccades.test.lat.roi%>%
  reframe(n_cond = length(unique(condition)), .by=c(subject))

saccades.test.lat.roi <- saccades.test.lat.roi%>%
  complete(subject, condition) %>% 
  mutate(latency = replace_na(latency, 2000)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>% 
  mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>% 
  mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))

# saccades.test.lat.roi.wide <- saccades.test.lat.roi %>%
#   select(c(subject, condition, latency)) %>%
#   pivot_wider(names_from = condition, values_from = latency)
# 
# saccades.test.lat.roi.wide <- saccades.test.lat.roi.wide %>%
#   left_join(scores, by="subject") %>%
#   rename_with(~ gsub(",\n", "_", .x, fixed = TRUE))
# 
# write.csv2(saccades.test.lat.roi.wide, file.path(path, "first_sacc_lat_test_wide_novelty.csv"), row.names=FALSE, quote=FALSE)

saccades.test.lat.roi.summary <- saccades.test.lat.roi %>%
  summarise(Mean = mean(latency), SD = sd(latency), .by=c(condition, condition_social, condition_threat, condition_novelty)) %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

saccades.test.lat.roi <- saccades.test.lat.roi %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

plot_lat_sacc <- ggplot(saccades.test.lat.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.test.lat.roi, aes(y = latency), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
  # labs(title = paste("Latency to First Saccade to the Stimulus (N = ", n_distinct(saccades.test.lat.roi$subject), ")", sep=""), x = "Conditions", y = "Latency [ms]") +
  labs(title = "Latency to First Saccade", x = NULL, y = "Latency [ms]") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

plot_lat_sacc <- plot_lat_sacc + facet_grid(cols = vars(condition_novelty))
ggsave(file.path(path, "Plots", "Gaze", "test", "saccades_latency_test.png"), width=2800, height=2000, units="px")

anova <- saccades.test.lat.roi %>%
  # # Either filter for complete cases, default: replace missing conditions with max. latency (2 s)
  # left_join(saccades.test.lat.roi.filter, by=c("subject")) %>%
  # filter(n_cond == 8) %>% 
  group_by(subject, condition_social, condition_threat, condition_novelty) %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat), condition_novelty = as.factor(condition_novelty)) %>%
  ez::ezANOVA(dv=.(latency),
              wid=.(subject),
              within=.(condition_threat, condition_social, condition_novelty),
              # between_covariates=.(SPAI), observed=.(SPAI),
              detailed=T, type=3)
anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # novelty
print(anova$ANOVA[5,] %>% partial_eta_squared_ci()) # threat * social
print(anova$ANOVA[6,] %>% partial_eta_squared_ci()) # threat * novelty
print(anova$ANOVA[7,] %>% partial_eta_squared_ci()) # social * novelty
print(anova$ANOVA[8,] %>% partial_eta_squared_ci()) # threat * social * novelty


saccades.test.lat.roi %>%
  filter(condition_social == "social") %>%
  apa::t_test(data=., latency ~ condition_threat, paired = T, alternative = "two.sided") %>% 
  apa::t_apa()

saccades.test.lat.roi %>%
  filter(condition_threat == "neg") %>%
  apa::t_test(data=., latency ~ condition_social, paired = T, alternative = "two.sided") %>% 
  apa::t_apa()

# ezDesign(data=saccades.test.lat.roi, x=condition_social, y=condition_threat, col=condition_novelty)

### FIXATIONS
fixations.test.valid = eye.valid.trial %>% filter(valid > validFixTime.trial) %>% 
  left_join(fixations %>% mutate(trial = as.numeric(trial)), by=c("subject", "trial")) %>% select(-valid)
# Baseline OK?
fixations.test.valid <- fixations.test.valid %>% 
  left_join(baseline.test.trials %>% mutate(trial = trial + acq2End) %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial")) %>% 
  filter(trial > acq2End) %>% 
  filter(!subject %in% excluded_subjects) %>% 
  ungroup

# Correction for baseline deviations 
fixations.test.valid <- fixations.test.valid %>% 
  mutate(x_corr = x + x_divergence,
         y_corr = y + y_divergence)

# Correct for picture onset
end.testtrial = 2000
fixations.test.valid = fixations.test.valid %>%
  mutate(end.testtrial = end.testtrial, # realign such that 0 = picture onset
         start = start - picOnset, end = end - picOnset, # realign such that 0 = picture onset
         start = ifelse(start < 0, 0, start), # discard fraction of fixation before stimulus
         end = ifelse(end > end.testtrial, end.testtrial, end), # discard fraction of fixation after end
         picOnset = picOnset - picOnset,
         dur = end - start) %>% 
  select(-feedbackOnset)

fixations.test.analysis = fixations.test.valid %>% filter(dur > 0) # drop saccades before picture onset

# ROIs
fixations.test.analysis <- fixations.test.analysis %>% 
  left_join(roisTest, by=c("subject", "trial"))

fixations.test.analysis <- fixations.test.analysis %>% 
  mutate(ROI1 = checkRoi(x_corr, y_corr, roi1.xleft, roi1.xright, roi1.ybottom, roi1.ytop),
         ROI2 = checkRoi(x_corr, y_corr, roi2.xleft, roi2.xright, roi2.ybottom, roi2.ytop),
         ROI3 = checkRoi(x_corr, y_corr, roi3.xleft, roi3.xright, roi3.ybottom, roi3.ytop),
         ROI4 = checkRoi(x_corr, y_corr, roi4.xleft, roi4.xright, roi4.ybottom, roi4.ytop)) %>% 
  ungroup

fixations.test.analysis <- fixations.test.analysis %>% 
  mutate(condition = ifelse(ROI1, stim1, ifelse(ROI2, stim2, ifelse(ROI3, stim3, ifelse(ROI4, stim4, NA))))) %>%
  mutate(across('condition', str_replace_all, rep_str))

# Add scores and write saccades to CSV
fixations.test.analysis <- fixations.test.analysis %>%
  left_join(scores, by="subject") %>% 
  ungroup

write.csv2(fixations.test.analysis, file.path(path, "Gaze", "fixations_test.csv"), row.names=FALSE, quote=FALSE)

fixations.test.analysis <- fixations.test.analysis %>%
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>% 
  mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>% 
  mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))

# Mean Dwell Time on stimuli
fixations.test.dwell <- fixations.test.analysis %>%
  # filter(motivation_points > 4) %>%  # Filter for Motivation
  filter(blok) %>%
  filter(!is.na(condition)) %>% 
  summarise(mean_dwell_time = mean(dur), .by=c(subject, SPAI, condition, condition_social, condition_threat, condition_novelty))

fixations.test.dwell.filter <- fixations.test.dwell %>%
  reframe(n_cond = length(unique(condition)), .by=c(subject))

fixations.test.dwell <- fixations.test.dwell %>%
  complete(subject, condition) %>% 
  mutate(mean_dwell_time = replace_na(mean_dwell_time, 0)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg")) %>% 
  mutate(condition_novelty = if_else(str_detect(condition, "new"), "novel", "familiar")) %>% 
  mutate(condition_novelty = factor(condition_novelty, levels=c("familiar", "novel")))

fixations.test.roi.summary <- fixations.test.dwell %>%
  summarise(Mean = mean(mean_dwell_time), SD = sd(mean_dwell_time), .by=c(condition, condition_social, condition_threat, condition_novelty)) %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

fixations.test.dwell <- fixations.test.dwell %>% 
  mutate(condition = str_remove(condition, ",\nnew"))

plot_fix_test <- ggplot(fixations.test.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=fixations.test.dwell, aes(y = mean_dwell_time), size = 2, shape = 21, color = "black", alpha=0.1, position = position_jitter(width=0.2, height=0.005)) +
  # labs(title = paste("Dwell Times on the Stimulus (N = ", n_distinct(fixations.test.dwell$subject), ")", sep=""), x = "Conditions", y = "Dwell time [ms]") +
  labs(title = "Dwell Times", x = NULL, y = "Dwell time [ms]") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

plot_fix_test <- plot_fix_test + facet_grid(cols = vars(condition_novelty))
ggsave(file.path(path, "Plots", "Gaze", "test", "dwell_time_test.png"), width=2800, height=2000, units="px")

anova <- fixations.test.dwell %>%
  # # Filter for complete cases, default: replace missing conditions with min. dwell time (0 s)
  # left_join(fixations.test.dwell.filter, by=c("subject")) %>%
  # filter(n_cond == 8) %>%
  group_by(subject, condition_social, condition_threat, condition_novelty) %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat), condition_novelty = as.factor(condition_novelty)) %>%
  ez::ezANOVA(dv=.(mean_dwell_time),
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

# ezDesign(data=fixations.test.dwell, x=condition_social, y=condition_threat, col=condition_novelty)

prop_sacc_novel <- plot_grid(plot_prop_sacc_novel, plot_prop_first_sacc_nov, ncol = 2, labels=c("A", "B"), align="vh")
ggsave(file.path(path, "Plots", "Gaze", "test", "prop_sacc_novel.png"), width=3500, height=1200, units="px")

prop_sacc_familiar <- plot_grid(plot_prop_sacc_familiar, plot_prop_first_sacc_familiar, ncol = 2, labels=c("A", "B"), align="vh")
ggsave(file.path(path, "Plots", "Gaze", "test", "prop_sacc_familiar.png"), width=3500, height=1200, units="px")

lat_dwell_plot <- plot_grid(plot_lat_sacc, plot_fix_test, ncol = 2, labels=c("A", "B"), align="vh")
ggsave(file.path(path, "Plots", "Gaze", "test", "lat_dwell_test.png"), width=3500, height=1200, units="px")
