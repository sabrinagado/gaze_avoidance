###############################################################################
# Gaze Contingent Avoidance Project
# Sabrina Gado & Yannik Stegmann
# Code adapted from Mario Reutter & Janna Teigeler
###############################################################################

theme_set(theme_minimal(base_size = 16))

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

theme_set(theme_minimal(base_size = 16))

{ # Variables ---------------------------------------------------------------
  # a priori exclusions for all variables
  exclusions = c() %>% unique() %>% sort()
  
  acq1End = 32 # no pain during first 20 trials
  acq2End = acq1End + 80 # acquisition for 80 trials
  testEnd = acq2End + 40
  # trials.n = 112 # number of trials that shall be analyzed (if more trials, last ones will be taken)
  
  sample.rate = 1000 #samples/second
  
  screen.height = 1080 #height of screen in pix
  screen.width  = 1920 # width of screen in pix
  
  image.acq.height = 450 # height of picture in acq phase in px
  image.acq.width = 426 # width of picture in acq phase in px
  
  image.test.height = 450 * 1.2 # height of picture in test phase in px
  image.test.width = 426 * 1.2 # width of picture in test phase in px
  
  vps = vps <- seq(1, 14)
  exclusions.eye.num = c() %>% c(exclusions) #a priori exclusions, e.g. calibration not successful
  # numToEye = function(nums) return(nums %>% formatC(width=2, format="d", flag="0") %>% paste0("gca_", .)) #add leading zeros and prefix "vp"
  # exclusions.eye = exclusions.eye.num %>% numToEye()
  vps <- vps[!(vps %in% exclusions.eye.num)]
  
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
  #load behavioral data (logs)
  path.logs = file.path("Study 1", "Experiment", "gaze_avoidance_task", "data")
  files.log.prefix = "gca_avoidance_task_"
  files.log.extension = ".log"
  
  files.rating.prefix = "gca_avoidance_task_"
  files.rating.extension = ".csv"
  
  path.eye = file.path("Study 1", "Experiment", "gaze_avoidance_task", "data") #eye tracking data
  path.pupil = file.path("Study 1", "Experiment", "gaze_avoidance_task", "data") #pupil data
  path.plots = file.path("Study 1", "Plots", "Gaze")
  
  # Get Scores
  path.scores = file.path("Study 1", "Questionnaires", "demo_scores.csv")
  scores = read_delim(path.scores, delim=";", locale=locale(decimal_mark=","), na=".", show_col_types=F)
  scores$subject <- scores$VP
  scores <- scores %>%
    select(subject, gender, age, digitimer, temperature, humidity, SPAI, SIAS, STAI_T, UI, motivation, tiredness)
  
  # Files
  files.rating = list.files(path.logs, pattern=paste0("^", files.rating.prefix, ".*", files.rating.extension, "$"), full.names=TRUE)
  files.log = list.files(path.logs, pattern=paste0("^", files.log.prefix, ".*", files.log.extension, "$"), full.names=TRUE)
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
  
  loadConditions = function(files.log.path) {
    conds <- data.frame(
      subject = integer(0),
      trial = integer(0),
      phase = character(0),
      condition = character(0),
      look = numeric(0),
      outcome = character()
    )
    
    cat("\n", "Read in subject csv files for conditions:", "\n")
    pb <- txtProgressBar(min = 1, max = length(files.log.path), style = 3)
    i = 1
    for (file.log.path in files.log.path) {
      # file.log.path = files.log[1]
      log = read_delim(file.log.path, delim="\t", col_names=F, locale=locale(decimal_mark=","), na=".", show_col_types=F)
      names(log) = c("time","type","log_message")
      log = log %>% filter(type == "INFO ")
      log <- log %>%
        group_by(log_message) %>%
        mutate(trial = row_number()) %>%
        ungroup()
      log <- log %>% 
        mutate(trial = ifelse(log_message == "FixationCross", trial, NA)) %>% 
        fill(trial)
      log <- log %>%
        slice(min(which(log_message == "FixationCross")):n())
      log <- log %>%
        filter(!grepl("Rating", log_message)) %>% 
        filter(!grepl("saved", log_message) )%>% 
        filter(!grepl("Fixation", log_message))
      log <- log %>%
        mutate(condition = log_message %>% sub("ImageOnset_", "", .)) %>% 
        mutate(condition = if_else(str_detect(log_message, "Feedback"), log_message, condition) %>% sub("FeedbackOnset_", "", .)) %>% 
        mutate(condition = if_else(str_detect(log_message, "Test"), log_message, condition) %>% sub("TestImageOnset_", "", .)) %>% 
        mutate(event = "condition") %>%
        mutate(event = if_else(str_detect(log_message, "Feedback"), "outcome", event)) %>% 
        mutate(phase = "acquisition") %>%
        mutate(phase = if_else(str_detect(log_message, "Test"), "test", phase))
      log <- log %>%
        pivot_wider(names_from = event, values_from = condition) %>% 
        fill(outcome, .direction = "up")
      log <- log %>%
        filter(!grepl("Feedback", log_message))
      log <- log %>% 
        mutate(look = if_else(str_detect(outcome, "lookedat"), 1, if_else(str_detect(outcome, "no_look"), 0, NA)))
      log <- log %>% 
        mutate(subject = file.log.path %>% sub(".*gca_avoidance_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer())
      log <- log %>%
        select(subject, trial, phase, condition, look)
      log <- log %>% 
        mutate(outcome = ifelse(look == 1 & str_detect(condition, "plus"), "shock", ifelse(look == 1 & str_detect(condition, "minus"), "reward", NA)))
      
      conds <- rbind(conds, log)
      
      setTxtProgressBar(pb, i)
      i = i+1
    }
    close(pb)
    
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
      mutate(subject = subject %>% sub("gca_avoidance_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer(),
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
      mutate(subject = subject %>% sub("gca_avoidance_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer(),
             start_y = screen.height - start_y, end_y = screen.height - end_y)
    saccades = saccades[c("subject","trial","contains_blink", "start_time", "end_time", "start_x", "start_y", "end_x", "end_y")]
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
      mutate(subject = subject %>% sub("gca_avoidance_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer())
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
    
    cat("\n", "Validating baselines: ", "\n")
    pb <- txtProgressBar(min = 1, max = vpn.n, style = 3)
    for (vpIndex in 1:vpn.n) { #TODO make inner of loop enclosed function and call it with apply?
      # vpIndex = 1
      code = vpn[vpIndex]
      
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
      
      setTxtProgressBar(pb, vpIndex)
    }
    close(pb)
    
    result_list <- list(baseline.summary = baselines.summary, baseline.trials = baselines.trial %>% select(subject, trial, condition, x, y, x_divergence, y_divergence, blok))
    return(result_list)
  }
  
  loadRois = function(files.roi.path) {
    # files.roi.path = files.rating
    roi <- data.frame(
      subject = integer(0),
      trial = integer(0),
      roi.xleft = integer(0),
      roi.xright = integer(0),
      roi.ybottom = integer(0),
      roi.ytop = integer(0),
      quadrant.xleft = integer(0),
      quadrant.xright = integer(0),
      quadrant.ybottom = integer(0),
      quadrant.ytop = integer(0)
    )
    
    cat("\n", "Read in subject csv files for ROIs: ", "\n")
    pb <- txtProgressBar(min = 1, max = length(files.roi.path), style = 3)
    i = 1
    for (file.roi.path in files.roi.path) {
      # file.roi.path = files.roi.path[1]
      roi_acq = suppressMessages(read_csv(file.roi.path, show_col_types=F))
      roi_acq <- roi_acq %>%
        select(position, trialtype)
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
        mutate(subject = file.roi.path %>% sub(".*gca_avoidance_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer())
      roi_acq <- roi_acq %>%
        select(subject, trial, roi.xleft, roi.xright, roi.ybottom, roi.ytop, quadrant.xleft, quadrant.xright, quadrant.ybottom, quadrant.ytop)
      
      roi <- rbind(roi, roi_acq)
      
      roi_test <- data.frame(
        subject = rep(file.roi.path %>% sub(".*gca_avoidance_task_", "", .) %>% sub("_20.*", "", .) %>% as.integer(), 40),
        trial = seq(113, 152),
        roi.xleft = rep(screen.width / 2 - image.test.width / 2, 40),
        roi.xright = rep(screen.width / 2 + image.test.width / 2, 40),
        roi.ytop = rep(screen.height / 2 - image.test.height / 2, 40),
        roi.ybottom = rep(screen.height / 2 + image.test.height / 2, 40),
        quadrant.xleft = NA,
        quadrant.xright = NA,
        quadrant.ytop = NA,
        quadrant.ybottom = NA)
      
      roi <- rbind(roi, roi_test)
      
      setTxtProgressBar(pb, i)
      i = i+1
    }
    close(pb)
    
    return(roi)
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

rep_str = c('cs_minus_ns'='CSpos,\nnon-social','cs_minus_s'='CSpos,\nsocial',
            'cs_plus_ns'='CSneg,\nnon-social', 'cs_plus_s'='CSneg,\nsocial')

###############################################################################
# Eye-Tracking
###############################################################################
conditions = loadConditions(files.log) %>%  # conditions
  mutate(across('condition', str_replace_all, rep_str))
write.csv2(conditions, file.path("Study 1", "Physio", "Trigger", "conditions.csv"), row.names=FALSE, quote=FALSE)

rois = loadRois(files.rating) # rois

path.eye = file.path("Study 1", "Experiment", "gaze_avoidance_task", "data") # eye tracking data
messages = loadMessages(file.path(path.eye, "messages.txt")) %>% # messages for timepoints
  mutate(across('event', str_replace_all, rep_str)) %>% 
  filter(time != 0)
write.csv2(messages, file.path("Study 1", "Gaze", "messages.csv"), row.names=FALSE, quote=FALSE)

fixations = loadFixations(file.path(path.eye, "fixations.txt"), screen.height) %>% # fixations
  left_join(conditions, by=c("subject", "trial")) %>% # get conditions & other variables
  left_join(messages %>% filter(grepl("ImageOnset", event)) %>% # get image onset
              rename(picOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  left_join(messages %>% filter(grepl("FeedbackOnset", event)) %>% # get feedback onset
              rename(feedbackOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial)
vpn.eye = fixations$subject %>% unique() %>% setdiff(exclusions.eye.num) %>% sort()

saccades = loadSaccades(file.path(path.eye, "saccades.txt"), screen.height) %>% 
  left_join(conditions, by=c("subject", "trial")) %>% # get conditions & other variables
  left_join(messages %>% filter(grepl("ImageOnset", event)) %>% # get image onset
              rename(picOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  left_join(messages %>% filter(grepl("FeedbackOnset", event)) %>% # get feedback onset
              rename(feedbackOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial)

###############################################################################
# Acquisition
###############################################################################
cat("\n", "ACQUISITION PHASE", "\n")
# determine valid acquisition trials by fixation time (invalid trials need not be evaluated by their baseline)
eye.valid.trial = fixations %>% filter(grepl("acquisition", phase)) %>% # filter for acq trials
  mutate(feedbackOnset = feedbackOnset - picOnset, # realign such that 0 = picture onset
         start = start - picOnset, end = end - picOnset, # realign such that 0 = picture onset
         start = ifelse(start < 0, 0, start), # discard fraction of fixation before stimulus
         end = ifelse(end > feedbackOnset + 200, feedbackOnset + 200, end), # discard fraction of fixation after feedback onset
         dur = end - start) %>% filter(dur > 0) %>% 
  summarise(valid = sum(dur) / mean(feedbackOnset + 200), .by=c(subject, trial))
with(eye.valid.trial, hist(valid, breaks=seq(0, 1, length.out=20+1), main=paste0("Valid Fixation Time"))); abline(v=validFixTime.trial, col="red", lwd=3, lty=2)

#exclude trials with insufficient valid fixations (need not be validated for their baseline)
fixations.acq.valid = eye.valid.trial %>% filter(valid > validFixTime.trial) %>% 
  left_join(fixations %>% mutate(trial = as.numeric(trial)), by=c("subject", "trial")) %>% select(-valid)

# Baseline Validation
baseline.acq.validation = validateBaselines(fixations.acq.valid, messages, exclusions, maxDeviation_rel, maxSpread, saveBaselinePlots, postfix="acq")
baseline.acq.summary <- baseline.acq.validation$baseline.summary
baseline.acq.trials <- baseline.acq.validation$baseline.trials
baseline.acq.trials <- baseline.acq.trials %>% 
  left_join(baseline.acq.summary %>% select(subject, invalid) %>% mutate(subject = subject %>% as.numeric()), by="subject") %>% 
  mutate(blok = ifelse(invalid > outlierLimit.eye, FALSE, blok))

baseline.acq.summary = baseline.acq.summary %>% 
  mutate(included = invalid <= outlierLimit.eye & range_x <= maxSpread & range_y <= maxSpread)
with(baseline.acq.summary, hist(invalid, breaks=max(ntrials), main=paste0("Valid Baselines"))); abline(v = outlierLimit.eye, col="red", lwd=2, lty=2)

baseline.acq.summary %>% summarise(totalN = n(), includedN = sum(included), includedP = mean(included))
baseline.acq.summary %>% summarise(invalid = mean(invalid), .by=subject) %>% arrange(desc(invalid))
print(baseline.acq.summary %>% summarise(mean_bl = mean(invalid)*100, sd_bl = sd(invalid)*100))
excluded_subjects = baseline.acq.summary %>% filter(included == FALSE)
excluded_subjects = as.numeric(excluded_subjects$subject)
cat("\n", "Number of excluded subjects (Acquisition): ", length(excluded_subjects), "\n")

### SACCADES
# Baseline OK?
saccades.acq.valid <- saccades %>% 
  left_join(baseline.acq.trials %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial")) %>% 
  filter(trial <= 112) %>% 
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
  left_join(rois, by=c("subject", "trial"))

saccades.acq.valid <- saccades.acq.valid %>% 
  mutate(ROI = checkRoi(end_x_corr, end_y_corr, roi.xleft, roi.xright, roi.ybottom, roi.ytop),
         Quadrant = checkRoi(end_x_corr, end_y_corr, quadrant.xleft, quadrant.xright, quadrant.ybottom, quadrant.ytop))

saccades.acq.valid <- saccades.acq.valid %>% 
  mutate(outcome_corr = ifelse(ROI & str_detect(condition, "neg"), "shock", ifelse(ROI & str_detect(condition, "pos"), "reward", NA)))

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
  mutate(block = ifelse(trial <= 32, 0, ifelse(trial <= 32 + 40, 1, 2))) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

saccades.acq.analysis <- saccades.acq.analysis %>% 
  mutate(correct_reaction = ifelse(str_detect(condition_threat, "neg") & str_detect(outcome_corr, "none"), TRUE, ifelse(str_detect(condition_threat, "pos") & str_detect(outcome_corr, "reward"), TRUE, FALSE)))

# Proportion of wrong outcome and invalid baselines
print(saccades.acq.analysis %>% 
  summarise(mean_outcome_trial = mean(outcomeok), .by=c(subject, trial)) %>% 
  summarise(mean_outcome_subj = (1-mean(mean_outcome_trial)) * 100, .by=c(subject)) %>% 
  summarise(mean_outcome = mean(mean_outcome_subj), sd_outcome = sd(mean_outcome_subj)))

# Add scores and write saccades to CSV
saccades.acq.analysis <- saccades.acq.analysis %>%
  left_join(scores, by="subject")

write.csv2(saccades.acq.analysis, file.path("Study 1", "Gaze", "saccades_acq.csv"), row.names=FALSE, quote=FALSE)


# Percentage of saccades going towards the stimuli
saccades.acq.prop <- saccades.acq.analysis %>%
  # filter (trial >= 32) %>%
  filter(blok) %>%
  filter(outcomeok) %>%
  filter(!contains_blink) %>%
  summarise(absolute_frequency_ROI = sum(angle >= 1 & ROI & start_time < feedbackOnset), relative_frequency_ROI = mean(angle >= 1 & ROI & start_time < feedbackOnset), relative_frequency_Quadrant = mean(Quadrant & start_time < feedbackOnset), .by=c(subject, SPAI, condition, condition_social, condition_threat))

saccades.acq.roi.summary <- saccades.acq.prop %>%
  summarise(Mean = mean(relative_frequency_ROI), SD = sd(relative_frequency_ROI), .by=c(condition, condition_social, condition_threat))

plot_saccades_acq_exp1 <- ggplot(saccades.acq.roi.summary %>% 
                                   mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                                          condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
                                   mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))), aes(x = condition_social, y = Mean, group = condition_threat, color=condition_threat, fill = condition_threat)) +
  geom_point(data = saccades.acq.prop %>%
               mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                      condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
               mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))),
             aes(x = condition_social, y = relative_frequency_ROI, group = condition_threat, color = condition_threat), 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.005, dodge.width = 0.6), alpha = 0.1, size=2, shape = 21) +
  geom_point(position = position_dodge(width = 0.6), shape = "square", size=3) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), position = position_dodge(width = 0.6), width = 0.25) +
  labs(title = paste("Proportion of Approach Trials", sep=""), x = NULL, y = "Proportion") +
  # theme_minimal() +
  theme(legend.position="right") + 
  theme(legend.title=element_blank()) +
  scale_fill_viridis_d("Condition", end = 0.15, begin = 0.85) + 
  scale_color_viridis_d("Condition", end = 0.15, begin = 0.85)

# ggsave(file.path("Study 1", "Plots", "Gaze", "acquisition", "saccades_proportion_roi.png"), width=1800, height=2000, units="px")

cat("\n\n", "Analyses of the proportion of trials with a saccade towards the stimuli in the acquisition phase", "\n")
anova <- saccades.acq.prop %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
  ez::ezANOVA(dv=.(relative_frequency_ROI),
              wid=.(subject),
              within=.(condition_threat, condition_social),
              between=.(SPAI),
              detailed=T, type=3)

anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # SPAI
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[5,] %>% partial_eta_squared_ci()) # interaction: SPAI x threat
print(anova$ANOVA[6,] %>% partial_eta_squared_ci()) # interaction: SPAI x social
print(anova$ANOVA[7,] %>% partial_eta_squared_ci()) # interaction: threat x social
print(anova$ANOVA[8,] %>% partial_eta_squared_ci()) # interaction: SPAI x threat x social

# Perform a paired samples t-test
cat("\n\n", "T-test of Social vs. Non-social in CSneg trials", "\n")
saccades.acq.prop.test <- saccades.acq.prop %>%
  filter(condition_threat == "neg") %>%
  mutate(condition_social = str_remove_all(condition_social, "-")) %>% 
  pivot_wider(id_cols = subject, names_from = condition_social, values_from = relative_frequency_ROI)
apa::t_test(saccades.acq.prop.test$social, saccades.acq.prop.test$nonsocial, alternative = "two.sided", paired=TRUE) %>% apa::t_apa()

cat("\n\n", "Correlation of proportion of trials with a saccade towards the stimuli in CSneg, social trials and SPAI", "\n")
saccades.acq.prop %>%
  filter(condition_threat == "neg") %>%
  filter(condition_social == "social") %>%
  cor.test(x=.$relative_frequency_ROI, y=.$SPAI, method="pearson", alternative = "two.sided") %>% 
  apa::cor_apa()

# Correlation Correct Avoidance and Discrimination
discrimination <- read.csv(file.path("Study 1", "Behavior", "discrimination.csv"), sep=';', dec=",")

avoidance_corr.acq.prop <- saccades.acq.analysis %>%
  filter(blok) %>%
  filter(outcomeok) %>%
  summarise(absolute_frequency_corr = sum(correct_reaction), relative_frequency_corr = mean(correct_reaction), .by=c(subject, condition_social))

avoidance_corr.acq.prop <- avoidance_corr.acq.prop %>%
  left_join(discrimination, by=c("subject", "condition_social")) %>%
  rename(Category = condition_social)

cat("\n\n", "Correlation between performance in task (correct avoidance) and discrimination performance", "\n")
cor.test(avoidance_corr.acq.prop$relative_frequency_corr, avoidance_corr.acq.prop$discrimination) %>% apa::cor_apa()

ggplot(avoidance_corr.acq.prop, aes(x = discrimination, y = relative_frequency_corr, color = Category, fill = Category)) +
  geom_point(shape = 21, size = 2, color = "black", alpha=0.5) +
  geom_smooth(method = "lm", aes(fill=Category), alpha = 0.1) +
  # labs(title = paste("Correlation of Discrimination and Avoidance Performance (N = ", n_distinct(avoidance_corr.acq.prop$subject), ")", sep=""), x = "Discrimination Score", y = "Performance Score") +
  labs(title = "Correlation of Discrimination and Performance", x = "Discrimination Score", y = "Performance Score") +
  # theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme(legend.position="top", legend.box = "horizontal", legend.title=element_blank())
ggsave(file.path("Study 1", "Plots", "Gaze", "acquisition", "corr_discrimination_avoidance.png"), width=1600, height=1600, units="px")


# Latency to first saccade going towards the stimuli
saccades.acq.lat.roi <- saccades.acq.analysis %>%
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

plot_latencies_acq_exp1 <- ggplot(saccades.acq.lat.roi.summary %>% 
                                   mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                                          condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
                                   mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))), aes(x = condition_social, y = Mean, group = condition_threat, color=condition_threat, fill = condition_threat)) +
  geom_point(data = saccades.acq.lat.roi %>%
               mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                      condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
               mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))),
             aes(x = condition_social, y = latency, group = condition_threat, color = condition_threat), 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.005, dodge.width = 0.6), alpha = 0.1, size=2, shape = 21) +
  geom_point(position = position_dodge(width = 0.6), shape = "square", size=3) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), position = position_dodge(width = 0.6), width = 0.25) +
  labs(title = paste("Latency of First Saccade", sep=""), x = NULL, y = "Latency [ms]") +
  # theme_minimal() +
  theme(legend.position="right") + 
  theme(legend.title=element_blank()) +
  scale_fill_viridis_d("Condition", end = 0.15, begin = 0.85) + 
  scale_color_viridis_d("Condition", end = 0.15, begin = 0.85)

# ggsave(file.path("Study 1", "Plots", "Gaze", "acquisition", "saccades_latency_roi.png"), width=1800, height=2000, units="px")

cat("\n\n", "Analyses of the latencies of the first saccade towards the stimuli in the acquisition phase", "\n")
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
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # interaction: threat x social

# Average lengths of saccades going towards the stimuli
saccades.acq.len.roi <- saccades.acq.analysis %>%
  # filter (trial >= 32) %>%
  filter(blok) %>%
  filter(outcomeok) %>%
  filter(!contains_blink) %>%
  # filter(angle >= 1) %>%
  filter(ROI & start_time < feedbackOnset) %>%
  reframe(length = mean(angle), .by=c(subject, SPAI, condition, condition_social, condition_threat))

saccades.acq.len.roi.filter <- saccades.acq.len.roi%>%
  reframe(n_cond = length(unique(condition)), .by=c(subject))

saccades.acq.len.roi.summary <- saccades.acq.len.roi %>%
  summarise(Mean = mean(length), SD = sd(length), .by=c(condition, condition_social, condition_threat))

plot_lengths_acq_exp1 <- ggplot(saccades.acq.len.roi.summary %>% 
                                    mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                                           condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
                                    mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))), aes(x = condition_social, y = Mean, group = condition_threat, color=condition_threat, fill = condition_threat)) +
  geom_point(data = saccades.acq.len.roi %>%
               mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                      condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
               mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))),
             aes(x = condition_social, y = length, group = condition_threat, color = condition_threat), 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.005, dodge.width = 0.6), alpha = 0.1, size=2, shape = 21) +
  geom_point(position = position_dodge(width = 0.6), shape = "square", size=3) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), position = position_dodge(width = 0.6), width = 0.25) +
  labs(title = paste("Average Length of Saccades", sep=""), x = NULL, y = "Length [degree visual angle]") +
  # theme_minimal() +
  theme(legend.position="right") + 
  theme(legend.title=element_blank()) +
  scale_fill_viridis_d("Condition", end = 0.15, begin = 0.85) + 
  scale_color_viridis_d("Condition", end = 0.15, begin = 0.85)

# ggsave(file.path("Study 1", "Plots", "Gaze", "acquisition", "saccades_latency_roi.png"), width=1800, height=2000, units="px")

cat("\n\n", "Analyses of the lengths of the saccades towards the stimuli in the acquisition phase", "\n")
anova <- saccades.acq.len.roi %>%
  left_join(saccades.acq.len.roi.filter, by=c("subject")) %>%
  filter(n_cond == 4) %>%
  group_by(subject, condition_social, condition_threat) %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
  ez::ezANOVA(dv=.(length),
              wid=.(subject),
              within=.(condition_threat, condition_social),
              # between=.(SPAI),
              detailed=T, type=3)
anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # interaction: threat x social


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
  left_join(rois, by=c("subject", "trial"))

fixations.acq.analysis <- fixations.acq.analysis %>% 
  mutate(ROI = checkRoi(x_corr,y_corr, roi.xleft, roi.xright, roi.ybottom, roi.ytop),
         Quadrant = checkRoi(x_corr, y_corr, quadrant.xleft, quadrant.xright, quadrant.ybottom, quadrant.ytop))

# Add scores and write saccades to CSV
fixations.acq.analysis <- fixations.acq.analysis %>%
  left_join(scores, by="subject") %>% 
  ungroup

write.csv2(fixations.acq.analysis, file.path("Study 1", "Gaze", "fixations_acq.csv"), row.names=FALSE, quote=FALSE)

fixations.acq.analysis <- fixations.acq.analysis %>% 
  mutate(block = ifelse(trial <= acq1End, 0, 1)) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

# Mean Dwell Time on stimuli
fixations.acq.dwell <- fixations.acq.analysis %>%
  # filter (trial >= 32) %>%
  filter(blok) %>%
  filter(ROI) %>% 
  summarise(mean_dwell_time = mean(dur), .by=c(subject, SPAI, condition, condition_social, condition_threat))

fixations.acq.roi.summary <- fixations.acq.dwell %>%
  summarise(Mean = mean(mean_dwell_time), SD = sd(mean_dwell_time), .by=c(condition, condition_social, condition_threat))

plot_fixations_acq_exp1 <- ggplot(fixations.acq.roi.summary %>% 
                                  mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                                         condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
                                  mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))), aes(x = condition_social, y = Mean, group = condition_threat, color=condition_threat, fill = condition_threat)) +
  geom_point(data = fixations.acq.dwell %>%
               mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                      condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
               mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))),
             aes(x = condition_social, y = mean_dwell_time, group = condition_threat, color = condition_threat), 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.005, dodge.width = 0.6), alpha = 0.1, size=2, shape = 21) +
  geom_point(position = position_dodge(width = 0.6), shape = "square", size=3) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), position = position_dodge(width = 0.6), width = 0.25) +
  labs(title = paste("Dwell Time", sep=""), x = NULL, y = "Dwell time [ms]") +
  # theme_minimal() +
  theme(legend.position="right") + 
  theme(legend.title=element_blank()) +
  scale_fill_viridis_d("Condition", end = 0.15, begin = 0.85) + 
  scale_color_viridis_d("Condition", end = 0.15, begin = 0.85)

# ggsave(file.path("Study 1", "Plots", "Gaze", "acquisition", "dwell_time_roi.png"), width=1800, height=2000, units="px")

cat("\n\n", "Analyses of the dwell time on the stimuli in the acquisition phase", "\n")
anova <- fixations.acq.dwell %>%
  filter(subject != 42) %>% 
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
  ez::ezANOVA(dv=.(mean_dwell_time),
              wid=.(subject),
              within=.(condition_threat, condition_social),
              # between=.(SPAI),
              detailed=T, type=3)
anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # interaction: threat x social

# attention_acq_plot <- plot_grid(plot_saccades_acq, plot_latencies_acq, plot_fixations_acq, ncol = 3, labels=c("A", "B", "C"))
# ggsave(file.path("Study 1", "Plots", "Gaze", "acquisition", paste0("gaze_acq.png")), width=3500, height=1200, units="px")


###############################################################################
# Test
###############################################################################
cat("\n", "TEST PHASE", "\n")
# Determine valid test trials by fixation time (invalid trials need not be evaluated by their baseline)
end.testtrial = 11000
eye.valid.trial = fixations  %>% filter(grepl("test", phase)) %>% # filter for test trials
  mutate(end.testtrial = end.testtrial - picOnset, # realign such that 0 = picture onset
         start = start - picOnset, end = end - picOnset, # realign such that 0 = picture onset
         start = ifelse(start < 0, 0, start), # discard fraction of fixation before stimulus
         end = ifelse(end > end.testtrial, end.testtrial, end), # discard fraction of fixation after end
         dur = end - start) %>% filter(dur > 0) %>% 
  summarise(valid = sum(dur) / mean(end.testtrial), .by=c(subject, trial))
with(eye.valid.trial, hist(valid, breaks=seq(0, 1, length.out=20+1), main=paste0("Valid Fixation Time"))); abline(v=validFixTime.trial, col="red", lwd=3, lty=2)

#exclude trials with insufficient valid fixations (need not be validated for their baseline)
fixations.test.valid = eye.valid.trial %>% filter(valid > validFixTime.trial) %>% 
  left_join(fixations %>% mutate(trial = as.numeric(trial)), by=c("subject", "trial")) %>% select(-valid)

# Baseline Validation
baseline.test.validation = validateBaselines(fixations.test.valid, messages, exclusions, maxDeviation_rel, maxSpread, saveBaselinePlots, postfix="test")
baseline.test.summary <- baseline.test.validation$baseline.summary
baseline.test.trials <- baseline.test.validation$baseline.trials
baseline.test.trials <- baseline.test.trials %>% 
  left_join(baseline.test.summary %>% select(subject, invalid) %>% mutate(subject = subject %>% as.numeric()), by="subject") %>% 
  mutate(blok = ifelse(invalid > outlierLimit.eye, FALSE, blok))

baseline.test.summary = baseline.test.summary %>% 
  mutate(included = invalid <= outlierLimit.eye & range_x <= maxSpread & range_y <= maxSpread)
with(baseline.test.summary, hist(invalid, breaks=max(ntrials), main=paste0("Valid Baselines"))); abline(v = outlierLimit.eye, col="red", lwd=2, lty=2)


### FIXATIONS
baseline.test.summary %>% summarise(totalN = n(), includedN = sum(included), includedP = mean(included))
baseline.test.summary %>% summarise(invalid = mean(invalid), .by=subject) %>% arrange(desc(invalid))
print(baseline.test.summary %>% summarise(mean_bl = mean(invalid)*100, sd_bl = sd(invalid)*100))
excluded_subjects = baseline.test.summary %>% filter(included == FALSE)
excluded_subjects = as.numeric(excluded_subjects$subject)
cat("\n", "Number of excluded subjects (Test): ", length(excluded_subjects), "\n")

# Baseline OK?
fixations.test.valid <- fixations.test.valid %>% 
  left_join(baseline.test.trials %>% mutate(trial = trial + 112) %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial")) %>% 
  filter(!subject %in% excluded_subjects)

# Correction for baseline deviations 
fixations.test.valid <- fixations.test.valid %>% 
  mutate(x_corr = x + x_divergence,
         y_corr = y + y_divergence)

# Correct for picture onset
fixations.test.valid = fixations.test.valid %>%
  mutate(end.testtrial = end.testtrial - picOnset,
         start = start - picOnset, end = end - picOnset, # realign such that 0 = picture onset
         start = ifelse(start < 0, 0, start), # discard fraction of fixation before stimulus
         end = ifelse(end > end.testtrial, end.testtrial, end), # discard fraction of fixation after end
         picOnset = picOnset - picOnset,
         dur = end - start)

fixations.test.valid = fixations.test.valid %>% filter(dur > 0)

# Add ROIs
fixations.test.valid <- fixations.test.valid %>% 
  left_join(rois, by=c("subject", "trial"))

fixations.test.valid <- fixations.test.valid %>% 
  mutate(ROI = checkRoi(x_corr, y_corr, roi.xleft, roi.xright, roi.ybottom, roi.ytop)) %>% 
  ungroup

fixations.test.valid <- fixations.test.valid %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

# # Add scores and write saccades to CSV
fixations.test.valid <- fixations.test.valid %>%
  left_join(scores, by="subject")

write.csv2(fixations.test.valid, file.path("Study 1", "Gaze", "fixations_test.csv"), row.names=FALSE, quote=FALSE)

# Mean Dwell Time off stimuli
fixations.test.dwell <- fixations.test.valid %>%
  filter(blok) %>%
  summarise(total_dur = sum(dur), off_stim_dur = sum(dur[ROI == FALSE]), .by=c(subject, SPAI, condition, condition_social, condition_threat)) %>% 
  mutate(off_stim_prop_dur = off_stim_dur / total_dur)

fixations.test.dwell.summary <- fixations.test.dwell %>%
  summarise(Mean = mean(off_stim_prop_dur), SD = sd(off_stim_prop_dur), .by=c(condition, condition_social, condition_threat))

plot_fixations_test_exp1 <- ggplot(fixations.test.dwell.summary %>% 
                                    mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                                           condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
                                    mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))), aes(x = condition_social, y = Mean, group = condition_threat, color=condition_threat, fill = condition_threat)) +
  geom_point(data = fixations.test.dwell %>%
               mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                      condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
               mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))),
             aes(x = condition_social, y = off_stim_prop_dur, group = condition_threat, color = condition_threat), 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.005, dodge.width = 0.6), alpha = 0.1, size=2, shape = 21) +
  geom_point(position = position_dodge(width = 0.6), shape = "square", size=3) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), position = position_dodge(width = 0.6), width = 0.25) +
  labs(title = paste("Proportional Dwell Time", sep=""), x = NULL, y = "Proportion") +
  # theme_minimal() +
  theme(legend.position="right") + 
  theme(legend.title=element_blank()) +
  scale_fill_viridis_d("Condition", end = 0.15, begin = 0.85) + 
  scale_color_viridis_d("Condition", end = 0.15, begin = 0.85)

# ggsave(file.path("Study 1", "Plots", "Gaze", "acquisition", "dwell_time_roi.png"), width=1800, height=2000, units="px")

cat("\n\n", "Analyses of the dwell time off the stimuli in the test phase", "\n")
anova <- fixations.test.dwell %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
  ez::ezANOVA(dv=.(off_stim_prop_dur),
              wid=.(subject),
              within=.(condition_threat, condition_social),
              # between=.(SPAI),
              detailed=T, type=3)
anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # interaction: threat x social


### SACCADES
# Baseline OK?
saccades.test.valid <- saccades %>% 
  left_join(baseline.test.trials %>% mutate(trial = trial + 112) %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial")) %>% 
  filter(trial > 112) %>% 
  filter(!subject %in% excluded_subjects)

# Correction for baseline deviations 
saccades.test.valid <- saccades.test.valid %>% 
  mutate(start_x_corr = start_x + x_divergence, end_x_corr = end_x + x_divergence,
         start_y_corr = start_y + y_divergence, end_y_corr = end_y + y_divergence)
         
# Correct for picture onset
saccades.test.valid = saccades.test.valid %>%
  mutate(end.testtrial = end.testtrial - picOnset, # realign such that 0 = picture onset
         start_time = start_time - picOnset, end_time = end_time - picOnset, # realign such that 0 = picture onset
         start_time = ifelse(start_time < 0, 0, start_time), # discard fraction of fixation before stimulus
         end_time = ifelse(end_time > end.testtrial, end.testtrial, end_time), # discard fraction of fixation after end
         picOnset = picOnset - picOnset,
         dur = end_time - start_time)

saccades.test.valid = saccades.test.valid %>% filter(dur > 0)

# Add ROIs
saccades.test.valid <- saccades.test.valid %>% 
  left_join(rois, by=c("subject", "trial"))

saccades.test.valid <- saccades.test.valid %>% 
  mutate(ROI = checkRoi(end_x_corr, end_y_corr, roi.xleft, roi.xright, roi.ybottom, roi.ytop))

# Calculate lengths of saccades
saccades.test.analysis = saccades.test.valid %>% 
  mutate(start_x_corr_cm = pixToCm(start_x_corr, screen.width, screen.width.cm),
         end_x_corr_cm = pixToCm(end_x_corr, screen.width, screen.width.cm),
         start_y_corr_cm = pixToCm(start_y_corr, screen.height, screen.height.cm),
         end_y_corr_cm = pixToCm(end_y_corr, screen.height, screen.height.cm)) %>% 
  mutate(angle = visangle(start_x_corr_cm, start_y_corr_cm, end_x_corr_cm, end_y_corr_cm, distance))

saccades.test.analysis <- saccades.test.analysis %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

# Add scores and write saccades to CSV
saccades.test.analysis <- saccades.test.analysis %>%
  left_join(scores, by="subject")

write.csv2(saccades.test.analysis, file.path("Study 1", "Gaze", "saccades_test.csv"), row.names=FALSE, quote=FALSE)

# Percentage of saccades going away from the stimuli
saccades.test.prop <- saccades.test.analysis %>%
  filter(blok) %>%
  filter(!contains_blink) %>%
  summarise(relative_frequency_ROI = mean(!ROI), .by=c(subject, condition))

saccades.test.prop <- saccades.test.analysis %>%
  filter(blok) %>%
  filter(!contains_blink) %>%
  summarise(absolute_frequency_ROI = sum(angle >= 1 & !ROI), relative_frequency_ROI = mean(angle >= 1 & !ROI), .by=c(subject, SPAI, condition, condition_social, condition_threat))

saccades.test.roi.summary <- saccades.test.prop %>%
  summarise(Mean = mean(relative_frequency_ROI), SD = sd(relative_frequency_ROI), .by=c(condition, condition_social, condition_threat))

saccades.test.roi.filter <- saccades.test.prop %>%
  reframe(n_cond = length(unique(condition)), .by=c(subject))

plot_sacc_test_exp1 <- ggplot(saccades.test.roi.summary %>% 
                                     mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                                            condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
                                     mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))), aes(x = condition_social, y = Mean, group = condition_threat, color=condition_threat, fill = condition_threat)) +
  geom_point(data = saccades.test.prop %>%
               mutate(condition_social = recode(condition_social, "non-social" = "Non-Social", "social" = "Social"),
                      condition_threat = recode(condition_threat, "neg" = "CSneg", "pos" = "CSpos")) %>% 
               mutate(condition_threat = factor(condition_threat, level=c("CSpos", "CSneg"))),
             aes(x = condition_social, y = relative_frequency_ROI, group = condition_threat, color = condition_threat), 
             position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.005, dodge.width = 0.6), alpha = 0.1, size=2, shape = 21) +
  geom_point(position = position_dodge(width = 0.6), shape = "square", size=3) +
  geom_line(position = position_dodge(width = 0.6), linewidth = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), position = position_dodge(width = 0.6), width = 0.25) +
  labs(title = paste("Prop. of Sacc. Away from Stimulus", sep=""), x = NULL, y = "Proportion") +
  # theme_minimal() +
  theme(legend.position="right") + 
  theme(legend.title=element_blank()) +
  scale_fill_viridis_d("Condition", end = 0.15, begin = 0.85) + 
  scale_color_viridis_d("Condition", end = 0.15, begin = 0.85)


cat("\n\n", "Analyses of the proportion of saccades away from the stimuli during the test trials", "\n")
anova <- saccades.test.prop %>%
  left_join(saccades.test.roi.filter, by=c("subject")) %>%
  filter(n_cond == 4) %>%
  group_by(subject, condition_social, condition_threat) %>%
  mutate(subject = as.factor(subject), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
  ez::ezANOVA(dv=.(relative_frequency_ROI),
              wid=.(subject),
              within=.(condition_threat, condition_social),
              # between=.(SPAI),
              detailed=T, type=3)
anova %>% apa::anova_apa()
print(anova$ANOVA[2,] %>% partial_eta_squared_ci()) # threat
print(anova$ANOVA[3,] %>% partial_eta_squared_ci()) # social
print(anova$ANOVA[4,] %>% partial_eta_squared_ci()) # interaction: threat x social

print(saccades.test.roi.summary)

