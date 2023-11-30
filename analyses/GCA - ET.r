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

{ #install packages needed
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
}


{ # Variables ---------------------------------------------------------------
  # a priori exclusions for all variables
  exclusions = c(
    1, 2, 3, 4, 5, 6 #old paradigm
  ) %>% unique() %>% sort()
  
  acq1End = 32 #no pain during first 20 trials
  acq2End = acq1End + 80 #acquisition for 40 trials
  testEnd = acq2End + 40
  trials.n = 112 #number of trials that shall be analyzed (if more trials, last ones will be taken)
  
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
  baseline = c(700, 1000) #Baseline in ms relative to stimulus onset; min(baseline) = start; max(baseline) = end
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
  # 
  screen.width <- 1920; screen.height <- 1080 #screen resolution (24" ASUS VG248QE)
  screen.width.cm <- 53.136; screen.height.cm <- 29.889
  pixsize_in_cm <- screen.width.cm / screen.width  # Pixel size in cm
  distance <- 560 # Distance Camera - Eye 560 mm (chin rest with remote tracking)
  # 
  # #statistical analysis
  # z.max = 2 #winsorize dependent variables to a z value of 2
  # q.max = pnorm(z.max * c(-1, 1)) #needed for DescTools::Winsorize function
  # q.max = c(0, 1) #switch off Winsorizing (comment out to switch on)
}

# Paths -------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')
path = getwd()

#load behavioral data (logs)
path.logs = file.path(path, "gaze_avoidance_task", "data")
files.log.prefix = "gca_"
files.log.extension = ".log"

files.rating.prefix = "gca_"
files.rating.extension = ".csv"

path.eye = file.path(path, "gaze_avoidance_task", "data") #eye tracking data
path.pupil = file.path(path, "gaze_avoidance_task", "data") #pupil data
path.plots = file.path(path, "plots")

# path.rois = "2 Material/Look At Me 3 Anx/Stimuli/mask/" %>% paste0(path2, .)

# Files -------------------------------------------------------------------
files.rating = list.files(path.logs, pattern=paste0("^", files.rating.prefix, ".*", files.rating.extension, "$"), full.names=TRUE)
files.log = list.files(path.logs, pattern=paste0("^", files.log.prefix, ".*", files.log.extension, "$"), full.names=TRUE)


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
        slice(which(log_message == "FixationCross"):n())
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
        mutate(subject = file.log.path %>% sub(".*gca_", "", .) %>% sub("_2023.*", "", .) %>% as.integer())
      log <- log %>%
        select(subject, trial, phase, condition, look)
      log <- log %>% 
        mutate(outcome = ifelse(look == 1 & str_detect(condition, "plus"), "shock", ifelse(look == 1 & str_detect(condition, "minus"), "reward", NA)))
      
      conds <- rbind(conds, log)
    }
    
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
      # file.rating.path = files.rating.path[1]
      rating = read_csv(file.rating.path)
      rating <- rating %>%
        select(stimtype, sliderStim.started, sliderStim.response)
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
        mutate(subject = file.rating.path %>% sub(".*gca_", "", .) %>% sub("_2023.*", "", .) %>% as.integer())
      rating <- rating %>%
        select(subject, phase, stimtype, sliderStim.response)
      names(rating) = c("subject","phase","condition", "rating")
      
      ratings <- rbind(rating, ratings)
    }
    
    return(ratings)
  }
  
  loadFixations = function(filePath, screen.height) {
    fixations = read_delim(filePath, delim="\t", col_names=F, skip=1, locale=locale(decimal_mark=","), na=".", show_col_types=F)
    names(fixations) = c("subject","trial","start","end","x","y", "eye")
    fixations = fixations %>%
      mutate(subject = subject %>% sub("gca_", "", .) %>% sub("_2023.*", "", .) %>% as.integer(),
             trial = trial %>% sub("Trial: ","", .) %>% as.numeric(),
             y = screen.height - y)
    return(fixations)
  }
  
  loadSaccades = function(filePath, screen.height) {
    saccades = read_delim(filePath, delim="\t", col_names=F, skip=1, locale=locale(decimal_mark=","), na=".", show_col_types=F)
    names(saccades) = c("subject","trial","contains_blink", "start_time", "end_time", "start_x", "start_y","end_x","end_y")
    saccades = saccades %>%
      mutate(subject = subject %>% sub("gca_", "", .) %>% sub("_2023.*", "", .) %>% as.integer(),
             trial = trial %>% sub("Trial: ","", .) %>% as.numeric(),
             start_y = screen.height - start_y, end_y = screen.height - end_y)
    return(saccades)
  }
  
  loadMessages = function(filePath) {
    messages = read_delim(filePath, delim="\t", col_names=F, skip=1, locale=locale(decimal_mark=","), na=".", show_col_types=F)
    names(messages) = c("subject", "trial","time","event")
    messages = messages %>%
      mutate(subject = subject %>% sub("gca_", "", .) %>% sub("_2023.*", "", .) %>% as.integer(),
             trial = trial %>% sub("Trial: ","", .) %>% as.numeric())
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
        
        # Determine onset (in ms)
        onset = 1000 # Stop Fixation Cross
        condition = fix.vp.trial$condition %>% unique()
    
        # Subtract onset from startamps
        fix.vp.trial$start  <- fix.vp.trial$start - onset
        fix.vp.trial$end <- fix.vp.trial$end - onset
        
        # Caluculate baseline as weighted average of fixations
        fix.vp.trial.bl <- fix.vp.trial[fix.vp.trial$end>min(baseline) & fix.vp.trial$start<max(baseline),]
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
        filename = file.path(path.plots, "baseline", sprintf("gca_%s_%s.png", code, postfix))
          
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
  
  loadRois = function(files.roi.path) {
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
    for (file.roi.path in files.roi.path) {
      # file.roi.path = files.roi.path[1]
      roi_acq = read_csv(file.roi.path)
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
        mutate(subject = file.roi.path %>% sub(".*gca_", "", .) %>% sub("_2023.*", "", .) %>% as.integer())
      roi_acq <- roi_acq %>%
        select(subject, trial, roi.xleft, roi.xright, roi.ybottom, roi.ytop, quadrant.xleft, quadrant.xright, quadrant.ybottom, quadrant.ytop)
      
      roi <- rbind(roi, roi_acq)
      
      roi_test <- data.frame(
        subject = rep(file.roi.path %>% sub(".*gca_", "", .) %>% sub("_2023.*", "", .) %>% as.integer(), 40),
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
    }
    
    return(roi)
  }
  
  checkRoi <- function(x, y, rect_xleft, rect_xright, rect_ybottom, rect_ytop) {
    return(x >= rect_xleft & x <= rect_xright & y <= rect_ybottom & y >= rect_ytop)
  }
}


conditions = loadConditions(files.log) # condition
ratings = loadRatings(files.rating, files.log) # ratings
rois = loadRois(files.rating) # rois
path.eye = file.path(path, "data") #eye tracking data
messages = loadMessages(file.path(path.eye, "messages.txt")) # %>% mutate(time = ifelse(time==0, 0, time + 1000)) #correct for duplicated pre-stimulus-baseline
fixations = loadFixations(file.path(path.eye, "fixations.txt"), screen.height) %>% 
  left_join(conditions, by=c("subject", "trial")) %>% # get conditions & other variables
  left_join(messages %>% filter(grepl("ImageOnset", event)) %>% # get image onset
              rename(picOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  left_join(messages %>% filter(grepl("FeedbackOnset", event)) %>% # get feedback onset
              rename(feedbackOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial)
saccades = loadSaccades(file.path(path.eye, "saccades.txt"), screen.height) %>% 
  left_join(conditions, by=c("subject", "trial")) %>% # get conditions & other variables
  left_join(messages %>% filter(grepl("ImageOnset", event)) %>% # get image onset
              rename(picOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  left_join(messages %>% filter(grepl("FeedbackOnset", event)) %>% # get feedback onset
              rename(feedbackOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial)
vpn.eye = fixations$subject %>% unique() %>% setdiff(exclusions.eye.num) %>% sort()

###############################################################################
# Acquisition
###############################################################################
#determine valid acquisition trials by fixation time (invalid trials need not be evaluated by their baseline)
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


### SACCADES
saccades = loadSaccades(file.path(path.eye, "saccades.txt"), screen.height) %>% 
  left_join(conditions, by=c("subject", "trial")) %>% # get conditions & other variables
  left_join(messages %>% filter(grepl("ImageOnset", event)) %>% # get image onset
              rename(picOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial) %>% 
  left_join(messages %>% filter(grepl("FeedbackOnset", event)) %>% # get feedback onset
              rename(feedbackOnset = time) %>% select(-event),
            by=c("subject", "trial")) %>% arrange(subject, trial)

saccades.acq.valid <- saccades %>% 
  left_join(baseline.acq.trials %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial"))

saccades.acq.valid = saccades.acq.valid %>%
  mutate(feedbackOnset = feedbackOnset - picOnset, # realign such that 0 = picture onset
         start_time = start_time - picOnset, end_time = end_time - picOnset, # realign such that 0 = picture onset
         start_time = ifelse(start_time < 0, 0, start_time), # discard fraction of fixation before stimulus
         picOnset = picOnset - picOnset,
         dur = end_time - start_time) %>% filter(dur > 0)

saccades.acq.valid = saccades.acq.valid %>% filter(start_time < feedbackOnset)

saccades.acq.valid <- saccades.acq.valid %>% 
  left_join(rois, by=c("subject", "trial"))

saccades.acq.valid <- saccades.acq.valid %>% 
  mutate(start_x_corr = start_x - x_divergence, end_x_corr = end_x - x_divergence,
         start_y_corr = start_y - y_divergence, end_y_corr = end_y - y_divergence,  
         ROI = checkRoi(end_x_corr, end_y_corr, roi.xleft, roi.xright, roi.ybottom, roi.ytop),
         Quadrant = checkRoi(end_x_corr, end_y_corr, quadrant.xleft, quadrant.xright, quadrant.ybottom, quadrant.ytop))

saccades.acq.valid <- saccades.acq.valid %>% 
  mutate(outcome_corr = ifelse(ROI & str_detect(condition, "plus"), "shock", ifelse(ROI & str_detect(condition, "minus"), "reward", NA)))

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

saccades.acq.analysis = saccades.acq.valid %>% 
  mutate(start_x_corr_cm = pixToCm(start_x_corr, screen.width, screen.width.cm),
         end_x_corr_cm = pixToCm(end_x_corr, screen.width, screen.width.cm),
         start_y_corr_cm = pixToCm(start_y_corr, screen.height, screen.height.cm),
         end_y_corr_cm = pixToCm(end_y_corr, screen.height, screen.height.cm)) %>% 
  mutate(angle = visangle(start_x_corr_cm, start_y_corr_cm, end_x_corr_cm, end_y_corr_cm, distance)) 

# Percentage of avoidance trials
avoidance.acq.prop <- saccades.acq.analysis %>% 
  filter (trial >= 32) %>%
  filter(blok) %>% 
  filter(outcomeok) %>%
  summarise(absolute_frequency_av = sum(outcome_corr=="none"), relative_frequency_av = mean(outcome_corr=="none"), .by=c(subject, condition)) 

avoidance.acq.prop.summary <- avoidance.acq.prop %>% 
  summarise(Mean = mean(relative_frequency_av), SD = sd(relative_frequency_av), .by=condition)

ggplot(avoidance.acq.prop.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD), #ymin = Mean - SD, 
    position = position_dodge(width = 0.7),
    width = 0.25
  ) +
  geom_line(data=avoidance.acq.prop, aes(y = relative_frequency_av, group = subject), alpha=0.2) +
  geom_point(data=avoidance.acq.prop, aes(y = relative_frequency_av), size = 2, shape = 21, color = "black", alpha=0.5) + # , position=position_jitter(width=0.05)) +
  labs(title = "Proportion of Avoidance-Trials", x = "Conditions", y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "avoidance_proportion_roi.png"), width=1800, height=2400, units="px")


# Percentage of saccades going towards the stimuli
saccades.acq.prop <- saccades.acq.analysis %>% 
  filter (trial >= 32) %>%
  filter(blok) %>% 
  filter(outcomeok) %>%
  filter(!contains_blink) %>%
  summarise(absolute_frequency_ROI = sum(ROI & start_time < feedbackOnset & angle >= 1), relative_frequency_ROI = mean(ROI & start_time < feedbackOnset & angle >= 1), relative_frequency_Quadrant = mean(Quadrant & start_time < feedbackOnset & angle >= 1), .by=c(subject, condition)) 
  
saccades.acq.roi.summary <- saccades.acq.prop %>% 
  summarise(Mean = mean(relative_frequency_ROI), SD = sd(relative_frequency_ROI), .by=condition)

ggplot(saccades.acq.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.acq.prop, aes(y = relative_frequency_ROI), size = 2, shape = 21, color = "black", alpha=0.5, position=position_jitter(width=0.05)) +
  labs(title = "Proportion of Trials with a Saccade to the Stimulus", x = "Conditions", y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "saccades_proportion_roi.png"), width=1800, height=2400, units="px")

# Percentage of saccades going towards the quadrant
saccades.acq.quad.summary <- saccades.acq.prop %>% 
  summarise(Mean = mean(relative_frequency_Quadrant), SD = sd(relative_frequency_Quadrant), .by=condition)

ggplot(saccades.acq.quad.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.acq.prop, aes(y = relative_frequency_Quadrant), size = 2, shape = 21, color = "black", alpha=0.5, position=position_jitter(width=0.05)) +
  labs(title = "Proportion of Trials with a Saccade to the Quadrant of the Stimulus", x = "Conditions", y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "saccades_proportion_quad.png"), width=1800, height=2400, units="px")

# Latency to first saccade going towards the stimuli
saccades.acq.lat.roi <- saccades.acq.analysis %>% 
  filter (trial >= 32) %>%
  filter(blok) %>% 
  filter(outcomeok) %>%
  filter(!contains_blink) %>%
  filter(ROI & start_time < feedbackOnset & angle >= 1) %>% 
  summarise(latency = start_time, .by=c(subject, condition)) 

saccades.acq.lat.roi.summary <- saccades.acq.lat.roi %>% 
  summarise(Mean = mean(latency), SD = sd(latency), .by=condition)

ggplot(saccades.acq.lat.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.acq.lat.roi, aes(y = latency), size = 2, shape = 21, color = "black", alpha=0.5, position=position_jitter(width=0.05)) +
  labs(title = "Latency to First Saccade to the Stimulus", x = "Conditions", y = "Latency [ms]") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "saccades_latency_roi.png"), width=1800, height=2400, units="px")

# Latency to first saccade going towards the quadrant of stimuli
saccades.acq.lat.quad <- saccades.acq.analysis %>% 
  filter (trial >= 32) %>%
  filter(blok) %>% 
  filter(outcomeok) %>%
  filter(!contains_blink) %>%
  filter(ROI & start_time < feedbackOnset & angle >= 1) %>% 
  summarise(latency = start_time, .by=c(subject, condition)) 

saccades.acq.lat.quad.summary <- saccades.acq.lat.quad %>% 
  summarise(Mean = mean(latency), SD = sd(latency), .by=condition)

ggplot(saccades.acq.lat.quad.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.acq.lat.quad, aes(y = latency), size = 2, shape = 21, color = "black", alpha=0.5, position=position_jitter(width=0.05)) +
  labs(title = "Latency to First Saccade to the Quadrant of the Stimulus", x = "Conditions", y = "Latency [ms]") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "saccades_latency_quad.png"), width=1800, height=2400, units="px")

# Length of saccade going towards the stimuli
saccades.acq.len.roi <- saccades.acq.analysis %>% 
  filter (trial >= 32) %>%
  filter(blok) %>% 
  filter(outcomeok) %>%
  filter(!contains_blink) %>%
  filter(ROI & start_time < feedbackOnset & angle >= 1) %>% 
  summarise(length = angle, .by=c(subject, condition)) 

saccades.acq.len.roi.summary <- saccades.acq.len.roi %>% 
  summarise(Mean = mean(length), SD = sd(length), .by=condition)

ggplot(saccades.acq.len.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.acq.len.roi, aes(y = length), size = 2, shape = 21, color = "black", alpha=0.5, position=position_jitter(width=0.05)) +
  labs(title = "Length of Saccade to the Stimulus", x = "Conditions", y = "Length [degree visual angle]") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "saccades_length_roi.png"), width=1800, height=2400, units="px")

# Latency to first saccade going towards the quadrant of stimuli
saccades.acq.len.quad <- saccades.acq.analysis %>% 
  filter (trial >= 32) %>%
  filter(blok) %>% 
  filter(outcomeok) %>%
  filter(!contains_blink) %>%
  filter(ROI & start_time < feedbackOnset & angle >= 1) %>% 
  summarise(length = angle, .by=c(subject, condition)) 

saccades.acq.len.quad.summary <- saccades.acq.len.quad %>% 
  summarise(Mean = mean(length), SD = sd(length), .by=condition)

ggplot(saccades.acq.len.quad.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.acq.len.quad, aes(y = length), size = 2, shape = 21, color = "black", alpha=0.5, position=position_jitter(width=0.05)) +
  labs(title = "Length of Saccade to the Quadrant of the Stimulus", x = "Conditions", y = "Length [degree visual angle]") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "saccades_length_quad.png"), width=1800, height=2400, units="px")

# ### FIXATIONS
# fixations.acq.valid <- fixations.acq.valid %>% 
#   left_join(baseline.acq.trials %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial"))
# 
# fixations.acq.valid = fixations.acq.valid %>%
#   mutate(feedbackOnset = feedbackOnset - picOnset, # realign such that 0 = picture onset
#          start = start - picOnset, end = end - picOnset, # realign such that 0 = picture onset
#          start = ifelse(start < 0, 0, start), # discard fraction of fixation before stimulus
#          picOnset = picOnset - picOnset,
#          dur = end - start) %>% filter(dur > 0)
# 
# fixations.acq.valid <- fixations.acq.valid %>% 
#   left_join(rois, by=c("subject", "trial"))
# 
# fixations.acq.valid <- fixations.acq.valid %>% 
#   mutate(x_corr = x - x_divergence,
#          y_corr = y - y_divergence,
#          ROI = checkRoi(x_corr, y_corr, roi.xleft, roi.xright, roi.ybottom, roi.ytop),
#          Quadrant = checkRoi(x_corr, y_corr, quadrant.xleft, quadrant.xright, quadrant.ybottom, quadrant.ytop))
# 
# fixations.acq.valid <- fixations.acq.valid %>% 
#   left_join(saccades.acq.valid %>% select("subject", "trial", "outcomeok"), by=c("subject", "trial")) %>% 
#   mutate(outcomeok = coalesce(outcomeok, TRUE))
# 
# fixations.acq.valid = fixations.acq.valid %>% filter(0 < start & start < feedbackOnset + 200)
# 
# # Proportion of fixations on the stimuli
# fixations.acq.valid %>% 
#   # filter (trial > 32) %>% 
#   filter(blok) %>% 
#   filter(outcomeok) %>%
#   group_by(subject, condition) %>%
#   summarise(relative_frequency_ROI = mean(ROI), relative_frequency_Quadrant = mean(Quadrant))
# 
# # Dwell time of fixations in respective quadrant
# fixations.acq.valid %>% 
#   # filter (trial > 32) %>% 
#   group_by(subject, condition, Quadrant) %>%
#   summarise(duration = sum(dur) / 1000)


###############################################################################
# Test
###############################################################################
# Determine valid test trials by fixation time (invalid trials need not be evaluated by their baseline)
end.testtrial = 11000
eye.valid.trial = fixations  %>% filter(grepl("test", phase)) %>% # filter for test trials
  mutate(end.testtrial = end.testtrial - picOnset, # realign such that 0 = picture onset
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

fixations.test.valid <- fixations.test.valid %>% 
  left_join(baseline.test.trials %>% mutate(trial = trial + 112) %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial"))

fixations.test.valid = fixations.test.valid %>%
  mutate(end.testtrial = end.testtrial - picOnset,
         start = start - picOnset, end = end - picOnset, # realign such that 0 = picture onset
         start = ifelse(start < 0, 0, start), # discard fraction of fixation before stimulus
         end = ifelse(end > end.testtrial, end.testtrial, end), # discard fraction of fixation after end
         picOnset = picOnset - picOnset,
         dur = end - start) %>% filter(dur > 0)

fixations.test.valid <- fixations.test.valid %>% 
  left_join(rois, by=c("subject", "trial"))

fixations.test.valid <- fixations.test.valid %>% 
  mutate(x_corr = x - x_divergence,
         y_corr = y - y_divergence,
         ROI = checkRoi(x_corr, y_corr, roi.xleft, roi.xright, roi.ybottom, roi.ytop)) %>% 
  ungroup

# Proportional Dwell Time On The Stimulus
fixations.test.dwell <- fixations.test.valid %>% 
  filter(blok) %>% 
  mutate(dwell.time.total = sum(dur), .by=c(subject, condition)) %>% 
  filter(ROI) %>%
  mutate(dwell.time.roi = sum(dur), .by=c(subject, condition)) %>% 
  # mutate(dwell.time.prop = dwell.time.roi/dwell.time.total)
  summarise(dwell.time.prop = mean(dwell.time.roi/dwell.time.total), .by=c(subject, condition)) 

fixations.test.dwell.summary <- fixations.test.dwell %>% 
  summarise(Mean = mean(dwell.time.prop), SD = sd(dwell.time.prop), .by=condition)

ggplot(fixations.test.dwell.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=fixations.test.dwell, aes(y = dwell.time.prop), size = 2, shape = 21, color = "black", alpha=0.5, position=position_jitter(width=0.05)) +
  labs(title = "Proportional Dwell Time On the Stimulus", x = "Conditions", y = "Proportional Dwell Time") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "fixations_dwell_time_on_test.png"), width=1800, height=2400, units="px")

# Proportional Dwell Time Off The Stimulus
fixations.test.dwell <- fixations.test.valid %>% 
  filter(blok) %>% 
  mutate(dwell.time.total = sum(dur), .by=c(subject, condition)) %>% 
  filter(!ROI) %>%
  mutate(dwell.time.roi = sum(dur), .by=c(subject, condition)) %>% 
  # mutate(dwell.time.prop = dwell.time.roi/dwell.time.total)
  summarise(dwell.time.prop = mean(dwell.time.roi/dwell.time.total), .by=c(subject, condition)) 

fixations.test.dwell.summary <- fixations.test.dwell %>% 
  summarise(Mean = mean(dwell.time.prop), SD = sd(dwell.time.prop), .by=condition)

ggplot(fixations.test.dwell.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=fixations.test.dwell, aes(y = dwell.time.prop), size = 2, shape = 21, color = "black", alpha=0.5, position=position_jitter(width=0.05)) +
  labs(title = "Proportional Dwell Time Off the Stimulus", x = "Conditions", y = "Proportional Dwell Time") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "fixations_dwell_time_off_test.png"), width=1800, height=2400, units="px")

# Percentage of fixation going towards the stimuli
fixations.test.prop <- fixations.test.valid %>% 
  filter(blok) %>%
  summarise(absolute_frequency_ROI = sum(ROI), relative_frequency_ROI = mean(ROI), .by=c(subject, condition)) 

fixations.test.roi.summary <- fixations.test.prop %>% 
  summarise(Mean = mean(relative_frequency_ROI), SD = sd(relative_frequency_ROI), .by=condition)

ggplot(fixations.test.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=fixations.test.prop, aes(y = relative_frequency_ROI), size = 2, shape = 21, color = "black", alpha=0.5, position=position_jitter(width=0.05)) +
  labs(title = "Proportion of Fixations on the Stimulus", x = "Conditions", y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "fixations_proportion_on_test.png"), width=1800, height=2400, units="px")


### SACCADES
saccades.test.valid <- saccades %>% 
  left_join(baseline.test.trials %>% select(subject, trial, x_divergence, y_divergence, blok), by=c("subject", "trial"))

saccades.test.valid = saccades.test.valid %>%
  mutate(end.testtrial = end.testtrial - picOnset, # realign such that 0 = picture onset
         start_time = start_time - picOnset, end_time = end_time - picOnset, # realign such that 0 = picture onset
         start_time = ifelse(start_time < 0, 0, start_time), # discard fraction of fixation before stimulus
         end_time = ifelse(end_time > end.testtrial, end.testtrial, end_time), # discard fraction of fixation after end
         picOnset = picOnset - picOnset,
         dur = end_time - start_time) %>% filter(dur > 0)

saccades.test.valid <- saccades.test.valid %>% 
  left_join(rois, by=c("subject", "trial"))

saccades.test.valid <- saccades.test.valid %>% 
  mutate(start_x_corr = start_x - x_divergence, end_x_corr = end_x - x_divergence,
         start_y_corr = start_y - y_divergence, end_y_corr = end_y - y_divergence,  
         ROI = checkRoi(end_x_corr, end_y_corr, roi.xleft, roi.xright, roi.ybottom, roi.ytop))

saccades.test.analysis = saccades.test.valid %>% 
  mutate(start_x_corr_cm = pixToCm(start_x_corr, screen.width, screen.width.cm),
         end_x_corr_cm = pixToCm(end_x_corr, screen.width, screen.width.cm),
         start_y_corr_cm = pixToCm(start_y_corr, screen.height, screen.height.cm),
         end_y_corr_cm = pixToCm(end_y_corr, screen.height, screen.height.cm)) %>% 
  mutate(angle = visangle(start_x_corr_cm, start_y_corr_cm, end_x_corr_cm, end_y_corr_cm, distance)) 

# Latency to first saccade going towards the quadrant of stimuli
saccades.test.lat.roi <- saccades.test.analysis %>% 
  group_by(subject, trial) %>%
  filter(!ROI) %>%
  distinct(ROI, .keep_all = TRUE) %>% 
  ungroup %>%
  filter(blok) %>% 
  filter(!contains_blink) %>%
  filter(!ROI & angle >= 1) %>% 
  summarise(latency = start_time, .by=c(subject, condition)) 

saccades.test.lat.roi.summary <- saccades.test.lat.roi %>% 
  summarise(Mean = mean(latency), SD = sd(latency), .by=condition)

ggplot(saccades.test.lat.roi.summary, aes(x = condition, y = Mean, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.7),
    width = 0.25) +
  geom_point(data=saccades.test.lat.roi, aes(y = latency), size = 2, shape = 21, color = "black", alpha=0.5, position=position_jitter(width=0.05)) +
  labs(title = "Latency to First Saccade Away from the Stimulus", x = "Conditions", y = "Latency [ms]") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_d() + 
  scale_color_viridis_d()

ggsave(file.path(path, "plots", "saccades_latency_off_test.png"), width=1800, height=2400, units="px")


###############################################################################
# Messages: Jitter Check
###############################################################################
messages_image <- messages %>% 
  filter(grepl("ImageOnset", event)) %>% 
  filter(!grepl("Test", event))
with(messages_image, hist(time, main=paste0("Image Onset")))

###############################################################################
# Ratings
###############################################################################
for (p in c("Baseline", "Acquisition", "Test")) {
  ratings.phase <- ratings %>% 
    filter(phase == tolower(p))
  
  ratings.phase.summary <- ratings.phase %>% 
    summarise(Mean = mean(rating), SD = sd(rating), .by=condition)
  
  ggplot(ratings.phase.summary, aes(x = condition, y = Mean, fill = condition)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_errorbar(
      aes(ymin = Mean - SD, ymax = Mean + SD),
      position = position_dodge(width = 0.7),
      width = 0.25) +
    geom_line(data=ratings.phase, aes(y = rating, group = subject), alpha=0.2) +
    geom_point(data=ratings.phase, aes(y = rating), size = 2, shape = 21, color = "black", alpha=0.5) + #, position=position_jitter(width=0.05)) +
    labs(title = paste("Rating in ", p, "Phase", sep=" "), x = "Conditions", y = "Rating") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_viridis_d() + 
    scale_color_viridis_d()
  
  ggsave(file.path(path, "plots", paste("ratings_", tolower(p), ".png")), width=1800, height=2400, units="px")
}
