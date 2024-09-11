library(tidyverse)
library(afex)
library(apa)
library(lme4)
library(lmerTest)

#source("0 General.R" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/C10 - Gamer Fear/3 Diagnostic Generalization/3 Analysis/main Anx/", .))

#exclusions.hr = c(
#  40:43, #FFP2 mask (interacts with breathing more heavily)
#  54, #temporary exclusion until workaround for missing ratings :D
#  95 #too many ectopic beats
#) %>% c(exclusions) %>% unique() %>% sort()

markers = c(2,3,4,5,6,7,8,9) # according to markers after importing social vs nonsocial and cs information
sample.rate = 100
minhr =  45 #minimum plausible heart rate
maxhr = 120 #maximum plausible heart rate

hr_bins = T #calculate HR in 0.5s bins after CS onset
baselineWindow = c(-.5, 0) #correct for Baseline in this time window
step_plotting = 0.5
scaling.window = c(seq(-4, 12, by=step_plotting)) # Scoring bins in seconds (real time scaling; may be non-integer)
bin_width = 0.5
bin.window = c(seq(0, 10, by=bin_width)) # Scoring bins in seconds (real time scaling; may be non-integer)

{ # Functions ---------------------------------------------------------------
  scaleHR = function(hr_t, hr, st, en) {
    wsum = 0 #weighted sum (will be divided by valid ranges below)
    naRanges = 0
    i = which(hr_t > st)[1] #first r peak after marker start (might be NA)
    while (!is.na(i) && i <= length(hr_t) && hr_t[i] < en) { #build up weighted sum until last r peak in interval (exclusively! see below)
      #determine range
      lower = ifelse(i > 1 && hr_t[i-1] > st, hr_t[i-1], st) #find lower end of scoring interval (might be st if previous beat does not exist or out or scoring range)
      range = hr_t[i] - lower
      
      #determine value
      if (i > 1 && !is.na(hr[i])) value = hr[i]
      else {
        value = 0
        naRanges = naRanges + range #lower total range by invalid ranges
      }
      
      wsum = wsum + range * value
      i = i + 1
    } #wsum and naRanges are computed except for last beat
    
    #last beat
    if (!is.na(i) && i > 1 && i <= length(hr_t)) {
      range = en - max(c(hr_t[i-1], st)) #it might happen that there is only one beat before st and one after en => use hr[i] (implicitly done by setting range = en - st via max(c(...)))
      
      value = hr[i]
      if (is.na(value)) {
        value = 0
        naRanges = naRanges + range #lower total range by invalid ranges
      } else value = hr[i]
      
      wsum = wsum + range * value
    }
    
    result = wsum / (en - st - naRanges)
    if (is.na(result) || result ==  0) result = NA
    return(result)
  }
  
  pathToCode = function(path, path.sep="/", file.ext="\\.") {
    first = path %>% gregexpr(path.sep, .) %>% lapply(max) %>% unlist() %>% {. + 1}
    last = path %>% gregexpr(file.ext, .) %>% lapply(max) %>% unlist() %>% {. - 1}
    return(path %>% substring(first, last))
  }
}

path.physio = file.path("Study 1", "Physio")

trigger_mat <- read.csv2(file.path(path.physio, "Trigger", "conditions.csv")) %>%
  mutate(subject = sprintf("gca_%02d", subject),
         trigger = ifelse(phase == "acquisition" & condition == "CSneg, social",2,
                          ifelse(phase == "acquisition" & condition == "CSpos, social",3,
                                 ifelse(phase == "acquisition" & condition == "CSneg, non-social",4,
                                        ifelse(phase == "acquisition" & condition == "CSpos, non-social",5,
                                               ifelse(phase == "test" & condition == "CSneg, social",6,
                                                      ifelse(phase == "test" & condition == "CSpos, social",7,
                                                             ifelse(phase == "test" & condition == "CSneg, non-social",8,
                                                                    ifelse(phase == "test" & condition == "CSpos, non-social",9,0)))))))))

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


# Read & Score HR --------------------------------------------------------------------
vpn.ecg.rpeaks = list.files(file.path(path.physio, "Peak_Export"), pattern="*_rpeaks.csv", full.names=TRUE)
vpn.ecg.hr = list.files(file.path(path.physio, "HR"), pattern="*.txt", full.names=TRUE)

#ratings.all = read_rds("ratings.rds" %>% paste0(path.rds, .))
#rawfiles = vpn.ecg.rpeaks %>% gsub("rpeaks/", "", .) %>% gsub(path.rpeaks.postfix, "", .) %>% paste0(".txt")

hr.list = list() #vector("list", length(vpn.ecg.rpeaks))

for (vpi in seq(vpn.ecg.hr)) {
  # vpi = 22
  
  vp = vpn.ecg.rpeaks[vpi]
  #if (vp %>% pathToCode() %>% codeToNum() %in% exclusions.hr) next
  print(vp %>% pathToCode())
  
  code = vp %>% pathToCode() %>% pathToCode(file.ext = "_") %>% substr(1, 6)  
  #ratingfile = ratings.all %>% filter(subject == code %>% codeToNum())
  
  #load r peaks
  allrpeak = read.csv2(vp)[, 1] #only first column
  hr = 60/diff(allrpeak) #convert to bpm
  
  #load triggers
  trigger <-  read.table(file.path(path.physio, "HR", paste0(code,".txt"))) %>% .$V3 #lade triggerspalte von HR file mit 0er und marker
  trigger <- ifelse(trigger > max(markers), 0, trigger)

  # if (exclusions.phys.trials[[code]] %>% is.null() == F) { #exclude trials manually (cp. triggerCheck)
  #   triggers.only = trigger[trigger != 0]
  #   toExclude = exclusions.phys.trials[[code]]
  #   triggers.only[toExclude] = 0
  #   trigger[trigger != 0] = triggers.only
  # }
  triggers.n = sum(trigger!=0)
  trials.n = triggers.n
  
  conditions = trigger[trigger!=0]
  
  timeline = seq(trigger) / sample.rate
  marker = timeline[{trigger != 0} %>% which() %>% tail(trials.n)]
  
  missingBegin = min(marker) + min(scaling.window); missingBegin = missingBegin[missingBegin <= 0] # <= 0 because first sample point is 1 / sample.rate, i.e. 0 is out of range
  missingEnd = max(timeline) - (max(marker) + max(scaling.window)); missingEnd = missingEnd[missingEnd < 0] # < 0 because if difference exactly 0, then last sample point in rage
  
  #check if data is missing at the edges of data
  if (!is_empty(missingBegin)) warning(paste(code, ": Lacking data at BEGINNING", -missingBegin, "sec"))
  if (!is_empty(missingEnd)) warning(paste(code, ": Lacking data at END", -missingEnd, "sec"))
  
  #check for plausiblilty and issue warning
  if (min(hr) < minhr || max(hr) > maxhr) warning(paste(code, ": Implausible heart rate", hr %>% min() %>% round(1), "-", hr %>% max() %>% round(1)))
  
  hr = c(NA, hr) #matching allrpeak[i] to hr[i] (heart rate values only valid if PREVIOUS r peak exists)
  
  # Real time scoring
  allhr = numeric()
  for (trial in seq(marker)) {
    # trial = 1
    mtime = marker[trial]
    hr_t = allrpeak - mtime #heart rate time (relative to marker)
    
    # Real time scaling
    hrtrial = numeric()
    for (j in seq(scaling.window)[-1]) { #for all indices except for the first (due to j-1 indexing)
      # j = seq(scaling.window)[-1][[1]]
      current = ifelse(mtime + scaling.window[j] < 0, NA, #skip marker time points that refer to negative times (i.e. out of data)
                       scaleHR(hr_t, hr, scaling.window[j-1], scaling.window[j]))
      hrtrial = c(hrtrial, current)
    }
    
    allhr = rbind(allhr, hrtrial)
  }
  
  # Bins
  allhrbins = numeric()
  for (trial in seq(marker)) {
    # trial = 1
    mtime = marker[trial]
    hr_t = allrpeak - mtime #heart rate time (relative to marker)
    
    # Real time scaling
    hrtrial = numeric()
    for (j in seq(bin.window)[-1]) { #for all indices except for the first (due to j-1 indexing)
      # j = seq(scaling.window)[-1][[1]]
      current = ifelse(mtime + bin.window[j] < 0, NA, #skip marker time points that refer to negative times (i.e. out of data)
                       scaleHR(hr_t, hr, bin.window[j-1], bin.window[j]))
      hrtrial = c(hrtrial, current)
    }
    
    allhrbins = rbind(allhrbins, hrtrial)
  }
  
  # Baseline Correction using first 0.5 seconds before trial
  allhrbl = numeric()
  for (trial in seq(marker)) {
    # trial = 1
    mtime = marker[trial]
    hr_t = allrpeak - mtime #heart rate time (relative to marker)
    
    # Real time scaling
    blhrtrial = numeric()
    for (j in seq(baselineWindow)[-1]) { #for all indices except for the first (due to j-1 indexing)
      # j = seq(scaling.window)[-1][[1]]
      current = ifelse(mtime + baselineWindow[j] < 0, NA, #skip marker time points that refer to negative times (i.e. out of data)
                       scaleHR(hr_t, hr, baselineWindow[j-1], baselineWindow[j]))
      blhrtrial = c(blhrtrial, current)
    }
    
    allhrbl = rbind(allhrbl, blhrtrial)
  }
  
  # Add baseline column to dataframes
  allhr = cbind(allhrbl, allhr)
  allhrbins = cbind(allhrbl, allhrbins)
  
  
  deltaval = allhr - matrix(allhr[, 1], nrow=nrow(allhr), ncol=ncol(allhr))
  deltavalbins = allhrbins - matrix(allhrbins[, 1], nrow=nrow(allhrbins), ncol=ncol(allhrbins))
  
  out = data.frame(trial = 1:trials.n, 
                   condition = conditions,
                   hrbl = allhr[, 1], 
                   hr.bin = deltavalbins[, 2:(ncol(deltavalbins))],
                   hr = deltaval[, 2:ncol(deltaval)]
                   )
  
  hr.list[[code]] = out
  #write.csv2(out, paste(savepath.ecg, code,"_task.csv",sep=""), row.names=FALSE, quote=FALSE)
}

#list to one giant dataframe
heart.wide = hr.list %>% bind_rows(.id="subject") %>% 
  mutate(subject = subject %>% gsub("\\D", "", .) %>% as.integer()) %>% tibble() %>% 
  # group_by(subject, condition) %>% mutate(trial_condition = 1:n()) %>% ungroup() %>% 
  select(subject, trial, condition, everything()) %>% 
  left_join(trigger_mat %>% select(subject, trial, outcome) %>% mutate(subject = as.integer(substr(subject, 5, 6))), by=c("subject", "trial")) %>% 
  mutate(outcome = ifelse(is.na(outcome), "no outcome", outcome))

rm(hr.list); row.names(heart.wide) = NULL

saveRDS(heart.wide,file.path("Study 1", "Physio", "HR.RData"))
