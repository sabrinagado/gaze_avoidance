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

hr_bins = T #calculate HR in 1s bins after CS onset
baselineWindow = c(-.5, 0) #correct for Baseline in this time window
step_plotting = 0.1
scaling.window = c(seq(-.5, 12, by=step_plotting)) # Scoring bins in seconds (real time scaling; may be non-integer)
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

trigger_mat <- read.csv2("../Physio/Trigger/conditions.csv") %>%
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
vpn.ecg.rpeaks = list.files("../Physio/Peak_Export/", pattern="*_rpeaks.csv", full.names=TRUE)
vpn.ecg.hr = list.files("../Physio/HR/", pattern="*.txt", full.names=TRUE)

#ratings.all = read_rds("ratings.rds" %>% paste0(path.rds, .))
#rawfiles = vpn.ecg.rpeaks %>% gsub("rpeaks/", "", .) %>% gsub(path.rpeaks.postfix, "", .) %>% paste0(".txt")

hr.list = list() #vector("list", length(vpn.ecg.rpeaks))

for (vpi in seq(vpn.ecg.hr)) {
  # vpi = 1
  
  vp = vpn.ecg.rpeaks[vpi]
  #if (vp %>% pathToCode() %>% codeToNum() %in% exclusions.hr) next
  print(vp %>% pathToCode())
  
  code = vp %>% pathToCode() %>% pathToCode(file.ext = "_") %>% substr(1, 6)  
  #ratingfile = ratings.all %>% filter(subject == code %>% codeToNum())
  
  #load r peaks
  allrpeak = read.csv2(vp)[, 1] #only first column
  hr = 60/diff(allrpeak) #convert to bpm
  
  #load triggers
  trigger <-  read.table(paste0("../Physio/HR/",code,".txt")) %>% .$V3 #lade triggerspalte von HR file mit 0er und marker
  trigger <- ifelse(trigger > max(markers), 0, trigger)

  # if (exclusions.phys.trials[[code]] %>% is.null() == F) { #exclude trials manually (cp. triggerCheck)
  #   triggers.only = trigger[trigger != 0]
  #   toExclude = exclusions.phys.trials[[code]]
  #   triggers.only[toExclude] = 0
  #   trigger[trigger != 0] = triggers.only
  # }
  triggers.n = sum(trigger!=0)
  trials.n = triggers.n
  
  # exclude first wrong triggers for VP6
  # if(triggers.n > 180) {
  #    trigger[head(which(trigger > 0),16)] <- 0
  #    triggers.n = sum(trigger!=0)
  #  }
  
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
  allhr = cbind(allhrbl, allhr)
  deltaval = allhr - matrix(allhr[, 1], nrow=nrow(allhr), ncol=ncol(allhr))
  deltavalbins = allhrbins - matrix(allhrbins[, 1], nrow=nrow(allhrbins), ncol=ncol(allhrbins))
  
  out = data.frame(trial = 1:trials.n, 
                   condition = conditions,
                   hrbl = allhr[, 1], 
                   hr.bin = deltavalbins[, 2:ncol(deltavalbins)],
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

saveRDS(heart.wide,"HR.RData")

heart = heart.wide %>% gather(key="time", value="HRchange", matches("hr\\.\\d+")) %>% tibble() %>% 
  mutate(time = time %>% gsub("hr.", "", .) %>% as.integer() %>% {. * step_plotting + min(baselineWindow)} %>% as.factor(),
         condition = as.factor(condition)) %>% 
  select(-contains("hr.bin."))

heart = heart %>% 
  left_join(trigger_mat %>% select(subject, trial, outcome) %>% mutate(subject = as.integer(substr(subject, 5, 6))), by=c("subject", "trial")) %>% 
  mutate(outcome = ifelse(is.na(outcome), "no outcome", outcome))

# heart.data <- heart.wide %>% mutate(ID = paste0("gca",str_pad(subject,2,pad="0"))) %>% select(ID, trial, hrbl:hr.20)


# Plots -------------------------------------------------------------
# heart.wide <- readRDS("HR.RData")

# Acquisition
# Long format for statistical testing
hr_df_long_acq <- heart.wide %>%
  # filter(outcome == "no outcome") %>%
  #filter(ID %in% responder) %>%
  pivot_longer(hr.bin.1:hr.bin.16, names_to = "timebin", values_to ="HR", names_prefix="hr.bin.") %>%
  mutate(timebin = as.numeric(timebin)) %>%
  filter(condition %in% c(2,3,4,5)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>%
  mutate(ID = subject) %>% 
  select(ID, trial, condition, outcome, timebin, HR) %>%
  mutate(time = (timebin - 1) * 0.5) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

main_effect_threat = list()
main_effect_social = list()
interaction_effect = list()
alpha = .05 / length(unique(hr_df_long_acq$time))
for (timepoint in unique(hr_df_long_acq$time)) {
  # timepoint = 3
  data = hr_df_long_acq %>% 
    filter(time == timepoint) %>% 
    mutate(ID = as.factor(ID), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
    group_by(ID, condition_social, condition_threat)
  model = lmer(HR ~ condition_social + condition_threat + condition_social:condition_threat + (1|ID), data)
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
heart %>% 
  filter(condition %in% c(2,3,4,5)) %>%
  # filter(outcome != "no outcome") %>%
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  group_by(condition, time) %>% 
  summarise(HR.se = sd(HRchange, na.rm=T)/sqrt(n()), HR.mean = mean(HRchange, na.rm=T)) %>% 
  mutate(time = time %>% as.character() %>% as.numeric()) %>%
  {ggplot(., aes(x=time, y=HR.mean)) +
      geom_vline(xintercept=0, colour="black",linetype="solid") + #zero
      geom_line(aes(colour=condition)) +
      geom_ribbon(aes(ymin=HR.mean-HR.se, ymax=HR.mean+HR.se, colour=condition, fill=condition), color = NA, alpha=.2) +
      # geom_segment(data = main_effect_social, aes(x=times_start, xend = times_end, y=-2.7, yend=-2.7, size="Main Effect Social"), colour = "#e874ff", linewidth = 1, inherit.aes=FALSE) +
      geom_segment(data = main_effect_threat, aes(x=times_start, xend = times_end, y=-3, yend=-3, size="Main Effect Threat"), colour = "#ff8383", linewidth = 1, inherit.aes=FALSE) +
      # geom_segment(data = interaction_effect, aes(x=times_start, xend = times_end, y=-3.3, yend=-3.3, size="Social x Threat Interaction "), colour = "#ffdd74", linewidth = 1, inherit.aes=FALSE) +
      scale_x_continuous("Time [s]",limits=c(-0.5, 8), minor_breaks=c(0,1,2,3,4,5,6,7,8), breaks=c(0, 2, 4, 6, 8)) +
      scale_y_continuous("Heart Rate",limits=c(-3, 5.5)) +
      scale_color_viridis_d(aesthetics = c("colour", "fill")) +
      theme_bw() +
      scale_size_manual("effects", values=rep(1,3), guide=guide_legend(override.aes = list(colour=c("#ff8383")))) # "#e874ff", "#ff8383", "#ffdd74"
  }
ggsave(paste0("../plots/HR/cs_acq.png"), type="cairo-png", width=2500/400, height=1080/300, dpi=300)


# Test
# Long format for statistical testing
hr_df_long_test <- heart.wide %>%
  # filter(outcome == "no outcome") %>%
  #filter(ID %in% responder) %>%
  pivot_longer(hr.bin.1:hr.bin.19, names_to = "timebin", values_to ="HR", names_prefix="hr.bin.") %>%
  mutate(timebin = as.numeric(timebin)) %>%
  filter(condition %in% c(6,7,8,9)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>%
  mutate(ID = subject) %>% 
  select(ID, trial, condition, outcome, timebin, HR) %>%
  mutate(time = (timebin - 1) * 0.5) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

main_effect_threat = list()
main_effect_social = list()
interaction_effect = list()
alpha = .05 / length(unique(hr_df_long_test$time))
for (timepoint in unique(hr_df_long_test$time)) {
  # timepoint = 0
  data = hr_df_long_test %>% 
    filter(time == timepoint) %>% 
    mutate(ID = as.factor(ID), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
    group_by(ID, condition_social, condition_threat)
  
  model = lmer(HR ~ condition_social + condition_threat + condition_social:condition_threat + (1|ID), data)
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

heart %>% 
  filter(condition %in% c(6,7,8,9)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>%
  group_by(condition, time) %>% 
  summarise(HR.se = sd(HRchange, na.rm=T)/sqrt(n()), HR.mean = mean(HRchange, na.rm=T)) %>% 
  mutate(time = time %>% as.character() %>% as.numeric()) %>% 
  {ggplot(., aes(x=time, y=HR.mean)) +
      geom_vline(xintercept=0, color="black",linetype="solid") + #zero = picture onset
      geom_vline(xintercept=10, color="black",linetype="solid") + #picture offset
      geom_line(aes(colour=condition)) +
      geom_ribbon(aes(ymin=HR.mean-HR.se, ymax=HR.mean+HR.se, colour=condition, fill=condition), color = NA, alpha=.2) +
      geom_segment(data = main_effect_social, aes(x=times_start, xend = times_end, y=-2.7, yend=-2.7, size="Main Effect Social"), colour = "#e874ff", linewidth = 1, inherit.aes=FALSE) +
      geom_segment(data = main_effect_threat, aes(x=times_start, xend = times_end, y=-3, yend=-3, size="Main Effect Threat"), colour = "#ff8383", linewidth = 1, inherit.aes=FALSE) +
      geom_segment(data = interaction_effect, aes(x=times_start, xend = times_end, y=-3.3, yend=-3.3, size="Social x Threat Interaction "), colour = "#ffdd74", linewidth = 1, inherit.aes=FALSE) +
      scale_x_continuous("Time [s]",limits=c(-0.5, 11), minor_breaks=c(0,1,2,3,4,5,6,7,8,9,10), breaks=c(0, 2, 4, 6, 8, 10)) +
      scale_y_continuous("Heart Rate") + #, breaks=c(-80,-40, 0, 40), minor_breaks=c(-80, -60, -40, -20, 0, 20, 40, 60)) +
      scale_color_viridis_d(aesthetics = c("colour", "fill")) +
      theme_bw() +
      scale_size_manual("effects", values=rep(1,4), guide=guide_legend(override.aes = list(colour=c("#e874ff", "#ff8383", "#ffdd74"))))
  }
ggsave(paste0("../Plots/HR/cs_test.png"), type="cairo-png", width=2500/400, height=1080/300, dpi=300)
