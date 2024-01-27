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
bin.window = c(seq(-.5, 12, by=bin_width)) # Scoring bins in seconds (real time scaling; may be non-integer)

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

#ratings.all = read_rds("ratings.rds" %>% paste0(path.rds, .))
#rawfiles = vpn.ecg.rpeaks %>% gsub("rpeaks/", "", .) %>% gsub(path.rpeaks.postfix, "", .) %>% paste0(".txt")

hr.list = list() #vector("list", length(vpn.ecg.rpeaks))

for (vpi in seq(vpn.ecg.rpeaks)) {
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
  select(subject, trial, condition, everything())
rm(hr.list); row.names(heart.wide) = NULL

heart = heart.wide %>% gather(key="time", value="HRchange", matches("hr\\.\\d+")) %>% tibble() %>% 
  mutate(time = time %>% gsub("hr.", "", .) %>% as.integer() %>% {. * step_plotting + min(baselineWindow)} %>% as.factor(),
         condition = as.factor(condition)) %>% 
  select(-contains("hr.bin."))

heart = heart %>% 
  left_join(trigger_mat %>% select(subject, trial, outcome) %>% mutate(subject = as.integer(substr(subject, 5, 6))), by=c("subject", "trial")) %>% 
  mutate(outcome = ifelse(is.na(outcome), "no outcome", outcome))

heart.data <- heart.wide %>% mutate(ID = paste0("gca",str_pad(subject,2,pad="0"))) %>% select(ID, trial, hrbl:hr.20)

#write_rds(heart, "Heart_df.rds")

# Plots -------------------------------------------------------------
#heart = read_rds("heart_df.rds" )

#plot hr change over trial time
# heart.ga.gen.timeOnly = heart %>% filter(phase == "Gen") %>% group_by(time) %>% 
#   summarise(HRchange.se = se(HRchange, na.rm=T), HRchange = mean(HRchange, na.rm=T)) %>% 
#   mutate(time = time %>% as.character() %>% as.numeric()) %>% 
#   bind_rows(data.frame(time=0, HRchange=0, HRchange.se=0)) #add origin
# print(heart.trialtimeOnly.plot <- heart.ga.gen.timeOnly %>% ggplot(aes(x=time, y=HRchange)) + 
#         geom_hline(yintercept = 0, linetype="dashed") +
#         geom_ribbon(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96), color=NA, alpha=.1) +
#         #geom_errorbar(aes(ymin=HRchange-HRchange.se*1.96, ymax=HRchange+HRchange.se*1.96)) +
#         geom_point(size=3) + geom_line() + 
#         ylab("Heart Rate Change (bpm)") + xlab("Trial Time (sec)") + myGgTheme)

#plot hr change over trial time
heart %>% 
  filter(condition %in% c(2,3,4,5)) %>%
  # filter(outcome != "no outcome") %>%
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  group_by(condition, time) %>% 
  summarise(HR.se = sd(HRchange, na.rm=T)/sqrt(n()), HR.mean = mean(HRchange, na.rm=T)) %>% 
  mutate(time = time %>% as.character() %>% as.numeric()) %>%
  {ggplot(., aes(x=time, y=HR.mean, color=condition, group=condition, fill=condition)) +
      geom_vline(xintercept=0, color="black",linetype="solid") + #zero
      geom_line() +
      geom_ribbon(aes(ymin=HR.mean-HR.se, ymax=HR.mean+HR.se), color = NA, alpha=.2) +
      scale_x_continuous("Time [s]",limits=c(-0.5, 8), minor_breaks=c(0,1,2,3,4,5,6,7,8), breaks=c(0, 2, 4, 6, 8)) +
      scale_y_continuous("Heart Rate",limits=c(-3, 5.5)) +
      scale_color_viridis_d(aesthetics = c("colour", "fill")) +
      theme_bw()
  }
ggsave(paste0("../plots/HR/cs_acq.png"), type="cairo-png", width=2500/400, height=1080/300, dpi=300)

#Test phase
heart %>% 
  filter(condition %in% c(6,7,8,9)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>%
  group_by(condition, time) %>% 
  summarise(HR.se = sd(HRchange, na.rm=T)/sqrt(n()), HR.mean = mean(HRchange, na.rm=T)) %>% 
  mutate(time = time %>% as.character() %>% as.numeric()) %>% 
  {ggplot(., aes(x=time, y=HR.mean, color=condition, group=condition, fill=condition)) +
      geom_vline(xintercept=0, color="black",linetype="solid") + #zero
      geom_line() +
      geom_ribbon(aes(ymin=HR.mean-HR.se, ymax=HR.mean+HR.se), color = NA, alpha=.2) +
      scale_x_continuous("Time [s]",limits=c(-0.5, 11), minor_breaks=c(0,1,2,3,4,5,6,7,8,9,10), breaks=c(0, 2, 4, 6, 8, 10)) +
      scale_y_continuous("Heart Rate") + #, breaks=c(-80,-40, 0, 40), minor_breaks=c(-80, -60, -40, -20, 0, 20, 40, 60)) +
      scale_color_viridis_d(aesthetics = c("colour", "fill")) +
      theme_bw()
  }
ggsave(paste0("../Plots/HR/cs_test.png"), type="cairo-png", width=2500/400, height=1080/300, dpi=300)



# heart %>% 
#   filter(trial_condition <= 30) %>% 
#   group_by(condition, time) %>% 
#   summarise(HRchange.se = sd(HRchange, na.rm=T)/sqrt(n()), HRchange = mean(HRchange, na.rm=T)) %>% 
#   mutate(time = time %>% as.character() %>% as.numeric()) %>% 
#   bind_rows(data.frame(condition=unique(heart$condition), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
#   ggplot(aes(x=time-0.5,y=HRchange, fill = condition)) +
#   #geom_rect(aes(xmin=4, xmax=6, ymin=-2, ymax=1), fill = "deepskyblue1", alpha=0.03) +
#   geom_errorbar(aes(ymin=HRchange-HRchange.se,ymax=HRchange+HRchange.se), width=0.2)+
#   geom_line()+
#   geom_point(size=3, shape=21)+
#   scale_fill_manual(values=c("red","orange","darkgreen")) +
#   scale_x_continuous(name ="Time [s]", breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12)) + 
#   scale_y_continuous(name="Heart rate change [bpm]", expand= c(0,0))+
#   ggtitle("1st Half") +
#   theme_classic() +
#   theme(legend.position=c(0.1,0.18),
#         legend.title = element_blank(),
#         #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#   )
# 
# heart %>% 
#   filter(trial_condition > 30) %>% 
#   group_by(condition, time) %>% 
#   summarise(HRchange.se = sd(HRchange, na.rm=T)/sqrt(n()), HRchange = mean(HRchange, na.rm=T)) %>% 
#   mutate(time = time %>% as.character() %>% as.numeric()) %>% 
#   bind_rows(data.frame(condition=unique(heart$condition), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
#   ggplot(aes(x=time-0.5,y=HRchange, fill = condition)) +
#   #geom_rect(aes(xmin=4, xmax=6, ymin=-2, ymax=1), fill = "deepskyblue1", alpha=0.03) +
#   geom_errorbar(aes(ymin=HRchange-HRchange.se,ymax=HRchange+HRchange.se), width=0.2)+
#   geom_line()+
#   geom_point(size=3, shape=21)+
#   scale_fill_manual(values=c("red","orange","darkgreen")) +
#   scale_x_continuous(name ="Time [s]", breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12)) + 
#   scale_y_continuous(name="Heart rate change [bpm]", expand= c(0,0))+
#   ggtitle("2nd Half") +
#   theme_classic() +
#   theme(legend.position=c(0.1,0.18),
#         legend.title = element_blank(),
#         #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#   )
# 
# 
# # Interference statistics
# 
# anova <- heart %>%  
#   filter(!(time %in% c(10.5,11,11.5,12,12.5,13,13.5,14,14.5,15))) %>%
#   group_by(subject,time,condition) %>% 
#   summarize(HRchange = mean(HRchange)) %>% ungroup() %>%
#   mutate(subject = as.factor(subject), condition = as.factor(condition), time = as.factor(time)) %>%
#   ez::ezANOVA(dv=.(HRchange), wid=.(subject), 
#               within=.(condition, time), 
#               #between=.(pairs),
#               detailed=T, type=2)
# 
# 
# anova %>% anova_apa()
# anova$ANOVA[2,] %>% partial_eta_squared_ci()
# anova$ANOVA[3,] %>% partial_eta_squared_ci()
# anova$ANOVA[4,] %>% partial_eta_squared_ci()
# 
# 
# heart %>%  
#   filter(!(time %in% c(10.5,11,11.5,12,12.5,13,13.5,14,14.5,15))) %>%
#   group_by(subject,time,condition) %>% 
#   summarize(HRchange = mean(HRchange)) %>%
#   mutate(subject = as.factor(subject), condition = as.factor(condition), time = as.factor(time)) %>%
#   aov_ez(dv="HRchange", id="subject", . ,within = c("condition","time")) %>% 
#   emmeans::emmeans(~ condition|time) %>%
#   pairs() %>% summary(adjust = "FDR")
# 
# 
# 
# 
# 
# heart %>% 
#   filter(!(time %in% c(11,11.5,12,12.5,13,13.5,14,14.5,15))) %>%
#   group_by(condition, time) %>% 
#   summarise(HRchange.se = sd(HRchange, na.rm=T)/sqrt(n()), HRchange = mean(HRchange, na.rm=T)) %>% 
#   mutate(time = time %>% as.character() %>% as.numeric()) %>% 
#   #bind_rows(data.frame(condition=unique(heart$condition), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
#   ggplot(aes(x=time-0.5,y=HRchange, group = condition, color = condition)) +
#   #geom_rect(aes(xmin=4, xmax=6, ymin=-2, ymax=1), fill = "deepskyblue1", alpha=0.03) +
#   #geom_errorbar(aes(ymin=HRchange-HRchange.se,ymax=HRchange+HRchange.se), width=0.2)+
#   geom_path()+
#   #geom_point(size=3, shape=21)+
#   scale_color_manual(values=c("red","orange","darkgreen"), label = c("Threat","Flight","Safety")) +
#   scale_x_continuous(name ="Time [s]", limits = c(0,10), breaks = c(0,1,2,3,4,5,6,7,8,9,10)) + 
#   scale_y_continuous(name="Heart rate change [bpm]")+
#   theme_classic() +
#   theme(legend.position=c(0.15,0.18),
#         legend.title = element_blank(),
#         #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#   )
# 
# ggsave(paste0("../Plots/HR/acq_cs_hr_bl_lines.png"), type="cairo-png", width=1920/400, height=1080/300, dpi=300)
# 
# # with US 
# heart %>% 
#   #filter(!(time %in% c(11,11.5,12,12.5,13,13.5,14,14.5,15))) %>%
#   group_by(condition, time) %>% 
#   summarise(HRchange.se = sd(HRchange, na.rm=T)/sqrt(n()), HRchange = mean(HRchange, na.rm=T)) %>% 
#   mutate(time = time %>% as.character() %>% as.numeric()) %>% 
#   #bind_rows(data.frame(condition=unique(heart$condition), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
#   ggplot(aes(x=time-0.5,y=HRchange, group = condition, color = condition)) +
#   #geom_rect(aes(xmin=4, xmax=6, ymin=-2, ymax=1), fill = "deepskyblue1", alpha=0.03) +
#   geom_ribbon(aes(ymin=HRchange-HRchange.se,ymax=HRchange+HRchange.se, fill = condition), color = "transparent", alpha=0.2, show_guide = FALSE)+
#   geom_path()+
#   #geom_point(size=3, shape=21)+
#   scale_color_manual(values=c("red","orange","darkgreen"), label = c("Threat","Flight","Safety")) +
#   scale_fill_manual(values=c("red","orange","darkgreen"), label = c("Threat","Flight","Safety")) +
#   scale_x_continuous(name ="Time [s]") + 
#   scale_y_continuous(name="Heart rate change [bpm]")+
#   theme_classic() +
#   theme(legend.position=c(0.15,0.18),
#         legend.title = element_blank(),
#         #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#   )
# 
# ggsave(paste0("../Plots/HR/acq_cs_hr_bl_lines.png"), type="cairo-png", width=1920/400, height=1080/300, dpi=300)
# 
# 
# 
# 
# 
# 
# 
# hr_out <- heart %>% group_by(subject, condition, time) %>%
#   summarise(hr = mean(HRchange)) %>%
#   mutate(time = as.numeric(time),
#          ID = paste0("tfo",str_pad(subject,2,pad="0")),
#          condition = factor(condition, levels = c(1,2,3), labels = c("threat", "flight", "safety")),
#          names = paste0(condition,"_",time)) %>% select(-subject) %>%
#   pivot_wider(id_cols = ID, names_from = names, values_from = hr) 
# hr_out <- hr_out[order(hr_out$ID),]
# write.csv2(hr_out,"HR_Daten.csv", row.names = F)
# 
# 
# 
# heart_wue <- heart %>% group_by(subject, condition, time) %>%
#   summarise(hr = mean(HRchange)) %>% ungroup() %>%
#   mutate(bin = as.numeric(time),
#          ID = subject,
#          cue = factor(condition, levels = c(1,2,3), labels = c("threat", "flight", "safety"))) %>% 
#   select(ID, cue, bin, hr, -subject) %>% mutate(site = "wuerzburg")
# 
# 
# # Gainesville Data
# 
# library(R.matlab)
# 
# hr_gainesville <- readMat("../Florida Sample/HR/allsubjects.mat")
# 
# hr_gainesville1 <- data.frame(hr_gainesville[[1]]) %>% 
#   mutate(code = c(1:49,51,52),
#          cue = "threat") %>%
#   pivot_longer(X1:X22, names_to = "bin", values_to = "hr") %>%
#   mutate(bin = as.numeric(str_remove(bin,"X"))-1)
# 
# hr_gainesville2 <- data.frame(hr_gainesville[[2]]) %>% 
#   mutate(code = c(1:49,51,52),
#          cue = "flight") %>%
#   pivot_longer(X1:X22, names_to = "bin", values_to = "hr") %>%
#   mutate(bin = as.numeric(str_remove(bin,"X"))-1)
# 
# hr_gainesville3 <- data.frame(hr_gainesville[[3]]) %>% 
#   mutate(code = c(1:49,51,52),
#          cue = "safety") %>%
#   pivot_longer(X1:X22, names_to = "bin", values_to = "hr") %>%
#   mutate(bin = as.numeric(str_remove(bin,"X"))-1)
# 
# hr_ga <- rbind(hr_gainesville1,hr_gainesville2) %>% rbind(.,hr_gainesville3) %>%
#   mutate(code = code+100) 
# 
# 
# 
# 
# 
# hr_ga %>% 
#   mutate(cue = factor(cue, levels = c("threat","flight","safety"), labels = c("threat", "flight", "safety"))) %>%
#   group_by(cue, bin) %>% 
#   summarise(hr.se = sd(hr, na.rm=T)/sqrt(n()), hr = mean(hr, na.rm=T)) %>% 
#   #bind_rows(data.frame(cue=unique(heart$cue), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
#   ggplot(aes(x=bin/2,y=hr, group = cue, color = cue)) +
#   #geom_rect(aes(xmin=4, xmax=6, ymin=-2, ymax=1), fill = "deepskyblue1", alpha=0.03) +
#   geom_ribbon(aes(ymin=hr-hr.se,ymax=hr+hr.se, fill = cue), color = "transparent", alpha=0.1, show_guide = FALSE)+
#   geom_path()+
#   #geom_point(size=3, shape=21)+
#   scale_color_manual(values=c("red","orange","darkgreen"), label = c("Threat","Flight","Safety")) +
#   scale_fill_manual(values=c("red","orange","darkgreen"), label = c("Threat","Flight","Safety")) +
#   #scale_x_continuous(name ="Time [s]", breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)) + 
#   scale_y_continuous(name="Heart rate change [bpm]")+
#   theme_classic() +
#   theme(legend.position=c(0.15,0.88),
#         legend.title = element_blank(),
#         #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#   )
# 
# 
# anova <- hr_ga %>%  
#   filter(!(bin %in% c(0, 21))) %>%
#   filter((bin %in% c(15,16,17,18))) %>%
#   group_by(code,bin,cue) %>% 
#   summarize(hr = mean(hr)) %>% ungroup() %>%
#   mutate(subject = as.factor(subject), cue = as.factor(cue), bin = as.factor(bin)) %>%
#   ez::ezANOVA(dv=.(hr), wid=.(code), 
#               within=.(cue, bin), 
#               #between=.(pairs),
#               detailed=T, type=2)
# 
# 
# anova %>% anova_apa()
# anova$ANOVA[2,] %>% partial_eta_squared_ci()
# anova$ANOVA[3,] %>% partial_eta_squared_ci()
# anova$ANOVA[4,] %>% partial_eta_squared_ci()
# 
# 
# # All together
# 
# hr_all <- rbind(heart_wue, hr_ga %>% filter(time!=0) %>% mutate(site = "gainesville"))
# 
# 
# 
# hr_all %>% 
#   filter(time < 21) %>% 
#   group_by(cue, bin) %>% 
#   summarise(hr.se = sd(hr, na.rm=T)/sqrt(n()), hr = mean(hr, na.rm=T)) %>% 
#   #bind_rows(data.frame(condition=unique(heart$condition), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
#   ggplot(aes(x=bin/2,y=hr, group = cue, color = cue)) +
#   #geom_rect(aes(xmin=4, xmax=6, ymin=-2, ymax=1), fill = "deepskyblue1", alpha=0.03) +
#   geom_errorbar(aes(ymin=hr-hr.se,ymax=hr+hr.se), width=0.2)+
#   geom_path()+
#   #geom_point(size=3, shape=21)+
#   scale_color_manual(values=c("darkred","orange","blue")) +
#   #scale_x_continuous(name ="Time [s]", breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)) + 
#   scale_y_continuous(name="Heart rate change [bpm]")+
#   theme_classic() +
#   theme(legend.position=c(0.15,0.18),
#         legend.title = element_blank(),
#         #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#   )
# 
# 
# anova <- hr_all %>%  
#   filter(bin != 0) %>%
#   filter(bin < 21) %>%
#   group_by(ID,bin,cue) %>% 
#   summarize(hr = mean(hr)) %>% ungroup() %>%
#   mutate(ID = as.factor(ID), cue = as.factor(cue), bin = as.factor(bin)) %>%
#   ez::ezANOVA(dv=.(hr), wid=.(ID), 
#               within=.(cue, bin), 
#               #between=.(pairs),
#               detailed=T, type=2)
# 
# 
# anova %>% anova_apa()
# anova$ANOVA[2,] %>% partial_eta_squared_ci()
# anova$ANOVA[3,] %>% partial_eta_squared_ci()
# anova$ANOVA[4,] %>% partial_eta_squared_ci()
