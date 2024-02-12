###############################################################################
# Gaze Contingent Avoidance Project
# Sabrina Gado & Yannik Stegmann
# Code adapted from Mario Reutter and Matthias Gamer

# EDA Analyse

#packages 
library(signal)
library(tidyverse)
library(stringr)
library(apa)
library(afex)
library(lme4)
library(lmerTest)
options(dplyr.summarise.inform = FALSE)


#global parameters, options and plots

ucr_plots = T
cr_plots = T

downsampling = T
sample_rate_new = 50 #for downsampling
check_ds_plots = F

lowpass = T
lowPassFreq = 2 #low pass filter (Hz)
highpass = F
highPassFreq = 0.01 #high pass filter (Hz)

check_minmax_plots = F 
saveSclPlotsToFile = F


shock_triggers = c(10,11,12) #Vector with the names of the US triggers
condition_triggers = c(2,3,4,5,6,7,8,9)

sample_rate = 1000  #samplerate after export
trial_length = 12  #length of trial in s

#Min-Max-SCR Analysis:

ucrMinWindow = c(.5, 3) #within how many seconds after UCS must the detected minimum appear?
ucrMin = .1  #minimum amplitude of an individual unconditioned response (in µS)
ucrRiseMax = 6 #latest Peak of UCR in s

crMin = 0.02      #minimum amplitude of CS-SCR
crMinWindow = c(.8,3) #within how many seconds after CS must the detected minimum appear? (Bouton et al., 2012)
crRiseMax = 6 #latest Peak of CR in s (Bouton et al., 2012)

#SCL Analysis:

scl_bins = T #calculate SCL in bins after CS onset
baselineWindow = c(-.5,0) #correct for Baseline in this time window
sclWindow = c(0,10) #calculate SCL in this time window


#First and second responses (FIR & SIR) analysis: 

fir_algorithm = F

eirWindow = c(1,10) #entire interval response window
firWindow = c(1,5) #first interval response window
sirWindow = c(5,10) #second interval response window
eirMinWindow = c(-.5,3) #minimun interval window for baseline correction

{ # Functions ---------------------------------------------------------------
  sample.down = function(signal, conversionRate) {
    suppressWarnings(matrix(signal, ncol=conversionRate, byrow=T)) %>% apply(1, mean) %>% 
      head(-1) #discard last sample because it gets distorted by zero padding
  }
  
  closestRepresentative = function(x, values, returnIndices=F) {
    indices = x %>% sapply(function(x, data) which.min(abs(data-x)), data=values)
    if (returnIndices) return(indices)
    return(values[indices])
  }
  
  localMaxima = function(signal, integerMode=F) {
    result = signal %>% c(-ifelse(integerMode, .Machine$integer.max, Inf), .) %>% #prepend smallest value to allow that first sample point is also local maximum
      diff() %>% {. > 0} %>% #which sample points are rising compared to their predecessor?
      rle() %>% .$lengths %>% #for how long are segments rising consecutively?
      cumsum() %>% #cumulative sum to find indices of local extrema
      {.[seq.int(from=1, to=length(.), by=2)]} #only take every other point
    if (signal[[1]] == signal[[2]]) result = result[-1]
    
    return(result)
  }
  
  localMinima = function(signal, integerMode=F) {
    result = signal %>% c(ifelse(integerMode, .Machine$integer.max, Inf), .) %>% #prepend biggest value to allow that first sample point is also local maximum
      diff() %>% {. > 0} %>% #which sample points are rising compared to their predecessor?
      rle() %>% .$lengths %>% #for how long are segments rising consecutively?
      cumsum() %>% #cumulative sum to find indices of local extrema
      {.[seq.int(from=1, to=length(.), by=2)]} #only take every other point
    if (signal[[1]] == signal[[2]]) result = result[-1]
    
    return(result)
  }
  
  inflectionPoints = function(signal, integerMode=F, upOnly=T) {
    extrema = signal %>% c(ifelse(integerMode, .Machine$integer.max, Inf), .) %>% diff()
    curvature = extrema %>% c(ifelse(integerMode, .Machine$integer.max, Inf)) %>% #for curvature, append biggest value (first and last sample point have no defined curvature)
      diff()
    result = curvature %>% #{. == 0} %>% which() #not working for discrete graphs
      {. < 0} %>% #which sample points are curved downwards (compared to adjacent)?
      rle() #run length of consecutive curvature
    
    if (upOnly) keep = result$values
    
    result = result %>% .$lengths %>% #for how long are segments curved upwards consecutively?
      cumsum() #cumulative sum to find indices of inflection points
    
    if (upOnly) result = result[keep]
    
    return(result)
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


edas_cr_list = list()
edas_list = list()
edas_df_list = list()
eda_maxima_list = list()
eda_minima_list = list()
eda_inflection_list = list()
eda_df = tibble()
eda_unified = tibble()

filemat = list.files("../Physio/Raw/", pattern="*.txt") # Das sollte der Ordner sein, in dem deine ganzen Files liegen. Pro VP ein File.

eda_ucr = tibble(subject = filemat %>% sub("\\..*","", .), ucr = 0, valid = 0, lat = 0, rise = 0) %>% 
  separate(subject,c("Exp","ID","tech","tech2"),sep="_",remove=T) %>%
  select(ID,ucr,valid,lat,rise)
eda_cr = tibble(subject = filemat %>% sub("\\..*","", .), cr = 0, valid = 0, lat = 0, rise = 0) %>% 
  separate(subject,c("Exp","ID","tech","tech2"),sep="_",remove=T) %>%
  select(ID,cr,valid,lat,rise)


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
  

for (subject_inmat in filemat){ #Jetzt rechnet er den Spaß für jedes File durch, wenn du einzelne Files berechnen willst, mach es nicht mit der for loop sondern kommentier die nächste zeile wieder ein
  
  # subject_inmat = filemat[[5]]
  
  eda <- read.csv(paste0("../Physio/Raw/",subject_inmat), sep="\t", header = F) %>% #read export
    rename(EDA = V1, ECG = V2, ImgAcq = V3, ImgTest = V4, Shock = V5, Reward = V6, NoFeedback = V7) %>%
    mutate(sample = 1:n())
  
  firstcue <- eda %>% filter(ImgAcq == 5) %>% head(1) %>% .$sample
  lastcue <- eda %>% filter(ImgTest == 5) %>% tail(1) %>% .$sample
  recordend <- eda$sample %>% tail(1)
  
  eda <- eda %>% filter(sample >= firstcue - 2000 & sample < lastcue + 20000) %>%
    mutate(ImgAcq = ifelse(ImgAcq == 5,1,0),
           ImgTest = ifelse(ImgTest == 5,1,0),
           Shock = ifelse(Shock == 5,10,0),
           Reward = ifelse(Reward == 5,11,0),
           NoFeedback = ifelse(NoFeedback == 5,12,0),
           trigger = ImgAcq + ImgTest + Shock + Reward + NoFeedback,
           sample = 1:n(),
           time = sample /1000) %>% 
    select(sample,time, EDA, ECG, trigger) %>%
    mutate(trigger = ifelse(trigger == lag(trigger), 0, trigger)) %>%
    fill(trigger,.direction = "up") %>%
    fill(trigger,.direction = "down")
    
  filename =  filemat[filemat == subject_inmat] %>% str_remove(.,pattern=".txt")# %>% str_remove(.,pattern="SiCP_")
  
  #sum(eda$trigger == 2);sum(eda$trigger == 3);sum(eda$trigger == 4);sum(eda$trigger == 5);sum(eda$trigger == 6);sum(eda$trigger == 7);
  #sum(eda$trigger == 8);sum(eda$trigger == 9);sum(eda$trigger == 10);sum(eda$trigger == 11);sum(eda$trigger == 12); #count triggers and check
  
  print(paste0(filename, " data reading complete!"))
  
  
  if (downsampling) {
    #downsample (all columns)
    conversion = round(sample_rate / sample_rate_new)
    eda_downsampled = data.frame(time=sample.down(eda$time,conversion),
                                 EDA=sample.down(eda$EDA,conversion),
                                 trigger=0) %>% 
      mutate(sample=1:n()) %>% select(sample, everything())
    
    #calculate closest position for trigger onsets in downsampled time
    triggers_time_old = eda$time[eda$trigger != 0]
    triggers_indices_new = triggers_time_old %>% closestRepresentative(eda_downsampled$time, returnIndices = T) #for each old trigger time, find index of closest existing downsampled time
    eda_downsampled$trigger[triggers_indices_new] = eda$trigger[eda$trigger != 0] #inject conditions as triggers of onsets
    
    
    for (row in 1:(nrow(eda_downsampled)-1)){ if (eda_downsampled$trigger[row] == eda_downsampled$trigger[row+1]) {eda_downsampled$trigger[row+1]=0}} #remove extra shock-markers
    
    # eda_downsampled %>% ggplot(aes(x=time, y=eda)) + geom_line() + geom_vline(xintercept = eda$time[eda$trigger!=0])
    if(check_ds_plots){
      eda %>% ggplot(aes(x=time, y=EDA)) + geom_line() + geom_line(data=eda_downsampled, color="red",linetype=4) + 
        geom_vline(xintercept = eda$time[eda$trigger!=0],alpha=0.1)  %>% print()
    }  
    
    eda = eda_downsampled; rm(eda_downsampled)
    print("Downsampling complete!")
  }
  
  
  # Calculate time difference of triggers because of missing triggers in GCA study --> this way we can adapt the condition labels
  trigger_diff <- eda %>% filter(trigger == 1) %>%
    mutate(sdiff = lead(time, 1) - time) %>%
    mutate(use.trial = TRUE) %>%
    mutate(trial = row_number())
  
  if (nrow(trigger_diff) != 152) {
    
    second_highest = sort(trigger_diff$sdiff,partial=length(trigger_diff$sdiff)-2)[(length(trigger_diff$sdiff))-2]
    
    trigger_diff <- trigger_diff %>%
      mutate(rownum = row_number()) %>%
      bind_rows(., filter(., (sdiff > 16 ) & (sdiff < second_highest)) %>%
                  mutate(use.trial = FALSE, rownum = rownum+.5)) %>%
      arrange(rownum) %>%
      mutate(trial = row_number())
  }
  
  filename =  filemat[filemat == subject_inmat] %>% str_remove(.,pattern=".txt")
  trigger_mat_subj <- trigger_mat %>% filter(subject == filename)
  
  trigger_mat_subj <- trigger_mat_subj %>%
    left_join(trigger_diff %>% select(trial, use.trial), by=c("trial")) %>%
    mutate(trigger = if_else((filename == "gca_14") & (trial > 112), 0, trigger)) # exclude test trials in VP 14 because one is missing, but we do not know which one
  
  trigger_mat_subj <- trigger_mat_subj %>%
    filter(use.trial)
  
  # if also the corrected triggers are not correct, then drop this subject
  if ((nrow(trigger_diff) == 152) |  (filename == "gca_14")) {
    # Include trigger informationwar
    eda$trigger[eda$trigger == 1] <-  trigger_mat_subj %>% .$trigger
  } else {
    print(paste0(filename, ': ', sum((eda$trigger > 1) & (eda$trigger < 10)), ' trials'))
    next
  }
  
  print(paste0(filename, ': ', sum((eda$trigger > 1) & (eda$trigger < 10)), ' trials'))
  if (sum((eda$trigger > 1) & (eda$trigger < 10)) != 152) { warning(paste0("Warning: Too few/many Acq triggers"))}
  

  #smoothing / filtering
  if (lowpass) {
    reps = 100 #filters can produce artifacts on the edges => fill edges with first and last values
    eda_filtered = c(rep(first(eda$EDA), reps), eda$EDA, rep(last(eda$EDA), reps)) %>% 
      signal::filtfilt(signal::butter(2, lowPassFreq/(ifelse(downsampling, sample_rate_new, sample_rate)/2)), .) #low-pass filter
    
    if(highpass){
      eda_filtered = eda_filtered %>%  # high-pass filter
        signal::filtfilt(signal::butter(2, type="high", highPassFreq/(ifelse(downsampling, sample_rate_new, sample_rate)/2)), .)  
    }
    
    eda_filt = eda
    eda_filt$EDA = eda_filtered[(reps+1):(length(eda_filtered)-reps)]; rm(eda_filtered)
    eda = eda_filt
    print("lowpass filtering complete!")
    
  }
  
  eda = eda %>% mutate(us = case_when(trigger %in% shock_triggers ~ T, TRUE ~ F))
  eda_vp = eda %>% filter(trigger != 0) %>% select(-EDA) %>% #One trial per row per VP
    mutate(trial = 1:n(), condition = trigger,
           us = us, us_prior = c(FALSE, lag(us)[-1]),
           time_start = time, time_end = time + trial_length,
           sample_start = sample, 
           sample_end =  {sample+(trial_length*ifelse(downsampling,sample_rate_new,sample_rate))} %>% round()) %>% 
    select(-trigger, -time, -sample) %>% select(trial, condition, everything())
  
  
  edas_df_list[[filename]] = eda_vp
  edas_list[[filename]] = eda
  
  # no nesting anymore
  #eda_vp$eda = eda_vp$time_start %>% lapply(function(start) 
  #  eda %>% select(-us) %>%  filter(time >= start, time < start + trial_length)) 
  
  eda_maxima = tibble(sample = localMaxima(eda$EDA)) %>% 
    mutate(n = 1:n(), time=eda$time[sample], EDA= eda$EDA[sample]) %>% select(n, everything())
  eda_minima = tibble(sample = localMinima(eda$EDA)) %>% 
    mutate(n = 1:n(), time=eda$time[sample], EDA = eda$EDA[sample]) %>% select(n, everything())
  eda_inflection = tibble(sample = inflectionPoints(eda$EDA)) %>% 
    mutate(n = 1:n(), time=eda$time[sample], EDA = eda$EDA[sample]) %>% select(n, everything())
  
  eda_maxima_list[[filename]] = eda_maxima
  eda_minima_list[[filename]] = eda_minima
  eda_inflection_list[[filename]] = eda_inflection
  
  if (check_minmax_plots || saveSclPlotsToFile) {
    sclPlot = eda %>% ggplot(aes(x=time, y=EDA)) + 
      geom_vline(xintercept = eda$time[eda$trigger!=0], color=ifelse(eda_vp$us, "orange", "black")) +
      geom_line() + 
      #geom_point(data=eda.inflection, color="purple") + 
      geom_point(data=eda_maxima, color="red") + 
      geom_point(data=eda_minima, color="blue") + 
      ggtitle(filename)
    if (check_minmax_plots) sclPlot %>% print()
    if (saveSclPlotsToFile) sclPlot %>% ggsave(paste0("../plots/EDA/", filename, ".png"), plot=., device="png", width=1920/50, height=1080/300, dpi=300)
  }
  
  # UCR Scoring    
  eda_vp_ucr = eda_vp %>% filter(us == TRUE)
  
  for (t in 1:nrow(eda_vp_ucr)) {
    
    #t = 1
    start = eda_vp_ucr$time_start[t]
    eda_minima_trial = eda_minima %>% filter(time >= start + min(ucrMinWindow),
                                             time <= start + max(ucrMinWindow))
    eda_maxima_trial = eda_maxima %>% filter(time > suppressWarnings(min(eda_minima_trial$time))) %>% 
      head(nrow(eda_minima_trial)) #same amount of maxima as minima
    
    #if no minimum found within time range use deflection point as minimum
    if (nrow(eda_minima_trial)==0) { 
      eda_minima_trial = eda_inflection %>% filter(time >= start + min(ucrMinWindow),
                                                   time <= start + max(ucrMinWindow))
      eda_maxima_trial = eda_maxima %>% filter(time > suppressWarnings(min(eda_minima_trial$time))) %>% 
        head(nrow(eda_minima_trial)) #same amount of maxima as minima
      
      while (nrow(eda_maxima_trial)>1) {
        latestMaxTimes = eda_maxima_trial %>% tail(2) %>% .$time
        if (eda_minima_trial %>% filter(time < max(latestMaxTimes), #check if any inflection point is between last two maxima (may be false if two inflection points before earlier max)
                                        time > min(latestMaxTimes)) %>% 
            nrow() %>% {. > 0}) break
        else { #if no inflection point between last two maxima
          eda_minima_trial = eda_minima_trial %>% head(-1) #delete last inflection point
          eda_maxima_trial = eda_maxima_trial %>% head(nrow(eda_minima_trial)) #same amount of maxima as minima
        }
      }
      
      #check if true min occurs between inflection point and corresponding max (with true min being too late)
      toKeep = c()
      if (nrow(eda_minima_trial) > 0) { 
        for (i in 1:nrow(eda_minima_trial)) {
          if (eda_minima %>% filter(time > eda_minima_trial[i, "time"], time < eda_maxima_trial[i, "time"]) %>% 
              nrow() %>% {. == 0}) toKeep = toKeep %>% c(i)
        }
        eda_minima_trial = eda_minima_trial[toKeep,]; eda_maxima_trial = eda_maxima_trial[toKeep,] #only keep inflection points that don't have a true minimum between themself and their max
      }
    }
    
    #calculate SCRs
    eda_minima_trial_renamed = eda_minima_trial; names(eda_minima_trial_renamed) = names(eda_minima_trial_renamed) %>% paste0(".min")
    eda_maxima_trial_renamed = eda_maxima_trial; names(eda_maxima_trial_renamed) = names(eda_maxima_trial_renamed) %>% paste0(".max")
    scr = bind_cols(eda_minima_trial_renamed, eda_maxima_trial_renamed)
    #names(scr) = names(scr) %>% {ifelse(endsWith(., "1"), gsub("1", ".max", ., fixed=T), paste(., "min", sep="."))}
    scr = scr %>% mutate(SCR = EDA.max - EDA.min,
                         Lat = time.min - start,
                         Rise = time.max - start) %>% select(SCR, Lat, Rise, everything()) %>% 
      filter(SCR == max(SCR)) #take maximum SCR
    
    if (nrow(scr) > 1){scr <- scr[1,]}
    
    if (nrow(scr) == 0) {scr = scr %>% bind_rows(tibble(SCR = 0, Lat = NA, Rise = NA)) #if no valid SCR -> 0
    } else if (scr$SCR <= ucrMin) {scr = scr %>% mutate(SCR = 0, Lat = NA, Rise = NA)
    } else if (scr$Rise > ucrRiseMax ) {scr = scr %>% mutate(SCR = 0, Lat = NA, Rise = NA)}
    
    scr_out = scr %>% select(-c(n.min,sample.min,n.max,sample.max))
    eda_vp_ucr[t, names(scr_out)] = scr_out
    
  }  
  
  eda_ucr[eda_ucr$ID==filename, -1] = eda_vp_ucr %>% mutate(valid = SCR > ucrMin) %>% 
    summarise(scr = mean(SCR, na.rm=T), 
              valid = mean(valid),
              lat = mean(Lat, na.rm=T),
              rise = mean(Rise, na.rm=T))
  
  unified = eda_vp_ucr %>% mutate(EDA = time_start %>% lapply(function(start) 
    eda %>% filter(time >= start,
                   time <= start + 8) %>% select(-trigger)))
  
  for (t in 1:nrow(unified)) {
    unified$EDA[[t]] = unified$EDA[[t]] %>% 
      mutate(trial = t,
             sample = sample - min(sample) + 1,
             time = time - min(time),
             condition = unified$condition[[t]])
  }
  
  minmax_unified = unified %>% select(trial,EDA.min,time.min,EDA.max,time.max,time_start) %>% 
    mutate(time.min = time.min - time_start,time.max = time.max - time_start,trial = as.factor(1:n())) # %>% filter(time.max < 8)
  
  if (ucr_plots) {
    unified$EDA %>% bind_rows() %>% 
      mutate(trial = as.factor(trial)) %>%
      mutate(across('condition', str_replace_all, rep_str)) %>%
      {ggplot(., aes(x=time, y=EDA, color=trial)) + facet_wrap(vars(condition)) +
          geom_vline(xintercept=0, color="black",linetype="solid") + #zero 
          geom_vline(xintercept=ucrMinWindow, color="grey") + #borders of min/max scoring
          geom_path() + scale_color_viridis_d() + scale_fill_viridis_d() +
          # geom_point(data=minmax_unified,aes(x=time.min,y=EDA.min),shape=25,size=1) +
          # geom_point(data=minmax_unified,aes(x=time.max,y=EDA.max),shape=24,size=1) +
          ggtitle(as.integer(substr(filename, 5, 6))) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))} %>%
      ggsave(paste0("../plots/EDA/UCS/", filename, ".png"), plot=., device="png", width=1920/150, height=1080/150, dpi=300)
  }
  
  print(paste0("UCRs: ", as.integer(substr(filename, 5, 6)), " of ", length(filemat), " files processed!"))
  
  eda_vp_cr = eda_vp %>% filter(condition %in% condition_triggers) #filter(us == F)
  
  for (t in 1:nrow(eda_vp_cr)) {
    
    #t = 1
    start = eda_vp_cr$time_start[t]
    
    baseline_trial = eda %>% filter(time <= start + max(baselineWindow),
                                    time >= start + min(baselineWindow)) %>% 
      summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
    
    
    scl_trial = eda %>% filter(time <= start + max(sclWindow),
                               time >= start + min(sclWindow)) %>% 
      summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
 
    if(scl_bins){
      
      scl_0 <- eda %>% filter(time >= start, time <= start+0.5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_1 <- eda %>% filter(time >= start+0.5, time <= start+1) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_2 <- eda %>% filter(time >= start+1, time <= start+1.5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_3 <- eda %>% filter(time >= start+1.5, time <= start+2) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_4 <- eda %>% filter(time >= start+2, time <= start+2.5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_5 <- eda %>% filter(time >= start+2.5, time <= start+3) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_6 <- eda %>% filter(time >= start+3, time <= start+3.5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_7 <- eda %>% filter(time >= start+3.5, time <= start+4) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_8 <- eda %>% filter(time >= start+4, time <= start+4.5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_9 <- eda %>% filter(time >= start+4.5, time <= start+5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_10 <- eda %>% filter(time >= start+5, time <= start+5.5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_11 <- eda %>% filter(time >= start+5.5, time <= start+6) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_12 <- eda %>% filter(time >= start+6, time <= start+6.5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_13 <- eda %>% filter(time >= start+6.5, time <= start+7) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_14 <- eda %>% filter(time >= start+7, time <= start+7.5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_15 <- eda %>% filter(time >= start+7.5, time <= start+8) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_16 <- eda %>% filter(time >= start+8, time <= start+8.5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_17 <- eda %>% filter(time >= start+8.5, time <= start+9) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_18 <- eda %>% filter(time >= start+9, time <= start+9.5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_19 <- eda %>% filter(time >= start+9.5, time <= start+10) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_2_10 <- eda %>% filter(time >= start+2, time <= start+10) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      }
     
    if(fir_algorithm){
      eir_max_searching = T
      upperLimit <- start + max(eirWindow)
      lowerLimit <- start + min(eirWindow)
      
      while (eir_max_searching){ 
        
        eir_max = eda %>% 
          filter(time <= upperLimit &
                   time >= lowerLimit) %>%
          filter(EDA == max(EDA)) %>% head(1)
        
        
        if(eir_max$time > lowerLimit & eir_max$time < upperLimit){
          eir_max_searching <- F
        }    
        
        if (isTRUE(all.equal(eir_max$time,upperLimit))){
          upperLimit <- round(upperLimit - 1/sample_rate,2)
        } else if(isTRUE(all.equal(eir_max$time,lowerLimit))){
          lowerLimit <- round(lowerLimit + 1/sample_rate,2)
        }
        
        if(upperLimit == lowerLimit){
          eir_max$EDA <- NA
          eir_max_searching <- F
        }
      }
    
    
    
    fir_max_searching = T
    upperLimit <- start + max(firWindow)
    lowerLimit <- start + min(firWindow)
    
    while (fir_max_searching){ 
      
      fir_max = eda %>% 
        filter(time <= upperLimit &
                 time >= lowerLimit) %>%
        filter(EDA == max(EDA)) %>% head(1)
      
      
      if(fir_max$time > lowerLimit & fir_max$time < upperLimit){
        fir_max_searching <- F
      }    
      
      if (isTRUE(all.equal(fir_max$time,upperLimit))){
        upperLimit <- round(upperLimit - 1/sample_rate,2)
      } else if(isTRUE(all.equal(fir_max$time,lowerLimit))){
        lowerLimit <- round(lowerLimit + 1/sample_rate,2)
      }
      
      if(upperLimit == lowerLimit){
        fir_max$EDA <- NA
        fir_max_searching <- F
      }
    }
    
    
    sir_max_searching = T
    upperLimit <- start + max(sirWindow)
    lowerLimit <- start + min(sirWindow)
    
    while (sir_max_searching){ 
      
      sir_max = eda %>% 
        filter(time <= upperLimit &
                 time >= lowerLimit) %>%
        filter(EDA == max(EDA)) %>% head(1)
      
      
      if(sir_max$time > lowerLimit & sir_max$time < upperLimit){
        sir_max_searching <- F
      }    
      
      if (isTRUE(all.equal(sir_max$time,upperLimit))){
        upperLimit <- round(upperLimit - 1/sample_rate,2)
      } else if(isTRUE(all.equal(sir_max$time,lowerLimit))){
        lowerLimit <- round(lowerLimit + 1/sample_rate,2)
      }
      
      if(upperLimit == lowerLimit){
        sir_max$EDA <- NA
        sir_max_searching <- F
      }
    }
    
    eir_min = eda %>% 
      filter(time <= start + max(eirMinWindow) &
               time >= start + min(eirMinWindow)) %>%
      filter(EDA == min(EDA)) %>% head(1)
    
    }
    
    eda_minima_trial = eda_minima %>% filter(time >= start + min(crMinWindow),
                                             time <= start + max(crMinWindow)) 
    eda_maxima_trial = eda_maxima %>% filter(time > suppressWarnings(min(eda_minima_trial$time))) %>% 
      head(nrow(eda_minima_trial)) #same amount of maxima as minima
    
    if (nrow(eda_maxima_trial) < nrow(eda_minima_trial)) { #if recording ended before next max
      eda_maxima_trial = eda_maxima_trial %>% bind_rows(
        eda %>% filter(time >= max(eda_minima_trial$time)) %>% filter(EDA==max(EDA)) %>% mutate(n = 0) %>% select(n, everything()) %>% select(-trigger)
      )
    }
    
    
    #if no minimum found within time range use deflection point as minimum
    if (nrow(eda_minima_trial)==0) { 
      eda_minima_trial = eda_inflection %>% filter(time >= start + min(crMinWindow),
                                                   time <= start + max(crMinWindow))
      eda_maxima_trial = eda_maxima %>% filter(time > suppressWarnings(min(eda_minima_trial$time))) %>% 
        head(nrow(eda_minima_trial)) #same amount of maxima as minima
      if (nrow(eda_maxima_trial) < nrow(eda_minima_trial)) { #if recording ended before next max
        eda_maxima_trial = eda_maxima_trial %>% bind_rows(
          eda %>% filter(time >= max(eda_minima_trial$time)) %>% filter(EDA==max(EDA)) %>% mutate(n = 0) %>% select(n, everything()) %>% select(-trigger)
        )
      }
      
      #check if mins and maxs alternate (in case inflection points were used)
      while (nrow(eda_maxima_trial)>1) {
        latestMaxTimes = eda_maxima_trial %>% tail(2) %>% .$time
        if (eda_minima_trial %>% filter(time < max(latestMaxTimes), #check if any inflection point is between last two maxima (may be false if two inflection points before earlier max)
                                        time > min(latestMaxTimes)) %>% 
            nrow() %>% {. > 0}) break
        else { #if no inflection point between last two maxima
          eda_minima_trial = eda_minima_trial %>% head(-1) #delete last inflection point
          eda_maxima_trial = eda_maxima_trial %>% head(nrow(eda_minima_trial)) #same amount of maxima as minima
        }
      }
      
      #check if true min occurs between inflection point and corresponding max (with true min being too late)
      toKeep = c()
      if (nrow(eda_minima_trial) > 0) { 
        for (i in 1:nrow(eda_minima_trial)) {
          if (eda_minima %>% filter(time > eda_minima_trial[i, "time"], time < eda_maxima_trial[i, "time"]) %>% 
              nrow() %>% {. == 0}) toKeep = toKeep %>% c(i)
        }
        eda_minima_trial = eda_minima_trial[toKeep,]; eda_maxima_trial = eda_maxima_trial[toKeep,] #only keep inflection points that don't have a true minimum between themself and their max
      }
    }
    
    #calculate SCRs
    eda_minima_trial_renamed = eda_minima_trial; names(eda_minima_trial_renamed) = names(eda_minima_trial_renamed) %>% paste0(".min")
    eda_maxima_trial_renamed = eda_maxima_trial; names(eda_maxima_trial_renamed) = names(eda_maxima_trial_renamed) %>% paste0(".max")
    scr = bind_cols(eda_minima_trial_renamed, eda_maxima_trial_renamed)
    #names(scr) = names(scr) %>% {ifelse(endsWith(., "1"), gsub("1", ".max", ., fixed=T), paste(., "min", sep="."))}
    scr = scr %>% mutate(SCR = EDA.max - EDA.min,
                         Lat = time.min - start,
                         Rise = time.max - start,
                         Baseline = baseline_trial,
                         SCL = scl_trial) %>% select(SCR, Lat, Rise, SCL, Baseline, everything()) %>% 
      filter(SCR == max(SCR)) %>% #take maximum SCR
      head(1) #if two or more equal SCRs, take only the first
    
    if (nrow(scr) == 0) {scr = scr %>% bind_rows(tibble(SCR = 0, Lat = NA, Rise = NA, Baseline = baseline_trial, SCL = scl_trial))
    } else if (scr$SCR <= crMin) {scr = scr %>% mutate(SCR = 0, Lat = NA, Rise = NA, Baseline = baseline_trial, SCL = scl_trial)
    } else if (scr$Rise > crRiseMax) {scr = scr %>% mutate(SCR = 0, Lat = NA, Rise = NA, Baseline = baseline_trial, SCL = scl_trial)}
    
    
    scr_out = scr %>% select(1:5,8,9,12,13)

    if(scl_bins){
      scr_out <- scr_out %>%
        mutate(
          scl_0 = scl_0 - Baseline,
          scl_1 = scl_1 - Baseline,
          scl_2 = scl_2 - Baseline,
          scl_3 = scl_3 - Baseline,
          scl_4 = scl_4 - Baseline,
          scl_5 = scl_5 - Baseline,
          scl_6 = scl_6 - Baseline,
          scl_7 = scl_7 - Baseline,
          scl_8 = scl_8 - Baseline,
          scl_9 = scl_9 - Baseline,
          scl_10 = scl_10 - Baseline,
          scl_11 = scl_11 - Baseline,
          scl_12 = scl_12 - Baseline,
          scl_13 = scl_13 - Baseline,
          scl_14 = scl_14 - Baseline,
          scl_15 = scl_15 - Baseline,
          scl_16 = scl_16 - Baseline,
          scl_17 = scl_17 - Baseline,
          scl_18 = scl_18 - Baseline,
          scl_19 = scl_19 - Baseline,
          scl_2_10 = scl_2_10 - Baseline,
        )
    }
    
    if(fir_algorithm){
      scr_out = scr %>% select(1:5,8,9,12,13) %>%
        mutate(EDA.eir = eir_max$EDA,
               EDA.eir.lat = eir_max$time - start,
               EDA.fir = fir_max$EDA,
               EDA.fir.lat = fir_max$time - start,
               EDA.sir = sir_max$EDA,
               EDA.sir.lat = sir_max$time - start,
               EDA.eir.min = eir_min$EDA,
               EDA.eir.min.lat = eir_min$time - start
        )
    }
    
    eda_vp_cr[t, names(scr_out)] = scr_out
    #readline(prompt="Press [enter] to continue")
  } 
  
  eda_vp_cr = eda_vp_cr %>% group_by(condition) %>% mutate(trial_condition = row_number()) %>% ungroup()
  edas_cr_list[[filename]] = eda_vp_cr
  eda_df  = rbind(eda_df,eda_vp_cr %>% mutate(ln_cr = log(SCR+1),ID=filename))
 
  print("trialwise parameterization complete!")
  
  unified = eda_vp_cr %>% 
    mutate(EDA = time_start %>% lapply(function(start)
      eda %>% filter(time >= start + min(baselineWindow),
                     time <= start + trial_length) %>% select(-trigger)))
  
  
  for (t in 1:nrow(unified)) {
    # t = 1
    unified$EDA[[t]] = unified$EDA[[t]] %>% 
      mutate(subject = as.integer(substr(filename, 5, 6)),
             trial = t,
             sample = sample - unified$sample_start[[t]],
             time = time - unified$time_start[[t]],
             samplepoint = sample-min(sample),
             condition = unified$condition[[t]])
    }
  
  minmax_unified = unified %>% select(trial,EDA.min,time.min,EDA.max,time.max,time_start, condition) %>% 
    mutate(time.min = time.min - time_start,time.max = time.max - time_start,trial = as.factor(1:n())) %>% 
    filter(time.max < 8) %>% 
    mutate(across('condition', str_replace_all, rep_str))
  
  if (cr_plots) {
    unified$EDA %>% bind_rows() %>% 
      mutate(across('condition', str_replace_all, rep_str)) %>%
      mutate(trial = as.factor(trial)) %>% 
      {ggplot(., aes(x=time, y=EDA, color=trial, fill=trial)) + facet_wrap(vars(condition)) +
          geom_vline(xintercept=c(ucrMinWindow,crRiseMax), color="blue",linetype="dashed") + #borders of min/max scoring
          geom_path() + scale_color_viridis_d() + scale_fill_viridis_d() +
          geom_point(data=minmax_unified,aes(x=time.min,y=EDA.min),shape=24,size=1)+
          geom_point(data=minmax_unified,aes(x=time.max,y=EDA.max),shape=24,size=1)+
          xlab("Time") + ylab("EDA") +
          ggtitle(as.integer(substr(filename, 5, 6))) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))} %>% 
      #print()
      ggsave(paste0("../plots/EDA/CS/", filename, ".png"), plot=., device="png", width=1920/150, height=1080/150, dpi=300)
    
    print("CR plotting complete!")
  }
  
  for (t in 1:nrow(unified)) {
    # t = 1
    unified$EDA[[t]] = unified$EDA[[t]] %>% 
      mutate(EDA = EDA - unified$Baseline[[t]])
  }
  
  unified <- unified %>% mutate(ID = filename)
  eda_unified = rbind(eda_unified,unified)
  
  
  print(paste0("CRs: ", as.integer(substr(filename, 5, 6)), " of ", length(filemat), " files processed!"))
  
} #end inmat for loop



# for (t in 1:nrow(eda_unified)) {
#   eda_unified$EDA[[t]] = eda_unified$EDA[[t]] %>% 
#     mutate(EDA_bl = EDA - eda_unified$Baseline[[t]], ID = eda_unified$ID[[t]])
# }
# 
# rm(eda, eda_vp)

# Read or Save eda_unified and eda_df for further processing

saveRDS(eda_unified,"EDA_unified.RData")
saveRDS(eda_df,"EDA_df.RData")

eda_unified <- readRDS("EDA_unified.RData")
eda_df <- readRDS("EDA_df.Rdata")

responder <- eda_df %>%
 mutate(scl_bl = abs(SCL-Baseline)) %>%
 group_by(ID) %>%
 summarise(scl_avg = mean(scl_bl)) %>% 
 filter(scl_avg >= .02) %>% .$ID

non_responder <- eda_df %>%
 mutate(scl_bl = abs(SCL-Baseline)) %>%
 group_by(ID) %>%
 summarise(scl_avg = mean(scl_bl)) %>% 
 filter(scl_avg < .02) %>% .$ID
