###############################################################################
# Gaze Contingent Avoidance Project
# Sabrina Gado & Yannik Stegmann

# Pupil Analyse

#packages 
library(signal)
library(zoo)
library(tidyverse)
library(stringr)
library(afex)
library(lme4)
library(lmerTest)

options(dplyr.summarise.inform = FALSE)

# movemean kann man rausnehmen, genauso wie z standardisierung.
# anschließend ergebnisse vergleichen - matthais hat auch eine 1s baseline


# options & plots
cr_plots = T


sample_rate = 100  #samplerate after export
trial_length = 12  #length of trial in s

downsampling = TRUE
sample_rate_new = 50 #for downsampling
check_ds_plots = F

standardizing = F

lowpass = T
lowPassFreq = 2 #low pass filter (Hz) entweder keinen oder 2Hz wie bei Matthias

baselineWindow = c(-0.5, 0)
# responseWindow = c(0, 8)

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
}


data = read.csv2("../Physio/pupil.txt", sep ="\t", dec = ",", na.strings=".") 

data <- data %>% select(-TRIAL_LABEL) %>%
  mutate(diameter = ifelse(is.na(RIGHT_PUPIL_SIZE_BIN), LEFT_PUPIL_SIZE_BIN, RIGHT_PUPIL_SIZE_BIN)) %>% 
  select(-LEFT_PUPIL_SIZE_BIN, -RIGHT_PUPIL_SIZE_BIN, -BIN_END_TIME) 

codes <- unique(data$RECORDING_SESSION_LABEL)

trigger_mat <- read.csv2("../Physio/Trigger/conditions.csv") %>%
  mutate(subject = sprintf("gca_%02d", subject),
         ID = as.numeric(substr(subject,5,6)),
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

pupils_list = list()
pupils_df_list = list()
pupil_df = tibble()
ga_unified = tibble()
code_times = c()


messages = read.csv2("../Gaze/messages.csv", sep =";") %>% 
  filter(grepl("ImageOnset", event))


for (subject_inmat in codes){ 
  
  # subject_inmat = codes[[1]]
  start.time <- Sys.time()
  
  vp = ifelse(is.na(as.numeric(substr(subject_inmat,20,21))),as.numeric(substr(subject_inmat,20,20)),as.numeric(substr(subject_inmat,20,21)))
  
  pupil = data %>% filter(RECORDING_SESSION_LABEL == subject_inmat) %>%
    rename(ID = RECORDING_SESSION_LABEL, trial = TRIAL_INDEX, bin_index = BIN_INDEX, bin_start = BIN_START_TIME) %>%
    mutate(ID = vp) %>% 
    mutate(sample = 1:n()) %>% 
    mutate(time = sample*0.01) 
  
  conditions = trigger_mat %>% filter(ID == unique(pupil$ID)) %>%
    select(ID, trial, condition, trigger)
  
  messages_vp = messages %>% 
    filter(subject == vp) %>%
    rename(ID = subject, event_time = time) %>% 
    left_join(. ,conditions, by=c("ID","trial")) %>% 
    select(trial, event_time, trigger)
  
  pupil = pupil %>% 
    left_join(messages_vp, join_by(trial, closest(bin_start >= event_time))) %>% 
    mutate(trigger = ifelse(is.na(lag(trigger, 1)) & trigger == lead(trigger, 1), trigger, NA)) %>%
    mutate(trigger = ifelse(is.na(trigger), 0, trigger)) %>% 
    select(-event_time)
  
  filename = unique(pupil$ID)
  
  print("data reading successful!")
  
  # interpolation of missing values
  pupil = pupil %>% 
    mutate(diameter = ifelse((diameter > mean(diameter, na.rm = T) - 3 * sd(diameter, na.rm = T)) & (diameter < mean(diameter, na.rm = T) + 3 * sd(diameter, na.rm = T)), diameter, NA))
  pupil = pupil %>% rowwise() %>% mutate(interpolated = ifelse(is.na(diameter),T,F)) %>% ungroup() #create column for interpolated values
  pupil$diameter[[1]] = pupil$diameter[[min(which(!is.na(pupil$diameter)))]] # remove leading NA
  pupil$diameter[[length(pupil$diameter)]] = pupil$diameter[[max(which(!is.na(pupil$diameter)))]] # remove trailing NA
  pupil$diameter = na.approx(pupil$diameter) #interpolate NA
  
  print("interpolating missing values done!")
  
  # artefact correction
  if (downsampling) {
    
    #downsample (all columns)
    conversion = round(sample_rate / sample_rate_new)
    pupil_downsampled = data.frame(time=sample.down(pupil$time,conversion),
                                   diameter=sample.down(pupil$diameter,conversion),
                                   interpolated = sample.down(pupil$interpolated,conversion),
                                   trigger=0) %>% 
      mutate(sample=1:n()) %>% select(sample, everything()) %>%
      mutate(interpolated = ifelse(interpolated < 0.5,F,T))
    
    #calculate closest position for trigger onsets in downsampled time
    triggers_time_old = pupil$time[pupil$trigger != 0]
    triggers_indices_new = triggers_time_old %>% closestRepresentative(pupil_downsampled$time, returnIndices = T) #for each old trigger time, find index of closest existing downsampled time
    pupil_downsampled$trigger[triggers_indices_new] = pupil$trigger[pupil$trigger != 0] #inject conditions as triggers of onsets
    
    # eda_downsampled %>% ggplot(aes(x=time, y=eda)) + geom_line() + geom_vline(xintercept = eda$time[eda$trigger!=0])
    if(check_ds_plots){
      pupil %>% ggplot(aes(x=time, y=diameter)) + geom_line() + geom_line(data=pupil_downsampled, color="red",linetype=4) + 
        geom_vline(xintercept = pupil$time[pupil$trigger!=0],alpha=0.1) 
    }  
    pupil = pupil_downsampled; rm(pupil_downsampled)
    print("downsampling done!")
  }  
  
  #smoothing / filtering
  if (lowpass) {
    reps = 100 #filters can produce artifacts on the edges => fill edges with first and last values
    diameter_filtered = c(rep(first(pupil$diameter), reps), pupil$diameter, rep(last(pupil$diameter), reps)) %>% 
      signal::filtfilt(signal::butter(2, lowPassFreq/(ifelse(downsampling, sample_rate_new, sample_rate)/2)), .) #low-pass filter
    
    pupil_filt = pupil
    pupil_filt$diameter = diameter_filtered[(reps+1):(length(diameter_filtered)-reps)]; rm(diameter_filtered)
    
    pupil = pupil_filt
    print("low pass filtering done!")
  }  
  
  pupil_vp = pupil %>%  select(-diameter) %>%  filter(trigger != 0) %>% ungroup() %>% #One trial per row per VP
    mutate(trial = 1:n(), condition = trigger,
           time_start = time, time_end = time + trial_length,
           sample_start = sample, 
           sample_end =  {sample+(trial_length*ifelse(downsampling,sample_rate_new,sample_rate))} %>% round()) %>% 
    select(-trigger, -time, -sample, -interpolated) %>% select(trial, condition, everything())
  
  for (row in nrow(pupil_vp)){ #count interpolated samples per trial
    pupil_vp$interpolated = pupil_vp$time_start %>% lapply(function(start)
      pupil %>% filter(time >= start-abs(min(baselineWindow)), time < start + trial_length) %>% .$interpolated %>% mean(.) %>% round(.,3))
    
    exclude = rep(TRUE,ifelse(downsampling, sample_rate_new/2, sample_rate/2)) #exclude if more than (sampling rate / 2) interpolated samplepoints (0.5s) in a row
    
    pupil_vp$valid = pupil_vp$time_start %>% lapply(function(start)
      pupil %>% filter(time >= start-abs(min(baselineWindow)), time < start + trial_length) %>% {!grepl(paste(exclude,collapse=""), paste(.$interpolated,collapse=""))})
    
    pupil_vp = pupil_vp %>% 
      mutate(valid = ifelse(interpolated < .25, valid, FALSE))  #exclude if more than 25% interpolated data
  }
  
  
  pupils_df_list[[filename]] = pupil_vp
  pupils_list[[filename]] = pupil
  
  
  unified = pupil_vp %>% mutate(diameter = time_start %>% lapply(function(start) 
    pupil %>% filter(time >= start - abs(min(baselineWindow)),
                     time <= start + trial_length) %>% select(-trigger)))
  
  
  for (t in 1:nrow(unified)) {
    unified$diameter[[t]] = unified$diameter[[t]] %>% 
      mutate(trial = t, 
             #time = time - (min(time)+abs(min(baselineWindow))), 
             #samplepoint = sample-min(sample),
             condition = unified$condition[[t]], 
             #diameter = diameter - unified$baseline[[t]] #unify starting time to allow overlap
      )}
  
  if (standardizing) {  
    unified$diameter <- unified$diameter %>% bind_rows() %>% 
      mutate(diameter = (diameter - mean(diameter))/sd(diameter)) %>%
      #z-stand %>%
      split(.$trial)
    
    print("z-standardizing done!")
    
  }
  
  
  # baseline and level scoring  
  test.time <- Sys.time()  
  pupil_diameter_binded <- unified$diameter %>% bind_rows()
  
  pupil_diameter_vp <- tibble()
  
  for (t in 1:nrow(pupil_vp)) {
    # t = 1
    start_time = pupil_vp$time_start[t]
    start_sample = pupil_vp$sample_start[t]
    trial_idx = pupil_vp$trial[t]
    
    
    baseline_trial = pupil_diameter_binded %>% filter(time <= start_time + max(baselineWindow),
                                                      time >= start_time + min(baselineWindow)) %>% 
      summarize(mean = mean(diameter,na.rm=T))
    
    baseline_trial_diameter = baseline_trial$mean
    
    diameter_level_trial = pupil_diameter_binded %>% filter(time >= start_time - 0.5 & trial == trial_idx) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_0 = pupil_diameter_binded %>% filter(time >= start_time, time <= start_time + 0.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_1 = pupil_diameter_binded %>% filter(time >= start_time + 0.5, time <= start_time + 1) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_2 = pupil_diameter_binded %>% filter(time >= start_time + 1, time <= start_time + 1.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_3 = pupil_diameter_binded %>% filter(time >= start_time + 1.5,time <= start_time + 2) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_4 = pupil_diameter_binded %>% filter(time >= start_time + 2, time <= start_time + 2.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_5 = pupil_diameter_binded %>% filter(time >= start_time + 2.5,time <= start_time + 3) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_6 = pupil_diameter_binded %>% filter(time >= start_time + 3, time <= start_time + 3.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_7 = pupil_diameter_binded %>% filter(time >= start_time + 3.5,time <= start_time + 4) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_8 = pupil_diameter_binded %>% filter(time >= start_time + 4, time <= start_time + 4.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_9 = pupil_diameter_binded %>% filter(time >= start_time + 4.5, time <= start_time + 5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_10 = pupil_diameter_binded %>% filter(time >= start_time + 5,time <= start_time + 5.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_11 = pupil_diameter_binded %>% filter(time >= start_time + 5.5, time <= start_time + 6) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_12 = pupil_diameter_binded %>% filter(time >= start_time + 6,time <= start_time + 6.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_13 = pupil_diameter_binded %>% filter(time >= start_time + 6.5, time <= start_time + 7) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_14 = pupil_diameter_binded %>% filter(time >= start_time + 7,time <= start_time + 7.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_15 = pupil_diameter_binded %>% filter(time >= start_time + 7.5, time <= start_time + 8) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_16 = pupil_diameter_binded %>% filter(time >= start_time + 8,time <= start_time + 8.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_17 = pupil_diameter_binded %>% filter(time >= start_time + 8.5, time <= start_time + 8) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_18 = pupil_diameter_binded %>% filter(time >= start_time + 9,time <= start_time + 9.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_19 = pupil_diameter_binded %>% filter(time >= start_time + 9.5, time <= start_time + 10) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    
    dilations = tibble(dilation = diameter_level_trial - baseline_trial,
                       dilation_0 = diameter_level_trial_1 - baseline_trial,
                       dilation_1 = diameter_level_trial_1 - baseline_trial,
                       dilation_2 = diameter_level_trial_2 - baseline_trial,
                       dilation_3 = diameter_level_trial_3 - baseline_trial,
                       dilation_4 = diameter_level_trial_4 - baseline_trial,
                       dilation_5 = diameter_level_trial_5 - baseline_trial,
                       dilation_6 = diameter_level_trial_6 - baseline_trial,
                       dilation_7 = diameter_level_trial_7 - baseline_trial,
                       dilation_8 = diameter_level_trial_8 - baseline_trial,
                       dilation_9 = diameter_level_trial_9 - baseline_trial,
                       dilation_10 = diameter_level_trial_10 - baseline_trial,
                       dilation_11 = diameter_level_trial_11 - baseline_trial,
                       dilation_12 = diameter_level_trial_12 - baseline_trial,
                       dilation_13 = diameter_level_trial_13 - baseline_trial,
                       dilation_14 = diameter_level_trial_14 - baseline_trial,
                       dilation_15 = diameter_level_trial_15 - baseline_trial,
                       dilation_16 = diameter_level_trial_12 - baseline_trial,
                       dilation_17 = diameter_level_trial_13 - baseline_trial,
                       dilation_18 = diameter_level_trial_14 - baseline_trial,
                       dilation_19 = diameter_level_trial_15 - baseline_trial,
                       baseline = baseline_trial,
                       mean_diameter = diameter_level_trial) %>% select(dilation:dilation_19, baseline, mean_diameter, everything()) 
    
    
    pupil_vp[t, names(dilations)] = dilations
  } 
  end.time <- Sys.time() - test.time
  
  print("0.5s bin averaging done!")
  print(end.time)
  
  pupil_vp = pupil_vp %>% group_by(condition) %>% mutate(trial_condition = row_number()) %>% ungroup()
  pupils_df_list[[filename]] = pupil_vp
  
  
  
  for (t in 1:nrow(unified)) {
    unified$diameter[[t]] = unified$diameter[[t]] %>% 
      mutate(time = time - (min(time)+abs(min(baselineWindow))), 
             samplepoint = sample-min(sample),
             diameter = diameter - pupil_vp$baseline[[t]] #unify starting time to allow overlap
      )}
  
  if (cr_plots) {
    
    unified$diameter %>% bind_rows() %>% 
      mutate(trial = as.factor(trial)) %>% 
      mutate(across('condition', str_replace_all, rep_str)) %>% 
      {ggplot(., aes(x=time-0.5, y=diameter, color=trial)) + facet_wrap(vars(condition)) +
          geom_vline(xintercept=0, color="black",linetype="solid") + #zero 
          geom_path() + scale_color_viridis_d() +
          xlab("Time") + ylab("Diameter") +
          ggtitle(filename) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))} %>% 
      #print()
      ggsave(paste0("../plots/Pupil/subject_level/", filename, ".png"), plot=., device="png", width=1920/150, height=1080/150, dpi=300)
    
    
    
    unified$diameter %>% bind_rows() %>% 
      mutate(across('condition', str_replace_all, rep_str)) %>% 
      mutate(condition = as.factor(condition)) %>% 
      group_by(condition,samplepoint) %>% 
      summarise(diameter = mean(diameter), time = mean(time)) %>%
      {ggplot(., aes(x=time-0.5, y=diameter, color=condition, group=condition,linetype=condition)) +
          geom_vline(xintercept=0, color="black",linetype="solid") + #zero 
          geom_path() + scale_color_viridis_d() +
          xlab("Time") + ylab("Diameter") +
          ggtitle(filename) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))} %>% 
      #print()
      ggsave(paste0("../plots/Pupil/subject_level/", filename, "_average.png"), plot=., device="png", width=1920/150, height=1080/150, dpi=300)
    
    print("Cr plotting successful!")
    
  }
  
  unified = unified %>% mutate(ID = filename)
  ga_unified = rbind(ga_unified,unified)
  pupil_df  = rbind(pupil_df,pupil_vp %>% mutate(ID=filename))
  
  
  # stuff for calculating computing duration
  
  end.time <- Sys.time()  
  time.taken <- end.time - start.time
  code_times = c(code_times,round(as.numeric(time.taken),1))
  mean_code_time = mean(code_times)
  
  
  #if(match(subject_inmat,codes) == 1){
  #  print(paste0(round(match(subject_inmat,codes)/length(codes)*100,1), "% ..... ",
  #               match(subject_inmat,codes), " of ", length(codes), " files processed!"))}
  #else{
  #  print(paste0(round(match(subject_inmat,codes)/length(codes)*100,1), "% ..... ",
  #               match(subject_inmat,codes), " of ", length(codes), " files processed! ..... ", 
  #               round((mean_code_time*length(codes)-mean_code_time * match(subject_inmat,codes))/60,1), 
  #               " min remaining"))
  #}
  
} # end inmat for loop

# Read or Save ga_unified and pupil_df for further processing

# saveRDS(ga_unified,"ET_pupil_ga_unified.RData")
# saveRDS(pupil_df,"ET_pupil_df.RData")

ga_unified <- readRDS("ET_ga_unified.RData")
pupil_df <- readRDS("ET_pupil_df.Rdata")


ga_unified = ga_unified %>% 
  left_join(trigger_mat %>% select(subject, trial, outcome) %>% mutate(ID = as.integer(substr(subject, 5, 6))), by=c("ID", "trial")) %>% 
  mutate(outcome = ifelse(is.na(outcome), "no outcome", outcome))

pupil_df = pupil_df %>% 
  left_join(trigger_mat %>% select(subject, trial, outcome) %>% mutate(ID = as.integer(substr(subject, 5, 6))), by=c("ID", "trial")) %>% 
  mutate(outcome = ifelse(is.na(outcome), "no outcome", outcome))

# Signal quality check
mean(as.numeric(pupil_df$interpolated)) #30.4% interpolated sample points
sum(pupil_df$valid == T) #6291 valid trials
sum(pupil_df$valid == F) #2709 invalid trials
sum(pupil_df$valid == T) / (sum(pupil_df$valid == T) + sum(pupil_df$valid == F) ) #70% valid trials

pupil_df %>% group_by(condition) %>%
  summarise(valid_sp = mean(as.numeric(interpolated)), 
            valid_trials = sum(valid == T))  # no condition differences!

responder <- pupil_df %>% filter(condition %in% c(1,2,3)) %>% group_by(ID) %>% 
  summarise(dilators = mean(as.numeric(interpolated))) %>% filter(dilators < 0.5) %>% .$ID

non_responder <- pupil_df %>% filter(condition %in% c(1,2,3)) %>% group_by(ID) %>% 
  summarise(dilators = mean(as.numeric(interpolated))) %>% filter(dilators >= 0.5) %>% .$ID



# Plot Unified Data
# Acquisition
# Long format for statistical testing
pupil_df_long_acq <- pupil_df %>%
  # filter(outcome == "no outcome") %>%
  #filter(ID %in% responder) %>%
  pivot_longer(dilation_0:dilation_15, names_to = "timebin", values_to ="diameter") %>%
  separate(timebin,c("quark","timebin")) %>% mutate(timebin = as.numeric(timebin)) %>%
  # filter(timebin > 4 & timebin < 21) %>% #filter(valid == T) %>%
  filter(condition %in% c(2,3,4,5)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  select(ID, trial, condition, outcome, timebin, diameter) %>% 
  mutate(diameter = as.numeric(diameter[,1])) %>% 
  mutate(time = timebin * 0.5) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

main_effect_threat = list()
main_effect_social = list()
interaction_effect = list()
alpha = .05 / length(unique(pupil_df_long_acq$time))
for (timepoint in unique(pupil_df_long_acq$time)) {
  # timepoint = 3
  data = pupil_df_long_acq %>% 
    filter(time == timepoint) %>% 
    mutate(ID = as.factor(ID), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
    group_by(ID, condition_social, condition_threat)
  model = lmer(diameter ~ condition_social + condition_threat + condition_social:condition_threat + (1|ID), data)
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

ga_unified %>% filter(valid == TRUE) %>%
  # filter(outcome == "no outcome") %>%
  #filter(ID %in% responder) %>%
  .$diameter %>% bind_rows() %>%
  mutate(condition = as.factor(condition)) %>%
  filter(condition %in% c(2,3,4,5)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  group_by(condition,samplepoint) %>% 
  summarise(diameter.mean = mean(diameter), diameter.se = sd(diameter)/sqrt(length(diameter)), time = mean(time)) %>%
  {ggplot(., aes(x=time, y=diameter.mean)) +
      geom_vline(xintercept=0, colour="black",linetype="solid") + #zero
      geom_line(aes(colour=condition)) +
      geom_ribbon(aes(ymin=diameter.mean-diameter.se, ymax=diameter.mean+diameter.se, colour=condition, fill=condition), color = NA, alpha=.2) +
      geom_segment(data = main_effect_social, aes(x=times_start, xend = times_end, y=-70, yend=-70, size="Main Effect Social"), colour = "#e874ff", linewidth = 1, inherit.aes=FALSE) +
      geom_segment(data = main_effect_threat, aes(x=times_start, xend = times_end, y=-75, yend=-75, size="Main Effect Threat"), colour = "#ff8383", linewidth = 1, inherit.aes=FALSE) +
      geom_segment(data = interaction_effect, aes(x=times_start, xend = times_end, y=-80, yend=-80, size="Social x Threat Interaction "), colour = "#ffdd74", linewidth = 1, inherit.aes=FALSE) +
      scale_x_continuous("Time [s]",limits=c(-0.5, 8), minor_breaks=c(0,1,2,3,4,5,6,7,8), breaks=c(0, 2, 4, 6, 8)) +
      scale_y_continuous("Pupil Diameter",limits=c(-80, 150)) +
      scale_color_viridis_d(aesthetics = c("colour", "fill")) +
      theme_bw() + 
      # labs(title = paste("Pupil (N = ", n_distinct(pupil_df_long_acq$ID), ")", sep="")) +
      scale_size_manual("effects", values=rep(1,3), guide=guide_legend(override.aes = list(colour=c("#e874ff", "#ff8383", "#ffdd74")))) # 
  }
ggsave("../plots/Pupil/cs_acq.png",type="cairo-png", width=2500/400, height=1080/300, dpi=300)

# Test
# Long format for statistical testing
pupil_df_long_test <- pupil_df %>%
  #filter(ID %in% responder) %>%
  pivot_longer(dilation_0:dilation_19, names_to = "timebin", values_to ="diameter") %>%
  separate(timebin,c("quark","timebin")) %>% mutate(timebin = as.numeric(timebin)) %>%
  # filter(timebin > 4 & timebin < 21) %>% #filter(valid == T) %>%
  filter(condition %in% c(2,3,4,5)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  select(ID, trial, condition, timebin, diameter) %>% 
  mutate(diameter = as.numeric(diameter[,1])) %>% 
  mutate(time = timebin * 0.5) %>% 
  mutate(condition_social = if_else(str_detect(condition, "non-social"), "non-social", "social")) %>% 
  mutate(condition_threat = if_else(str_detect(condition, "pos"), "pos", "neg"))

main_effect_threat = list()
main_effect_social = list()
interaction_effect = list()
alpha = .05 / length(unique(pupil_df_long_test$time))
for (timepoint in unique(pupil_df_long_test$time)) {
  # timepoint = 0
  data = pupil_df_long_test %>% 
    filter(time == timepoint) %>% 
    mutate(ID = as.factor(ID), condition_social = as.factor(condition_social), condition_threat = as.factor(condition_threat)) %>%
    group_by(ID, condition_social, condition_threat)
  
  model = lmer(diameter ~ condition_social + condition_threat + condition_social:condition_threat + (1|ID), data)
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

ga_unified %>% filter(valid == TRUE) %>%
  #filter(ID %in% responder) %>%
  .$diameter %>% bind_rows() %>%
  mutate(condition = as.factor(condition)) %>% 
  # filter(trial < 60) %>%
  filter(condition %in% c(6,7,8,9)) %>%
  mutate(across('condition', str_replace_all, rep_str)) %>% 
  group_by(condition,samplepoint) %>% 
  summarise(diameter.mean = mean(diameter), diameter.se = sd(diameter)/sqrt(length(diameter)), time = mean(time)) %>%
  {ggplot(., aes(x=time, y=diameter.mean)) +
      geom_vline(xintercept=0, color="black",linetype="solid") + #zero = picture onset
      geom_vline(xintercept=10, color="black",linetype="solid") + #picture offset
      geom_line(aes(colour=condition)) +
      geom_ribbon(aes(ymin=diameter.mean-diameter.se, ymax=diameter.mean+diameter.se, colour=condition, fill=condition), color = NA, alpha=.2) +
      geom_segment(data = main_effect_social, aes(x=times_start, xend = times_end, y=-90, yend=-90, size="Main Effect Social"), colour = "#e874ff", linewidth = 1, inherit.aes=FALSE) +
      geom_segment(data = main_effect_threat, aes(x=times_start, xend = times_end, y=-93, yend=-93, size="Main Effect Threat"), colour = "#ff8383", linewidth = 1, inherit.aes=FALSE) +
      geom_segment(data = interaction_effect, aes(x=times_start, xend = times_end, y=-96, yend=-96, size="Social x Threat Interaction "), colour = "#ffdd74", linewidth = 1, inherit.aes=FALSE) +
      scale_x_continuous("Time [s]",limits=c(-0.5, 11), minor_breaks=c(0,1,2,3,4,5,6,7,8,9,10), breaks=c(0, 2, 4, 6, 8, 10)) +
      scale_y_continuous("Pupil Diameter", breaks=c(-80,-40, 0, 40), minor_breaks=c(-80, -60, -40, -20, 0, 20, 40, 60)) +
      scale_color_viridis_d(aesthetics = c("colour", "fill")) +
      theme_bw() +
      scale_size_manual("effects", values=rep(1,4), guide=guide_legend(override.aes = list(colour=c("#e874ff", "#ff8383", "#ffdd74"))))
  }
ggsave("../plots/Pupil/cs_test.png",type="cairo-png", width=2500/400, height=1080/300, dpi=300)
