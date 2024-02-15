###############################################################################
# Gaze Contingent Avoidance Project
# Sabrina Gado & Yannik Stegmann
###############################################################################

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
ucr_plots = T
cr_plots = T


sample_rate = 50  #samplerate after export
trial_length = 12  #length of trial in s

downsampling = F
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
  
  transfer_saccade <- function(blink, saccade) {
    blink_new <- logical(length(blink))
    
    rle_saccade <- rle(saccade)
    rle_saccade$indices <- c(1, cumsum(rle_saccade$lengths) + 1)[1:length(rle_saccade$values)]
    saccade_indices <- rle_saccade$indices[rle_saccade$values]
    saccade_lengths <- rle_saccade$lengths[rle_saccade$values]
    
    for (index in 1:length(saccade_indices)) {
      # index = 7
      start <- saccade_indices[index]
      end <- saccade_indices[index] + saccade_lengths[index] - 1
      if(any(blink[start:end])) {
        blink_new[start:end] <- TRUE
      }
    }
    
    return(blink_new)
  }
}


data = read.csv2("../Physio/pupil.txt", sep ="\t", dec = ",", na.strings=".")

data <- data %>%
  mutate(diameter = ifelse(is.na(RIGHT_PUPIL_SIZE), LEFT_PUPIL_SIZE, RIGHT_PUPIL_SIZE)) %>% 
  select(-LEFT_PUPIL_SIZE, -RIGHT_PUPIL_SIZE)

data <- data %>%
  group_by(RECORDING_SESSION_LABEL, TRIAL_INDEX) %>%
  mutate(timestamp = (TIMESTAMP - first(TIMESTAMP))) %>% 
  ungroup() %>% 
  mutate(saccade = as.logical(RIGHT_IN_SACCADE) | as.logical(LEFT_IN_SACCADE),
         blink = as.logical(RIGHT_IN_BLINK) | as.logical(LEFT_IN_BLINK)) %>% 
  mutate(saccade = replace_na(saccade, FALSE),
         blink = replace_na(blink, FALSE))

data$blink <- transfer_saccade(data$blink, data$saccade)
data <- data %>% select(-c(RIGHT_IN_SACCADE, LEFT_IN_SACCADE, RIGHT_IN_BLINK, LEFT_IN_BLINK, TIMESTAMP))

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
  filter(grepl("Onset", event))


for (subject_inmat in codes){ 
  
  # subject_inmat = codes[[1]]
  start.time <- Sys.time()
  
  vp = ifelse(is.na(as.numeric(substr(subject_inmat,20,21))),as.numeric(substr(subject_inmat,20,20)),as.numeric(substr(subject_inmat,20,21)))
  
  pupil = data %>% filter(RECORDING_SESSION_LABEL == subject_inmat) %>%
    rename(ID = RECORDING_SESSION_LABEL, trial = TRIAL_INDEX, time = timestamp) %>%
    mutate(ID = vp, 
           sample = 1:n(),
           time = time / 1000)
  
  conditions = trigger_mat %>% filter(ID == unique(pupil$ID)) %>%
    select(ID, trial, condition, trigger, outcome)
  
  messages_vp = messages %>% 
    filter(subject == vp) %>%
    rename(ID = subject, event_time = time) %>% 
    left_join(. ,conditions, by=c("ID","trial")) %>% 
    mutate(trigger = ifelse(grepl("Feedback", event), 100, trigger)) %>% 
    select(trial, event_time, trigger) %>% 
    mutate(event_time = event_time / 1000)
  
  pupil = pupil %>% 
    left_join(messages_vp, join_by(trial, closest(time >= event_time))) %>% 
    mutate(trigger = ifelse((is.na(lag(trigger, 1))| (trigger != lag(trigger, 1))) & trigger == lead(trigger, 1), trigger, NA)) %>%
    mutate(trigger = ifelse(is.na(trigger), 0, trigger)) %>% 
    select(-event_time)
  
  filename = unique(pupil$ID)
  
  print(paste0(vp, " data reading successful!"))
  
  # interpolation of missing values
  pupil = pupil %>% 
    mutate(diameter = ifelse(blink, NA, diameter))
  
  pupil <- pupil %>% 
    group_by(trial) %>% 
    mutate(diameter = ifelse((diameter > mean(diameter, na.rm = T) - 3 * sd(diameter, na.rm = T)) & (diameter < mean(diameter, na.rm = T) + 3 * sd(diameter, na.rm = T)), diameter, NA)) %>% 
    ungroup()
  
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
  
  pupil_vp = pupil %>%  select(-diameter) %>% 
    filter(trigger != 0) %>% ungroup() %>% #One trial per row per VP
    mutate(timepoint = ifelse(trigger < 100, "start", "end")) %>% 
    pivot_wider(names_from=timepoint, values_from=time) %>% 
    mutate(end = lead(end)) %>% 
    filter(trigger < 100) %>% 
    mutate(trial = 1:n(),
           condition = trigger,
           time_start = start,
           time_end = ifelse(is.na(end), time_start + trial_length, end),
           sample_start = sample, 
           sample_end =  {sample+(trial_length*ifelse(downsampling,sample_rate_new,sample_rate))} %>% round()) %>% 
    select(-c(trigger, start, end, sample, interpolated, saccade, blink)) %>% select(ID, trial, condition, everything())
  
  evaluate_interpolation <- function(trial_identifier, start, end) {
    # Subset data based on start and end times
    subset_data <- pupil %>%
      filter(trial == trial_identifier, time >= start - abs(min(baselineWindow)), time < end)
    
    # Calculate mean of the 'interpolated' column
    mean_value <- mean(subset_data$interpolated)
    
    # Round the mean to three decimal places
    rounded_mean <- round(mean_value, 3)
    
    return(rounded_mean)
  }
  
 
  pupil_vp$interpolated <- numeric(nrow(pupil_vp))
  for (t in 1:nrow(pupil_vp)){ #count interpolated samples per trial
    pupil_vp$interpolated[[t]] <- evaluate_interpolation(trial_identifier=pupil_vp$trial[[t]], start= pupil_vp$time_start[[t]], end=pupil_vp$time_end[[t]])
    
    exclude = rep(TRUE,ifelse(downsampling, sample_rate_new/2, sample_rate/2)) #exclude if more than (sampling rate / 2) interpolated samplepoints (0.5s) in a row
    }
  
  pupil_vp = pupil_vp %>% 
    mutate(valid = ifelse(interpolated < .25, TRUE, FALSE))  #exclude if more than 25% interpolated data
  
  pupil <- pupil %>% 
    left_join(pupil_vp %>% select(c(ID, trial, time_end, time_start)), by = c("ID", "trial"))
  
  pupil <- pupil %>% 
    filter(time >= time_start - abs(min(baselineWindow)),
           time <= time_end) %>% 
    mutate(time = time - time_start) %>% select(-c(time_start, time_end, saccade, blink, trigger))
  
  pupils_df_list[[filename]] = pupil_vp
  pupils_list[[filename]] = pupil
  pupil_vp <- pupil_vp %>% mutate(trial_idx = trial)
  
  unified = pupil_vp %>% mutate(diameter = trial_idx %>% lapply(function(trial_idx) 
    pupil %>% filter(trial == trial_idx)))
  
  
  # for (t in 1:nrow(unified)) {
  #   unified$diameter[[t]] = unified$diameter[[t]] %>% 
  #     mutate(trial = t, 
  #            #time = time - (min(time)+abs(min(baselineWindow))), 
  #            #samplepoint = sample-min(sample),
  #            condition = unified$condition[[t]], 
  #            #diameter = diameter - unified$baseline[[t]] #unify starting time to allow overlap
  #     )}
  # 
  
  # baseline and level scoring  
  test.time <- Sys.time()  
  pupil_diameter_binded <- unified$diameter %>% bind_rows()
  
  pupil_diameter_vp <- tibble()
  
  for (t in 1:nrow(pupil_vp)) {
    # t = 1
    start_time = pupil_vp$time_start[t]
    start_sample = pupil_vp$sample_start[t]
    trial_idx = pupil_vp$trial[t]
    
    pupil_diameter_binded_trial = pupil_diameter_binded %>% filter(trial == t)
    
    baseline_trial = pupil_diameter_binded_trial %>% filter(time >= min(baselineWindow), time <= 0) %>% pull(diameter) %>% mean(na.rm=T)
    
    diameter_level_trial = pupil_diameter_binded_trial %>% filter(time >= 0) %>% pull(diameter) %>% mean(na.rm=T)
    
    diameter_level_trial_0 = pupil_diameter_binded_trial %>% filter(time >= 0, time <= 0.5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_1 = pupil_diameter_binded_trial %>% filter(time >= 0.5, time <= 1) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_2 = pupil_diameter_binded_trial %>% filter(time >= 1, time <= 1.5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_3 = pupil_diameter_binded_trial %>% filter(time >= 1.5,time <= 2) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_4 = pupil_diameter_binded_trial %>% filter(time >= 2, time <= 2.5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_5 = pupil_diameter_binded_trial %>% filter(time >= 2.5,time <= 3) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_6 = pupil_diameter_binded_trial %>% filter(time >= 3, time <= 3.5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_7 = pupil_diameter_binded_trial %>% filter(time >= 3.5,time <= 4) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_8 = pupil_diameter_binded_trial %>% filter(time >= 4, time <= 4.5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_9 = pupil_diameter_binded_trial %>% filter(time >= 4.5, time <= 5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_10 = pupil_diameter_binded_trial %>% filter(time >= 5,time <= 5.5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_11 = pupil_diameter_binded_trial %>% filter(time >= 5.5, time <= 6) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_12 = pupil_diameter_binded_trial %>% filter(time >= 6,time <= 6.5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_13 = pupil_diameter_binded_trial %>% filter(time >= 6.5, time <= 7) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_14 = pupil_diameter_binded_trial %>% filter(time >= 7,time <= 7.5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_15 = pupil_diameter_binded_trial %>% filter(time >= 7.5, time <= 8) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_16 = pupil_diameter_binded_trial %>% filter(time >= 8,time <= 8.5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_17 = pupil_diameter_binded_trial %>% filter(time >= 8.5, time <= 9) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_18 = pupil_diameter_binded_trial %>% filter(time >= 9,time <= 9.5) %>% pull(diameter) %>% mean(na.rm=T)
    diameter_level_trial_19 = pupil_diameter_binded_trial %>% filter(time >= 9.5, time <= 10) %>% pull(diameter) %>% mean(na.rm=T)
    
    dilations = tibble(dilation = diameter_level_trial - baseline_trial,
                       dilation_0 = diameter_level_trial_0 - baseline_trial,
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
                       dilation_16 = diameter_level_trial_16 - baseline_trial,
                       dilation_17 = diameter_level_trial_17 - baseline_trial,
                       dilation_18 = diameter_level_trial_18 - baseline_trial,
                       dilation_19 = diameter_level_trial_19 - baseline_trial,
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
    # t = 1
    unified$diameter[[t]] = unified$diameter[[t]] %>% 
      mutate(subject = vp,
             condition = pupil_vp$condition[[t]],
             time = time - (min(time)+abs(min(baselineWindow))), 
             samplepoint = sample-min(sample),
             diameter = diameter - pupil_vp$baseline[[t]])
    }
  
  if (cr_plots) {
    
    unified %>% filter(valid) %>% 
      .$diameter %>% bind_rows() %>% 
      mutate(trial = as.factor(trial)) %>% 
      mutate(across('condition', str_replace_all, rep_str)) %>% 
      {ggplot(., aes(x=time, y=diameter, color=trial)) + facet_wrap(vars(condition)) +
          geom_vline(xintercept=0, color="black",linetype="solid") + #zero 
          geom_path() + scale_color_viridis_d() +
          xlab("Time") + ylab("Diameter") +
          ggtitle(filename) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))} %>% 
      #print()
      ggsave(paste0("../plots/Pupil/subject_level/", filename, ".png"), plot=., device="png", width=1920/150, height=1080/150, dpi=300)
    
    unified %>% filter(valid) %>% 
      .$diameter %>% bind_rows() %>%
      mutate(across('condition', str_replace_all, rep_str)) %>% 
      mutate(condition = as.factor(condition)) %>% 
      group_by(condition, samplepoint) %>% 
      summarise(diameter = mean(diameter), time = mean(time)) %>%
      {ggplot(., aes(x=time, y=diameter, color=condition, group=condition,linetype=condition)) +
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

saveRDS(ga_unified,"ET_pupil_ga_unified.RData")
saveRDS(pupil_df,"ET_pupil_df.RData")

pupil_df <- pupil_df %>% 
  mutate(valid = as.logical(valid))

pupil_df_invalid <- pupil_df %>% filter(condition %in% c(6,7,8,9)) %>% summarise(invalid = mean(!valid), .by=c(ID)) %>% arrange(desc(invalid))
pupil_df_invalid %>% summarise(mean_invalid = mean(invalid)*100, sd_invalid = sd(invalid)*100)
