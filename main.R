###############################################################################
# Demographics and Questionnaires
###############################################################################
rm(list=(ls()))

source(file.path("Study 1", "Analyses - Study 1", "GCA - Scores.R"))
source(file.path("Study 2", "Analyses - Study 2", "GCA2 - Scores.R"))

rm(list=setdiff(ls(), c("plot_spai_exp1", "plot_spai_exp2", "plot_stai_exp1", "plot_stai_exp2")))

# Plots Supplements
# plots experiment 1
plot_scores_exp1 <- plot_grid(plot_spai_exp1, plot_stai_exp1,
                              labels=c("A", "B"), nrow = 1, align="vh")

# add title
title_exp1 <- ggdraw() + 
  draw_label("Experiment 1", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_scores_exp1 <- plot_grid(title_exp1, plot_scores_exp1, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# plots experiment 2
plot_scores_exp2 <- plot_grid(plot_spai_exp2, plot_stai_exp2,
                           nrow = 1, labels=c("C", "D"), align="vh")

# add title
title_exp2 <- ggdraw() + 
  draw_label("Experiment 2", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6), text = element_text(size = 26)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_scores_exp2 <- plot_grid(title_exp2, plot_scores_exp2, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# combine plot
plot_scores <- plot_grid(plot_scores_exp1, plot_scores_exp2, nrow = 2)

ggsave(file.path("Plots", "scores.svg"), width=12, height=7)

###############################################################################
# Discrimination Task
###############################################################################
rm(list=(ls()))
source(file.path("Study 1", "Analyses - Study 1", "GCA - Discrimination task.R"))
source(file.path("Study 2", "Analyses - Study 2", "GCA2 - Discrimination task.R"))

###############################################################################
# Ratings
###############################################################################
rm(list=(ls()))
source(file.path("Study 1", "Analyses - Study 1", "GCA - Ratings.R"))
source(file.path("Study 2", "Analyses - Study 2", "GCA2 - Ratings.R"))

rm(list=setdiff(ls(), c("plot_ratings_exp1", "plot_ratings_exp2")))

plot_ratings <- plot_grid(plot_ratings_exp1 + theme(legend.position="none") + scale_y_continuous(breaks=c(2, 4, 6, 8, 10), limits=c(1, 10)),
                          plot_ratings_exp2 + theme(legend.position="none") + scale_y_continuous(breaks=c(2, 4, 6, 8, 10), limits=c(1, 10)),
                          labels = c("A", "B"), nrow = 1, align = 'vh', rel_widths = c(3, 4))

# add title
title_exp1 <- ggdraw() + 
  draw_label("Experiment 1", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6), text = element_text(size = 14)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
title_exp2 <- ggdraw() + 
  draw_label("Experiment 2", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6), text = element_text(size = 14)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_titles <- plot_grid(title_exp1, title_exp2, nrow=1, rel_widths = c(3, 4))
plot_ratings <- plot_grid(plot_titles, plot_ratings, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# add legend
legend <- get_legend(plot_ratings_exp1)
plot_ratings <- plot_grid(plot_ratings, legend, nrow = 1, rel_widths = c(7, 1.2))

ggsave(file.path("Plots", "ratings.svg"), width=12, height=5)

################################################################################
# Gaze
###############################################################################
rm(list=(ls()))
source(file.path("Study 1", "Analyses - Study 1", "GCA - ET.R"))
source(file.path("Study 2", "Analyses - Study 2", "GCA2 - ET.R"))

rm(list=setdiff(ls(), c("plot_saccades_acq_exp1", "plot_latencies_acq_exp1", "plot_lengths_acq_exp1", "plot_fixations_acq_exp1",
                        "plot_fixations_test_exp1", "plot_sacc_test_exp1", 
                        "plot_saccades_acq_exp2", "plot_latencies_acq_exp2", "plot_lengths_acq_exp2", "plot_fixations_acq_exp2",
                        "plot_prop_sacc_familiar_exp2", "plot_prop_sacc_novel_exp2", 
                        "plot_prop_first_sacc_familiar_exp2", "plot_prop_first_sacc_nov_exp2",
                        "plot_lat_sacc_familiar_exp2", "plot_lat_sacc_novel_exp2",
                        "plot_fix_test_familiar_exp2", "plot_fix_test_novel_exp2",
                        "saccades.acq.prop1", "saccades.acq.prop2")))

# Set ylim
min_y_sacc = min(c(layer_scales(plot_saccades_acq_exp1)$y$range$range, layer_scales(plot_saccades_acq_exp2)$y$range$range))
max_y_sacc = max(c(layer_scales(plot_saccades_acq_exp1)$y$range$range, layer_scales(plot_saccades_acq_exp2)$y$range$range))
plot_saccades_acq_exp1 <- plot_saccades_acq_exp1 + ylim(min_y_sacc, max_y_sacc)
plot_saccades_acq_exp2 <- plot_saccades_acq_exp2 + ylim(min_y_sacc, max_y_sacc)

min_y_lat = min(c(layer_scales(plot_latencies_acq_exp1)$y$range$range, layer_scales(plot_latencies_acq_exp2)$y$range$range))
max_y_lat = max(c(layer_scales(plot_latencies_acq_exp1)$y$range$range, layer_scales(plot_latencies_acq_exp2)$y$range$range))
plot_latencies_acq_exp1 <- plot_latencies_acq_exp1 + ylim(min_y_lat, max_y_lat)
plot_latencies_acq_exp2 <- plot_latencies_acq_exp2 + ylim(min_y_lat, max_y_lat)

min_y_len = min(c(layer_scales(plot_lengths_acq_exp1)$y$range$range, layer_scales(plot_lengths_acq_exp2)$y$range$range))
max_y_len = max(c(layer_scales(plot_lengths_acq_exp1)$y$range$range, layer_scales(plot_lengths_acq_exp2)$y$range$range))
plot_lengths_acq_exp1 <- plot_lengths_acq_exp1 + ylim(min_y_len, max_y_len)
plot_lengths_acq_exp2 <- plot_lengths_acq_exp2 + ylim(min_y_len, max_y_len)

min_y_fix = min(c(layer_scales(plot_fixations_acq_exp1)$y$range$range, layer_scales(plot_fixations_acq_exp2)$y$range$range))
max_y_fix = max(c(layer_scales(plot_fixations_acq_exp1)$y$range$range, layer_scales(plot_fixations_acq_exp2)$y$range$range))
plot_fixations_acq_exp1 <- plot_fixations_acq_exp1 + ylim(min_y_fix, max_y_fix)
plot_fixations_acq_exp2 <- plot_fixations_acq_exp2 + ylim(min_y_fix, max_y_fix)

### Create subplots
# Plots Paper
# experiment 1
plot_acq_exp1 <- plot_grid(plot_saccades_acq_exp1 + theme(legend.position="none") + scale_y_continuous(breaks=c(0, 0.5, 1)), plot_latencies_acq_exp1 + theme(legend.position="none"),
                           labels=c("A", "B"), nrow = 1, align="vh")

# add title
title_exp1 <- ggdraw() + 
  draw_label("Experiment 1", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_acq_exp1 <- plot_grid(title_exp1, plot_acq_exp1, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# plots experiment 2
plot_acq_exp2 <- plot_grid(plot_saccades_acq_exp2 + theme(legend.position="none") + scale_y_continuous(breaks=c(0, 0.5, 1)), plot_latencies_acq_exp2 + theme(legend.position="none"),
                           nrow = 1, labels=c("C", "D"), align="vh")

# add title
title_exp2 <- ggdraw() + 
  draw_label("Experiment 2", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_acq_exp2 <- plot_grid(title_exp2, plot_acq_exp2, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# combine plot
plot_acq <- plot_grid(plot_acq_exp1, plot_acq_exp2, nrow = 2)

# add legend
legend <- get_legend(plot_saccades_acq_exp1)
plot_acq <- plot_grid(plot_acq, legend, nrow = 1, rel_widths = c(11, 1))

ggsave(file.path("Plots", "gaze_acq.svg"), width=12, height=8)


# Experiment 2
# familiar
plot_comp_fam <- plot_grid(plot_prop_sacc_familiar_exp2 + theme(legend.position="none"), plot_lat_sacc_familiar_exp2 + theme(legend.position="none"), plot_fix_test_familiar_exp2 + theme(legend.position="none"),
                      labels = c("A", "B", "C"), nrow = 1, align = 'vh')

# add title
title_fam <- ggdraw() + 
  draw_label("Trials Without a Novel Stimulus", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot

plot_comp_fam <- plot_grid(title_fam, plot_comp_fam, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# novel
plot_comp_nov <- plot_grid(plot_prop_sacc_novel_exp2 + theme(legend.position="none"), plot_lat_sacc_novel_exp2 + theme(legend.position="none"), plot_fix_test_novel_exp2 + theme(legend.position="none"),
                           labels = c("D", "E", "F"), nrow = 1, align = 'vh')

# add title
title_nov <- ggdraw() + 
  draw_label("Trials With a Novel Stimulus", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot

plot_comp_nov <- plot_grid(title_nov, plot_comp_nov, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# combine plot
plot_comp <- plot_grid(plot_comp_fam, plot_comp_nov, nrow = 2)

# add legend
legend <- get_legend(plot_prop_sacc_familiar_exp2)
plot_comp <- plot_grid(plot_comp, legend, nrow = 1, rel_widths = c(14, 1))

ggsave(file.path("Plots", "gaze_comp_exp2.svg"), width=14, height=8)


# Plots Supplements
# plots acquisition-phase
# plots experiment 1
plot_acq_exp1 <- plot_grid(plot_lengths_acq_exp1 + theme(legend.position="none"), plot_fixations_acq_exp1 + theme(legend.position="none"),
                           labels=c("A", "B"), nrow = 1, align="vh")

# add title
title_exp1 <- ggdraw() + 
  draw_label("Experiment 1", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_acq_exp1 <- plot_grid(title_exp1, plot_acq_exp1, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# plots experiment 2
plot_acq_exp2 <- plot_grid(plot_lengths_acq_exp2 + theme(legend.position="none"), plot_fixations_acq_exp2 + theme(legend.position="none"),
                           nrow = 1, labels=c("C", "D"), align="vh")

# add title
title_exp2 <- ggdraw() + 
  draw_label("Experiment 2", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_acq_exp2 <- plot_grid(title_exp2, plot_acq_exp2, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# combine plot
plot_acq <- plot_grid(plot_acq_exp1, plot_acq_exp2, nrow = 2)

# add legend
legend <- get_legend(plot_saccades_acq_exp1)
plot_acq <- plot_grid(plot_acq, legend, nrow = 1, rel_widths = c(11, 1))

ggsave(file.path("Plots", "gaze_acq_supplm.svg"), width=12, height=8)

# plots test-phase experiment 1
plot_test_exp1 <- plot_grid(plot_sacc_test_exp1 + theme(legend.position="none"), plot_fixations_test_exp1 + theme(legend.position="none"),
                      labels = c("A", "B"), nrow = 1, align = 'vh')

# add legend
legend <- get_legend(plot_fixations_test_exp1)
plot_test_exp1 <- plot_grid(plot_test_exp1, legend, nrow = 1, rel_widths = c(8, 1))

ggsave(file.path("Plots", "gaze_test_exp1_supplm.svg"), width=12, height=3.5)

# plots competition-phase experiment 2
# familiar
plot_comp_fam <- plot_grid(plot_prop_first_sacc_familiar_exp2 + theme(legend.position="none"),
                           labels = c("A"), nrow = 1, align = 'vh')
plot_comp_nov <- plot_grid(plot_prop_first_sacc_nov_exp2 + theme(legend.position="none"),
                           labels = c("B"), nrow = 1, align = 'vh')

plot_comp <- plot_grid(plot_comp_fam, plot_comp_nov, nrow = 1) # rel_heights values control vertical title margins

# add title
title_fam <- ggdraw() + 
  draw_label("Trials Without a Novel Stimulus", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
title_nov <- ggdraw() + 
  draw_label("Trials With a Novel Stimulus", fontface = 'bold', x = 0, hjust = 0, size = 18) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot

title <- plot_grid(title_fam, title_nov, nrow = 1) # rel_heights values control vertical title margins

plot_comp <- plot_grid(title, plot_comp, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# add legend
legend <- get_legend(plot_fix_test_familiar_exp2)
plot_comp <- plot_grid(plot_comp, legend, nrow = 1, rel_widths = c(12, 1))

ggsave(file.path("Plots", "gaze_comp_exp2_supplm.svg"), width=12, height=4)

# # Cross-Task Comparison regarding aversive US
# saccades.acq.prop1.av <- saccades.acq.prop1 %>% 
#   filter(condition_threat == "neg") %>% 
#   summarise(rel_freq = mean(relative_frequency_ROI), .by=c(subject))
# 
# saccades.acq.prop2.av <- saccades.acq.prop2 %>% 
#   filter(condition_threat == "neg") %>% 
#   summarise(rel_freq = mean(relative_frequency_ROI), .by=c(subject))
# 
# apa::t_test(saccades.acq.prop1.av$rel_freq, saccades.acq.prop2.av$rel_freq, alternative = "two.sided", paired=FALSE) %>% apa::t_apa()
# 
# saccades.acq.prop1.av <- saccades.acq.prop1.av %>% mutate(dataset = "Experiment 1") %>% select(c(rel_freq, dataset))
# saccades.acq.prop2.av <- saccades.acq.prop2.av %>% mutate(dataset = "Experiment 2") %>% select(c(rel_freq, dataset))
# saccades.acq.prop.av <- rbind(saccades.acq.prop1.av, saccades.acq.prop2.av)
# 
# saccades.acq.prop.av.summarise <- saccades.acq.prop.av %>% 
#   summarise(mean_rel_freq = mean(rel_freq), se_rel_freq = sd(rel_freq)/sqrt(n()), .by=c(dataset))
# 
# plot_US_av <- ggplot(saccades.acq.prop.av.summarise, aes(x = dataset, y = mean_rel_freq, group = dataset, color = dataset, fill=dataset)) +
#   geom_point(data = saccades.acq.prop.av, aes(x = dataset, y = rel_freq, group = dataset, color = dataset), 
#              position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.005, dodge.width = 0.6), alpha = 0.4, size=2, shape = 21) +
#   geom_point(position = position_dodge(width = 0.6), shape = "square", size=3) +
#   geom_line(position = position_dodge(width = 0.6), linewidth = 0.5) +
#   geom_errorbar(aes(ymin = mean_rel_freq - se_rel_freq, ymax = mean_rel_freq + se_rel_freq), position = position_dodge(width = 0.6), width = 0.25) +
#   labs(title = paste("Proportion of Trials with Saccades on CSneg", sep=""), x = NULL, y = "Proportion") +
#   # theme_minimal() +
#   theme(legend.position="none") + 
#   scale_fill_viridis_d("Condition", end = 0.15, begin = 0.85) + 
#   scale_color_viridis_d("Condition", end = 0.15, begin = 0.85)

##############################################################################
# Physiology
###############################################################################
rm(list=(ls()))
source(file.path("Study 1", "Analyses - Study 1", "GCA - Tests Physiology.R"))

rm(list=setdiff(ls(), c("plot_pupil_cluster", "plot_eda_cluster", "plot_hr_cluster")))

plot_physio <- plot_grid(plot_pupil_cluster + theme(legend.position="none"), plot_eda_cluster + theme(legend.position="none"), plot_hr_cluster + theme(legend.position="none"),
                         ncol = 3, labels=c("A", "B", "C"), align="vh")

# add legend
legend <- get_legend(plot_pupil_cluster)
plot_grid(plot_physio, legend, nrow = 1, reL_widths = c(20, 1))

ggsave(file.path("Plots", "physio.svg"), width=14, height=4)

ggsave(file.path("Plots", "pupil.svg"), plot_pupil_cluster, width=6, height=4)


##############################################################################
# Distance on Screen to Visual Ancle
###############################################################################

# Define monitor and viewing parameters
monitor_width_cm <- 53          # Monitor width in cm
resolution_width_px <- 1920     # Screen resolution width in pixels
viewing_distance_cm <- 50       # Viewing distance in cm

# Stimulus positions in pixels
x1 <- 1920/2 - (1920/2 * 0.7)
y1 <- 1080/2
x2 <- 1920/2
y2 <- 1080/2

# Step 1: Calculate pixel distance
pixel_distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)

# Step 2: Convert pixel distance to physical distance
physical_distance_cm <- pixel_distance * (monitor_width_cm / resolution_width_px)

# Step 3: Calculate the visual angle
visual_angle_rad <- 2 * atan(physical_distance_cm / (2 * viewing_distance_cm))
visual_angle_deg <- visual_angle_rad * (180 / pi)

cat("Visual angle (degrees):", visual_angle_deg, "Visual angle (rad):", visual_angle_rad, "\n")

# Stimulus positions in pixels
x1 <- 1920/2 - (1920/2 * 0.7)
y1 <- 1080/2
x2 <- 1920/2 + (1920/2 * 0.7)
y2 <- 1080/2

# Step 1: Calculate pixel distance
pixel_distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)

# Step 2: Convert pixel distance to physical distance
physical_distance_cm <- pixel_distance * (monitor_width_cm / resolution_width_px)

# Step 3: Calculate the visual angle
visual_angle_rad <- 2 * atan(physical_distance_cm / (2 * viewing_distance_cm))
visual_angle_deg <- visual_angle_rad * (180 / pi)

cat("Visual angle (degrees):", visual_angle_deg, "Visual angle (rad):", visual_angle_rad, "\n")

