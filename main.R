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
  draw_label("Experiment 1", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_scores_exp1 <- plot_grid(title_exp1, plot_scores_exp1, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# plots experiment 2
plot_scores_exp2 <- plot_grid(plot_spai_exp2, plot_stai_exp2,
                           nrow = 1, labels=c("C", "D"), align="vh")

# add title
title_exp2 <- ggdraw() + 
  draw_label("Experiment 2", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
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

plot_ratings <- plot_grid(plot_ratings_exp1 + theme(legend.position="none"), plot_ratings_exp2 + theme(legend.position="none"),
                          labels = c("A", "B"), nrow = 1, align = 'vh', rel_widths = c(3, 4))

# add title
title_exp1 <- ggdraw() + 
  draw_label("Experiment 1", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
title_exp2 <- ggdraw() + 
  draw_label("Experiment 2", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_titles <- plot_grid(title_exp1, title_exp2, nrow=1, rel_widths = c(3, 4))
plot_ratings <- plot_grid(plot_titles, plot_ratings, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# add legend
legend = get_legend(plot_ratings_exp1)
plot_ratings <- plot_grid(plot_ratings, legend, nrow = 1, rel_widths = c(7, 1.3))

ggsave(file.path("Plots", "ratings.svg"), width=10, height=4)

###############################################################################
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
                        "plot_lat_sacc_exp2", "plot_fix_test_exp2")))

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
plot_acq <- plot_grid(plot_saccades_acq_exp1 + theme(legend.position="none"), plot_saccades_acq_exp2 + theme(legend.position="none"),
                      labels = c("A", "B"), nrow = 1, align = 'vh')

# add title
title_exp1 <- ggdraw() + 
  draw_label("Experiment 1", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
title_exp2 <- ggdraw() + 
  draw_label("Experiment 2", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_titles <- plot_grid(title_exp1, title_exp2, nrow=1)
plot_acq <- plot_grid(plot_titles, plot_acq, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# add legend
legend = get_legend(plot_saccades_acq_exp1)
plot_acq <- plot_grid(plot_acq, legend, nrow = 1, rel_widths = c(8, 1))

ggsave(file.path("Plots", "gaze_acq.svg"), width=12, height=3.5)

# experiment 2
# familiar
plot_comp_fam <- plot_grid(plot_prop_sacc_familiar_exp2 + theme(legend.position="none"), plot_prop_first_sacc_familiar_exp2 + theme(legend.position="none"),
                      labels = c("A", "B"), nrow = 1, align = 'vh')

# add title
title_fam <- ggdraw() + 
  draw_label("Trials Without a Novel Stimulus", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot

plot_comp_fam <- plot_grid(title_fam, plot_comp_fam, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# novel
plot_comp_nov <- plot_grid(plot_prop_sacc_novel_exp2 + theme(legend.position="none"), plot_prop_first_sacc_nov_exp2 + theme(legend.position="none"),
                           labels = c("C", "D"), nrow = 1, align = 'vh')

# add title
title_nov <- ggdraw() + 
  draw_label("Trials With a Novel Stimulus", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot

plot_comp_nov <- plot_grid(title_nov, plot_comp_nov, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# combine plot
plot_comp <- plot_grid(plot_comp_fam, plot_comp_nov, nrow = 2)

# add legend
legend = get_legend(plot_prop_sacc_familiar_exp2)
plot_comp <- plot_grid(plot_comp, legend, nrow = 1, rel_widths = c(12, 1))

ggsave(file.path("Plots", "gaze_comp_exp2.svg"), width=12, height=7)


# Plots Supplements
# plots acquisition-phase
# plots experiment 1
plot_acq_exp1 <- plot_grid(plot_latencies_acq_exp1 + theme(legend.position="none"), plot_lengths_acq_exp1 + theme(legend.position="none"), plot_fixations_acq_exp1 + theme(legend.position="none"),
                           labels=c("A", "B", "C"), nrow = 1, align="vh")

# add title
title_exp1 <- ggdraw() + 
  draw_label("Experiment 1", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_acq_exp1 <- plot_grid(title_exp1, plot_acq_exp1, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# plots experiment 2
plot_acq_exp2 <- plot_grid(plot_latencies_acq_exp2 + theme(legend.position="none"), plot_lengths_acq_exp2 + theme(legend.position="none"), plot_fixations_acq_exp2 + theme(legend.position="none"),
                           nrow = 1, labels=c("D", "E", "F"), align="vh")

# add title
title_exp2 <- ggdraw() + 
  draw_label("Experiment 2", fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 6)) # add margin on the left of the drawing canvas, so title is aligned with left edge of first plot
plot_acq_exp2 <- plot_grid(title_exp2, plot_acq_exp2, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins

# combine plot
plot_acq <- plot_grid(plot_acq_exp1, plot_acq_exp2, nrow = 2)

# add legend
legend = get_legend(plot_saccades_acq_exp1)
plot_acq <- plot_grid(plot_acq, legend, nrow = 1, rel_widths = c(11, 1))

ggsave(file.path("Plots", "gaze_acq_supplm.svg"), width=12, height=7)

# plots test-phase experiment 1
plot_test_exp1 <- plot_grid(plot_sacc_test_exp1 + theme(legend.position="none"), plot_fixations_test_exp1 + theme(legend.position="none"),
                      labels = c("A", "B"), nrow = 1, align = 'vh')

# add legend
legend = get_legend(plot_fixations_test_exp1)
plot_test_exp1 <- plot_grid(plot_test_exp1, legend, nrow = 1, rel_widths = c(8, 1))

ggsave(file.path("Plots", "gaze_test_exp1_supplm.svg"), width=12, height=3.5)

# plots competition-phase experiment 2
plot_comp_exp2 <- plot_grid(plot_lat_sacc_exp2 + theme(legend.position="none"), plot_fix_test_exp2 + theme(legend.position="none"),
                            labels = c("A", "B"), nrow = 1, align = 'vh')

# add legend
legend = get_legend(plot_lat_sacc_exp2)
plot_comp_exp2 <- plot_grid(plot_comp_exp2, legend, nrow = 1, rel_widths = c(8, 1))

ggsave(file.path("Plots", "gaze_comp_exp2_supplm.svg"), width=12, height=3.5)

