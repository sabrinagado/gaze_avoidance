library(tidyverse)
library(mc2d) # for pempiricalD()

# For parallelization
library(parallel)
library(doSNOW)

load("rdata/labexpvr_physio_processed.RData")

# Functions --------------------------------------------------------------------

## General functions -----------------------------------------------------------

# Find clusters in the vector of t-statistics. `x` must be a call to
# pt_ttest_statistic_between() or pt_ttest_statistic_within() and `critical_t`
# must be a call to pt_critical_t_between() or pt_critical_t_between(),
# respectively.
# 
# Returns a tibble with the following columns:
# `start`  index of the start of the cluster(s)
# `end`    index of the end of the cluster(s)
# `length` length of the cluster(s)
# (all indices with respect to the series of time points)
clusters_over_time <- function(x, critical_t)
{
  # Find clusters in the t-statistic vector that are above/below the threshold
  x <- rle(abs(x$value) > critical_t)
  
  # Extract only clusters that are above the threshold
  cluster_lengths <- x$lengths[x$values]
  
  # Find the start and end of clusters via the cumsum of cluster lengths. For
  # cluster start, we'll put a 0/FALSE in front, so that we don't run into
  # indexing vec[0]
  cluster_start <- c(0L, cumsum(x$lengths))[which(c(FALSE, x$values)) - 1] + 1L
  cluster_end <- cumsum(x$lengths)[x$values]
  
  tibble(cluster = seq_along(cluster_lengths), start = cluster_start,
         end = cluster_end, length = cluster_lengths)
}

# clusters <- clusters_over_time(t, pt_critical_t_between(tmp, id = "id"))

## Functions for two-group between-subjects comparisons ------------------------

# tmp <-
#   eda |>
#   filter(phase == "habituation")

# This is an implementation of a cluster-based permutation test for a two-group
# between-subject comparison of 2d data (signal x time)
# 
# The data needs to be in a tibble with the following columns and must not
# include any NAs:
#
# (1) A column containing the dependent variable (dv) (e.g., heart rate, skin
#     conductance)
# (2) A column indicating the time points / samples (time)
# (3) A column containing the between group identifier (between), must be a
#     factor with two levels
# (4) A column containing the participant identifier (id)
#
# All functions expect column names as character.

# Function for calculating the t-score for a between-subject t-test. Used
# instead of `t.test(...)$statistic` for increased performance
t_statistic_between <- function (x, y) 
{
  nx <- length(x)
  mx <- mean(x)
  vx <- var(x)
  
  ny <- length(y)
  my <- mean(y)
  vy <- var(y)
  
  df <- nx + ny - 2
  
  v <- 0
  if (nx > 1) 
    v <- v + (nx - 1) * vx
  if (ny > 1) 
    v <- v + (ny - 1) * vy
  v <- v / df
  stderr <- sqrt(v * (1 / nx + 1 / ny))
  
  (mx - my) / stderr
}

# Calculate the t-statistic (between-subject comparison) for each time point
pt_ttest_statistic_between <- function(data, dv, between, time)
{
  data |>
    # For each time point ...
    group_by(!!sym(time)) |>
    # ... calculate the test statistic of a two-sample t-test
    summarize(value = t_statistic_between(
      # From the current subset of the tibble (i.e., per sample), extract the
      # dependent variable data of the first ...
      .data[[dv]][.data[[between]] == levels(.data[[between]])[1]],
      # ... and second level of the between-subjects factor
      .data[[dv]][.data[[between]] == levels(.data[[between]])[2]])) |>
    ungroup()
}

# Calculate the critical t-value for a two-sided test for the given sample size
pt_critical_t_between <- function(data, id, alpha = .05)
{
  qt(1 - alpha / 2, length(unique(data[[id]])) - 2)
}

# t <- pt_ttest_statistic_between(tmp, dv = "eda", between = "expviol",
#                                 time = "t")
# critical_t <- pt_critical_t_between(tmp, id = "id")

# # Plot
# t |>
#   ggplot(aes(x = name, y = value, group = 1)) +
#   geom_path() +
#   geom_hline(yintercept = c(-critical_t, critical_t))

# Calculation of the null distribution (H0) of cluster lengths for the
# between-subjects comparison
pt_null_distribution_between <- function(data, dv, between, time, id,
                                         nperm = 100)
{
  # Create a vector to store the null distribution of cluster lengths
  null_distribution <- vector("integer", nperm)
  
  critical_t <- pt_critical_t_between(data, id = id)
  
  # Number of samples per participant
  length_per_id <- length(unique(data[[time]]))
  
  # Number of participants per group
  n_per_group <- table(data[[between]]) / length_per_id
  
  # Keep factor levels of the between-subjects factor
  flevels <- levels(data[[between]])
  
  for (i in seq_len(nperm))
  {
    # Shuffle group identifier randomly
    data[[between]] <-
      rep(names(n_per_group), times = as.numeric(n_per_group)) |>
      sample() |>
      rep(each = length_per_id) |>
      factor(levels = flevels)
    
    # Find clusters in the random permutation
    t <- pt_ttest_statistic_between(data, dv = dv, between = between,
                                    time = time)
    clusters <- clusters_over_time(t, critical_t)
    
    # If clusters were found, take the one with the maximum length and save it
    null_distribution[i] <- if (nrow(clusters) > 0) max(clusters$length) else 0
    
  }
  
  null_distribution
}

# permutations <- pt_null_distribution_between(tmp, dv = "hr",
#                                              between = "expviol", time = "t",
#                                              id = "id", nperm = 5)

## Functions for within-subject comparisons (two conditions) -------------------

# This is an implementation of a cluster-based permutation test for within-
# subject comparisons (two conditions) of 2d data (signal x time)
# 
# The data needs to be in a tibble in long format with the following columns and
# must not have any missing data (NAs or missing rows, i.e., each subject needs
# to have 2 * "number of sample" rows):
# 
# (1) A column containing the dependent variable (dv) (e.g., heart rate, skin
#     conductance)
# (2) A column indicating the time points / samples (time)
# (3) A column containing the within-group identifier (within)
# (4) A column containing the participant identifier (id)
#
# All functions expect the column names as character.

# Function for calculating the t-statistic of a within-subject t-test. Used
# instead of `t.test(...)$statistic` for performance increase
t_statistic_within <- function(x, y)
{
  x <- x - y
  
  mean(x) / sqrt(var(x) / length(x))
}

# Calculate the t-statistic (within-subject comparison) for each time point
pt_ttest_statistic_within <- function(data, dv, within, time)
{
  data |>
    # Convert to wide format (i.e., values of the `dv` of both within-subject
    # conditions are stored in separate columns
    pivot_wider(names_from = !!within, values_from = !!dv) |>
    # For each sample in the signal ...
    group_by(!!sym(time)) |>
    # ... calculate the t-statistic
    summarize(value = t_statistic_within(
      # From the current subset of the tibble (i.e., per sample), extract the
      # column containing data of the first ...
      .data[[as.character(unique(data[[within]])[1])]],
      # ... and second level of the within-subjects condition
      .data[[as.character(unique(data[[within]])[2])]])) |>
    ungroup()
}

# Calculate the critical t-value for a paired t-test for the given sample size
pt_critical_t_within <- function(data, id, alpha = .05)
{
  qt(1 - alpha / 2, length(unique(data[[id]])) - 1)
}

# Calculation of the null distribution (H0) of cluster lengths for the
# within-subjects comparison
pt_null_distribution_within <- function(data, dv, within, time, id, nperm = 100)
{
  # Create a vector to store the null distribution of cluster lengths
  null_distribution <- vector("integer", nperm)
  
  critical_t <- pt_critical_t_within(data, id = id)
  
  for (i in seq_len(nperm))
  {
    # Make sure that `data` is arranged correctly
    data <- data |> arrange(!!!syms(c(id, within, time)))
    
    # Create a random permutation of the data (= permute the `within` column)
    data <-
      data |>
      # For each participant ...
      group_by(!!sym(id)) |>
      # ... shuffle the within-subject condition and repeat each for the number
      # of samples in the signal
      mutate("{within}" := rep(sample(unique(data[[within]])),
                               each = length(unique(data[[time]])))) |>
      ungroup()
    
    # Find clusters in the random permutation
    t <- pt_ttest_statistic_within(data, dv = dv, within = within, time = time)
    clusters <- clusters_over_time(t, critical_t)
    
    # If clusters were found, take the one with the maximum length and save it
    null_distribution[i] <- if (nrow(clusters) > 0) max(clusters$length) else 0
  }
  
  null_distribution
}

# Prepare data -----------------------------------------------------------------

# For between-subject comparisons, we prepare the data on the fly

# For within-subject comparisons, prepare data frames for comparisons between
# specific phases of the experiment

eda_within <- list(
  
  # Habituation vs. acquisition
  hab_vs_acq = eda |> filter(phase %in% c("habituation", "acquisition")),
  
  # Acquisition vs. extinction
  acq_vs_ext = eda |> filter(phase %in% c("acquisition", "extinction")),
  
  # Acquisition vs. spontaneous recovery
  acq_vs_spontrec = eda |> filter(phase %in% c("acquisition", "spontrec")),
  
  # Acquisition vs. reinstatement
  acq_vs_reinst = eda |> filter(phase %in% c("acquisition", "reinstatement"))
)

hr_within <- list(
  
  # Habituation vs. acquisition
  hab_vs_acq = hr |> filter(phase %in% c("habituation", "acquisition")),
  
  # Acquisition vs. extinction
  acq_vs_ext = hr |> filter(phase %in% c("acquisition", "extinction")),
  
  # Acquisition vs. spontaneous recovery
  acq_vs_spontrec = hr |> filter(phase %in% c("acquisition", "spontrec")),
  
  # Acquisition vs. reinstatement
  acq_vs_reinst = hr |> filter(phase %in% c("acquisition", "reinstatement"))
)

# Calculate null distributions -------------------------------------------------

# Use parallelization to run the permutations. As this was written for Windows,
# the doSNOW backend is used in the current implementation.

## Prepare parallelization ----------------------------------------------------

cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)

## EDA between -----------------------------------------------------------------

eda_null_dist_between <- setNames(vector("list", length(unique(eda$phase))),
                                  unique(eda$phase))

for (j in seq_along(eda_null_dist_between))
{
  # Extract the data of the current phase
  tmp <-
    eda |>
    filter(phase == names(eda_null_dist_between[j]))

  # Calculate null distribution
  x1 <- foreach(i = 1:7, .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_between(tmp, dv = "eda", between = "expviol",
                                      time = "t", id = "id", nperm = 143)
    return(y)
  }

  eda_null_dist_between[[j]] <- x1[1:1000]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
eda_critical_length_null_between <-
  tibble(phase = names(eda_null_dist_between),
         critical_length = map_dbl(eda_null_dist_between, quantile, .95))

## HR between ------------------------------------------------------------------

hr_null_dist_between <- setNames(vector("list", length(unique(hr$phase))),
                                 unique(hr$phase))

for (j in seq_along(hr_null_dist_between))
{
  # Extract the data of the current phase
  tmp <-
    hr |>
    filter(phase == names(hr_null_dist_between[j]))

  # Calculate null distribution
  x1 <- foreach(i = 1:7, .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_between(tmp, dv = "hr", between = "expviol",
                                      time = "t", id = "id", nperm = 143)
    return(y)
  }

  hr_null_dist_between[[j]] <- x1[1:1000]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
hr_critical_length_null_between <-
  tibble(phase = names(hr_null_dist_between),
         critical_length = map_dbl(hr_null_dist_between, quantile, .95))

## EDA within ------------------------------------------------------------------

eda_null_dist_within <- setNames(vector("list", length(eda_within)),
                                 names(eda_within))

for (j in seq_along(eda_within))
{
  # Extract the data of the current comparison
  tmp <- eda_within[[j]]
  
  # Calculate null distribution
  x1 <- foreach(i = 1:7, .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_within(tmp, dv = "eda", within = "phase",
                                     time = "t", id = "id", nperm = 143)
    return(y)
  }
  
  eda_null_dist_within[[j]] <- x1[1:1000]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
eda_critical_length_null_within <-
  tibble(phase = names(eda_null_dist_within),
         critical_length = map_dbl(eda_null_dist_within, quantile, .95))

## HR within -------------------------------------------------------------------

hr_null_dist_within <- setNames(vector("list", length(hr_within)),
                                names(hr_within))

for (j in seq_along(hr_within))
{
  # Extract the data of the current comparison
  tmp <- hr_within[[j]]
  
  # Calculate null distribution
  x1 <- foreach(i = 1:7, .combine = "c", .packages = "tidyverse") %dopar% {
    y <- pt_null_distribution_within(tmp, dv = "hr", within = "phase",
                                     time = "t", id = "id", nperm = 143)
    return(y)
  }
  
  hr_null_dist_within[[j]] <- x1[1:1000]
}

# Create a tibble with the 95% quantile of cluster lengths under H0 per phase
hr_critical_length_null_within <-
  tibble(phase = names(hr_null_dist_within),
         critical_length = map_dbl(hr_null_dist_within, quantile, .95))

### Save distributions and stop cluster ----------------------------------------

save(eda_null_dist_between, eda_critical_length_null_between,
     hr_null_dist_between, hr_critical_length_null_between,
     eda_null_dist_within, eda_critical_length_null_within,
     hr_null_dist_within, hr_critical_length_null_within,
     file = "RData/labexpvr_physio_pt_null_distributions.RData")

stopCluster(cl)

## Cluster-based permutation tests ---------------------------------------------

load("RData/labexpvr_physio_pt_null_distributions.RData")

### EDA ------------------------------------------------------------------------

# Find clusters of differences in the EDA signal in each phase of the experiment
eda_clusters_between <-
  eda |>
  (\(x) split(x, x$phase))() |>
  map(pt_ttest_statistic_between, dv = "eda", between = "expviol",
      time = "t") |>
  map(clusters_over_time, critical_t = pt_critical_t_between(eda, id = "id")) |>
  # Bind to a single tibble
  bind_rows(.id = "phase") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(eda_critical_length_null_between, by = "phase") |>
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, eda_null_dist_between[[phase]])) |>
  ungroup()

# No sign. clusters

### HR -------------------------------------------------------------------------

# Find clusters of differences in the HR signal in each phase of the experiment
hr_clusters_between <-
  hr |>
  (\(x) split(x, x$phase))() |>
  map(pt_ttest_statistic_between, dv = "hr", between = "expviol",
      time = "t") |>
  map(clusters_over_time, critical_t = pt_critical_t_between(hr, id = "id")) |>
  # Bind to a single tibble
  bind_rows(.id = "phase") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(hr_critical_length_null_between, by = "phase") |>
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, hr_null_dist_between[[phase]])) |>
  ungroup()

# No sign. clusters

## Save significant clusters ---------------------------------------------------

save(eda_clusters, hr_clusters, file = "RData/labexpvr_physio_clusters.RData")

# Within-group comparison (compare different phases) ---------------------------


## Cluster-based permutation tests ---------------------------------------------

load("RData/labexpvr_physio_cluster_mass_distributions.RData")

## Cluster-based permutation tests (within) ------------------------------------

### EDA ------------------------------------------------------------------------

# Find clusters of differences in the EDA signal in each phase of the experiment
eda_clusters_within <-
  eda_within |>
  map(pt_ttest_statistic_within, dv = "eda", within = "phase", time = "t") |>
  map(clusters_over_time, critical_t = pt_critical_t_within(eda, id = "id")) |>
  # Bind to a single tibble
  bind_rows(.id = "phase") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(eda_critical_length_null_within, by = "phase") |>
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, eda_null_dist_within[[phase]])) |>
  ungroup()

### HR -------------------------------------------------------------------------

# Find clusters of differences in the EDA signal in each phase of the experiment
hr_clusters_within <-
  hr_within |>
  map(pt_ttest_statistic_within, dv = "hr", within = "phase", time = "t") |>
  map(clusters_over_time, critical_t = pt_critical_t_within(hr, id = "id")) |>
  # Bind to a single tibble
  bind_rows(.id = "phase") |>
  # Add critical cluster length under H0 (95% quantile)
  left_join(hr_critical_length_null_within, by = "phase") |>
  # Calculate the p-value for each cluster
  rowwise() |>
  mutate(p = 1 - pempiricalD(length, hr_null_dist_within[[phase]])) |>
  ungroup()
