# =============================================================================
# Power Analysis
# study: GCA
# =============================================================================
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import statsmodels.stats.power as pwr
import math

"""
Hypotheses:
- More fixations towards the social CS+ compared to the non-social CS+
"""

colors = ['#1F82C0', '#F29400', '#E2001A', '#B1C800', '#179C7D']
dir_path = os.getcwd()
save_path = os.path.join(dir_path, 'analyses', 'plots')
if not os.path.exists(save_path):
    print('creating path for saving')
    os.makedirs(save_path)

# Factors for power analysis
alpha = 0.05
power = 0.8
N = 52
pwr_analysis = pwr.TTestPower()

np.random.seed(42)

# =============================================================================
# Number of Fixations
# =============================================================================
# parameters for power analysis
# TODO: Check for plausbility
sd_percentage = 10
n = N*2
prob_cs_plus_soc = 12
percentage_cs_plus_soc = np.random.normal(prob_cs_plus_soc, sd_percentage, n)
prob_cs_plus_nsoc = 10
percentage_cs_plus_nsoc = np.random.normal(prob_cs_plus_nsoc, sd_percentage, n)

values = np.concatenate([percentage_cs_plus_soc, percentage_cs_plus_nsoc])
std = np.std(values)
differences = np.arange(2, 15, 0.1)
effect_sizes = differences/std
sample_sizes = []

# calculate sample sizes based on effect sizes
for effect_size in effect_sizes:
    sample_sizes.append(pwr_analysis.solve_power(effect_size=effect_size, alpha=alpha, power=power, alternative="two-sided"))

# plot
dv = "Difference between Fixations on social CS+ and non-social CS+"
measure = "Fixation in % of Trials"
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5))
ax.plot(differences, sample_sizes, color=colors[0])
ax.grid(color='lightgrey', linestyle='-', linewidth=0.3)
ax.set_title(f"Power Analysis:\n {dv}", fontweight="bold", fontsize="x-large")
ax.set_xlabel(f"Difference [{measure}]")
ax.set_ylabel(f"Required Sample Size")
ax.set_ylim([np.min(sample_sizes), np.max(sample_sizes)])
ax.set_xlim([np.min(differences), np.max(differences)])

# # N = X:
effect = pwr_analysis.solve_power(nobs=N, alpha=alpha, power=power, alternative="two-sided") * std
ax.axhline(N, color=colors[4], linestyle="--")
ax.text(effect + 0.01 * np.max(differences), N + 0.01 * np.max(sample_sizes), f"N: {math.ceil(N)}", color=colors[4])
ax.axvline(effect, color=colors[4], linestyle="--")
ax.text(effect + 0.01 * np.max(differences), 0.65 * np.max(sample_sizes), f"Effect: Difference = {round(effect, 2)}%", color=colors[4])

ax.legend(
    [Line2D([0], [0], color="white"),
     Line2D([0], [0], color="white"),
     Line2D([0], [0], color="white"),
     Line2D([0], [0], color="white")],
    [f"Probability of social CS+ Fixation = {prob_cs_plus_soc}%",
     f"Probability of non-social CS+ Fixation = {prob_cs_plus_nsoc}%",
     f"SD of CS+ Fixations = {sd_percentage}%", f"N = {n}"], loc="upper right")

plt.tight_layout()
plt.savefig(os.path.join(save_path, f"pwr_fixations_{round(prob_cs_plus_soc, 2)}_{round(prob_cs_plus_nsoc, 2)}_{round(sd_percentage, 2)}_{n}.png"), dpi=300)
plt.close()
