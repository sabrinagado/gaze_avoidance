# =============================================================================
# Scores
# source: SosciSurvey
# study: Gaze Camouflage
# =============================================================================
import os
import pandas as pd
import numpy as np
import random
from scipy.stats import linregress, ttest_ind, chi2_contingency, zscore, pearsonr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns


dir_path = os.getcwd()
df_scores = pd.read_excel(os.path.join(dir_path, "questionnaires.xlsx"), usecols=list(np.arange(6, 165)))

# Recode variables from numeric to string
df_scores['gender'] = df_scores['gender'].replace({1: "male", 2: "female", 3: "diverse"})
df_scores['handedness'] = df_scores['handedness'].replace({1: "right", 2: "left"})
df_scores['smoking'] = df_scores['smoking'].replace({1: "smoker", 2: "no smoker"})

# % ===========================================================================
# SPAI
# =============================================================================
df_spai = df_scores.filter(like='spai')
# Adapt scaling (from 1-7 to 0-6)
df_spai = df_spai - 1

# Calculate means of nested items
columns_spai_multi = [col for col in df_spai.columns if '_' in col]
df_spai_multi = df_spai[columns_spai_multi]
columns_spai = [col for col in df_spai.columns if not '_' in col]
df_spai = df_spai[columns_spai]
items = list(dict.fromkeys([int(column.split('spai')[1].split('_')[0]) for column in df_spai_multi.columns]))
for item in items:
    # item = items[0]
    df_spai_multi_subset = df_spai_multi.filter(like=f'{item}_')
    df_spai[f'spai_{item}'] = df_spai_multi_subset.mean(axis=1)

# Calculate mean of scale
df_spai['SPAI'] = df_spai.mean(axis=1)
df_spai = df_spai[['SPAI']]

# % ===========================================================================
# SIAS
# =============================================================================
df_sias = df_scores.filter(like='sias')
# Adapt scaling (from 1-5 to 0-4)
df_sias = df_sias - 1

# Calculate sum of items
df_sias['SIAS'] = df_sias.sum(axis=1)
df_sias = df_sias[['SIAS']]

# % ===========================================================================
# STAI-T
# =============================================================================
df_stai = df_scores.filter(like='stai')

# Calculate sum of respective items
df_stai['STAI-T'] = df_stai.sum(axis=1)
df_stai = df_stai[['STAI-T']]

# % ===========================================================================
# Intolerance of Uncertainty
# =============================================================================
df_ui = df_scores.filter(like='ui')

# Calculate sum of respective items
df_ui['UI'] = df_ui.sum(axis=1)
df_ui = df_ui[['UI']]

# % ===========================================================================
# VAS: State Anxiety, Nervousness, Distress, Stress
# =============================================================================
df_vas = df_scores.filter(like='VAS')
df_vas["VAS_diff_anxiety"] = df_vas["VAS_end_anxiety"] - df_vas["VAS_start_anxiety"]
df_vas["VAS_diff_nervous"] = df_vas["VAS_end_nervous"] - df_vas["VAS_start_nervous"]
df_vas["VAS_diff_distress"] = df_vas["VAS_end_distress"] - df_vas["VAS_start_distress"]
df_vas["VAS_diff_stress"] = df_vas["VAS_end_stress"] - df_vas["VAS_start_stress"]

# % ===========================================================================
# Debriefing / Control Questions
# =============================================================================
df_debriefing = df_scores[['purpose', 'variables']]
df_debriefing['purpose'] = df_debriefing['purpose'].str.replace('-', '')
df_debriefing['purpose'] = df_debriefing['purpose'].str.replace('\n', ', ')
df_debriefing['variables'] = df_debriefing['variables'].str.replace('-', '')
df_debriefing['variables'] = df_debriefing['variables'].str.replace('\n', ', ')

# % ===========================================================================
# Labbook
# =============================================================================
df_labbook = df_scores[['labbook', 'digitimer', 'temperature', 'humidity']]

# % ===========================================================================
# Demographic Information
# =============================================================================
df_scores['VP'] = df_scores['ID'].astype('int32')
df_demo = df_scores[['VP', 'gender', 'age', 'motivation', 'tiredness', 'handedness']]

# % ===========================================================================
# Create Summary
# =============================================================================
df = pd.concat([df_demo, df_debriefing, df_labbook, df_vas, df_spai, df_sias, df_stai, df_ui], axis=1)
df.to_csv(os.path.join(dir_path, 'scores_summary.csv'), index=False, decimal='.', sep=';', encoding='utf-8-sig')

# % ===========================================================================
# Sample Description
# =============================================================================
print(f"N = {len(df)}")
print(f"Mean Age = {df['age'].mean()}, SD = {df['age'].std()}, Range = {df['age'].min()}-{df['age'].max()}")
print(df['gender'].value_counts(normalize=True))

# Plot Distributions
for scale in ["SPAI", "SIAS", "STAI-T"]:
    df_plot = df.copy()
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 1.6))
    colors = ['#FF5733', '#FFC300', '#183DB2']
    if scale == "SPAI":
        sns.histplot(df_plot[scale], color="#1B4C87", ax=ax, binwidth=0.2, binrange=(0, 5),
                     kde=True, line_kws={"linewidth": 1, "color": "#173d6a"}, edgecolor='#f8f8f8',)
    elif scale == "SIAS":
        sns.histplot(df_plot[scale], color="#1B4C87", ax=ax, binwidth=2.5, binrange=(0, 60),
                     kde=True, line_kws={"linewidth": 1, "color": "#173d6a"}, edgecolor='#f8f8f8',)
    elif scale == "STAI-T":

        sns.histplot(df_plot[scale], color="#1B4C87", ax=ax, binwidth=2.5, binrange=(20, 70),
                     kde=True, line_kws={"linewidth": 1, "color": "#173d6a"}, edgecolor='#f8f8f8',)
    ax.set_xlabel(scale)
    ax.axvline(x=df_plot[scale].median(), color="#FFC300")
    if scale == "SPAI":
        ax.axvline(x=df_plot["SPAI"].median(), color="#FFC300")
        ax.axvline(x=2.79, color="#FF5733")
        ax.legend(
                [Line2D([0], [0], color='#FFC300'), Line2D([0], [0], color='#FF5733')],
                ['Median', 'Remission Cut-Off'], fontsize='xx-small', loc="best", frameon=False)
    ax.set_facecolor('#f8f8f8')
    plt.tight_layout()
    plt.savefig(os.path.join(dir_path, "plots", f"Distribution_{scale}.png"), dpi=300)
    plt.close()
