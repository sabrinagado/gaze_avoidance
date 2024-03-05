# =============================================================================
# Create File with Group of Participants
# study: GCA
# =============================================================================
import os
import pandas as pd
import numpy as np

np.random.seed(42)

conditions = ["A", "B", "C", "D"]
N = 52
vps = np.arange(1, N+1)

groups = conditions * int(N/len(conditions))
np.random.shuffle(np.array(groups))

df = pd.DataFrame({'VP': list(vps), 'Group': list(groups)})
df.to_csv('condition_list.csv', decimal='.', sep=';', index=False)
