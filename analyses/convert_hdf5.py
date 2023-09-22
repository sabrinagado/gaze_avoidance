# =============================================================================
# HDF5
# study: GCA
# =============================================================================
import os
import h5py
import pandas as pd

dir_path = os.getcwd()
file_path = os.path.join(dir_path, 'data')
filenames = [item for item in os.listdir(file_path) if (item.endswith(".hdf5"))]

for filename in filenames:
    # filename = filenames[0]
    id = filename.split("_")[0]
    with h5py.File(os.path.join(file_path, filename), "r") as f:
        # get the list of eyetracker measures available in the hdf5
        eyetracker_measures = list(f['data_collection']['events']['eyetracker'])

        for measure in eyetracker_measures:
            print('Extracting events of type: ', measure)
            data_collection = list(f['data_collection']['events']['eyetracker'][measure])
            if len(data_collection) > 0:
                column_headers = data_collection[0].dtype.descr
                cols = []
                data_dict = {}
                for ch in column_headers:
                    cols.append(ch[0])
                    data_dict[ch[0]] = []

                for row in data_collection:
                    for i, col in enumerate(cols):
                        data_dict[col].append(row[i])
                pd_data = pd.DataFrame.from_dict(data_dict)
                pd_data.to_csv(os.path.join(file_path, id + '_' + measure + '.csv'), index=False, decimal=',', sep=';')
            else:
                print('No data for type', measure, ' moving on')


df = pd.read_csv(os.path.join(file_path, '0_gca_2023-09-21_14h46.53.285.csv'), decimal='.', sep=',')
start_grid = df["cross_1.started"].dropna().item()
end_grid = start_grid + 10

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

save_path = os.path.join(dir_path, 'analyses', 'plots')
if not os.path.exists(save_path):
    print('creating path for saving')
    os.makedirs(save_path)

filename = '0_MonocularEyeSampleEvent.csv'

# read as pandas dataframe
data = pd.read_csv(os.path.join(file_path, filename), decimal=',', sep=';')
data = data.loc[(data["time"] > start_grid) & (data["logged_time"] < end_grid)]

# convert pandas arrays to no arrays
x = data['gaze_x'].to_numpy()
y = data['gaze_y'].to_numpy()

# remove nan values
x = x[~np.isnan(x)]
y = y[~np.isnan(y)]

plt.plot(x, y, marker="+", linestyle="")
plt.savefig(os.path.join(save_path, f"grid_test.png"), dpi=300)
plt.close()

