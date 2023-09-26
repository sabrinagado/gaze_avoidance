import os
from pycasso import Canvas

project_dir = os.getcwd()
file_path = os.path.join(project_dir, 'stimuli', 'social')
file_path_ns = os.path.join(project_dir, 'stimuli', 'non-social')
files = [item for item in os.listdir(file_path) if (item.endswith(".jpg"))]
for file in files:
    # file = files[0]
    Canvas(os.path.join(file_path, file), (10, 10), 'seed').export('scramble', os.path.join(file_path_ns, file + "_scrambled"), 'jpeg')
