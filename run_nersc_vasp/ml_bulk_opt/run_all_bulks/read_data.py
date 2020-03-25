import sys
import pickle
import os
import pandas as pd

pd.set_option('display.max_rows', None)



# | - Get command line argument
filename = sys.argv[-1]

#__|

# path_i = os.path.join("df_dict.pickle")
path_i = os.path.join(filename)
with open(path_i, "rb") as fle:
    data = pickle.load(fle)


df_new_jobs = data["df_new_jobs"]
df_new_jobs = df_new_jobs.drop("pre_path", axis=1)

df_not_done = df_new_jobs[df_new_jobs["action"] != "ALL DONE! | ISIF 2"]

df_not_done_busy = df_not_done[df_not_done["action"] != "Job is busy, will skip"]
