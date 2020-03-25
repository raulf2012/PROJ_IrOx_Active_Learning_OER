import os
import sys


listdirs = os.listdir(".")

job_dirs = [i for i in listdirs if i.isnumeric()]
print("len(job_dirs):", len(job_dirs))

completed_ids = []
for dir_i in job_dirs:

    # | - Getting path of max revision
    rev_dirs_list = os.listdir(dir_i)

    rev_dirs_list = [i for i in rev_dirs_list if i[0] == "_"]

    if len(rev_dirs_list) == 0:
        print(dir_i)
        continue

    tmp0 = [int(i.split("_")[-1]) for i in rev_dirs_list]

    max_rev_i = max(tmp0)
    max_rev_path = os.path.join(dir_i, "_" + str(max_rev_i))
    #__|

    # | - Checking isif
    isif_2_bool = False

    with open(max_rev_path + "/model.py") as fle:
        model_lines = fle.readlines()

    isif_lines = [i.strip() for i in model_lines if "isif" in i]
    # isif_lines = [i.strip() for i in isif_lines]
    isif_lines_2 = [i for i in isif_lines if i[0:4] == "isif"]

    assert len(isif_lines_2) == 1, "sijfisd"
    isif_i = isif_lines_2[0].split("=")[-1].split(",")[0]
    isif_i = int(isif_i)

    if isif_i == 2:
        isif_2_bool = True
    #__|

    # | - Checking if .FINISHED file is present
    finished_file_path = os.path.join(max_rev_path, ".FINISHED")
    finished_bool = os.path.isfile(finished_file_path)
    #__|

    if isif_2_bool and finished_bool:
        completed_ids.append(dir_i)

print("len(completed_ids):", len(completed_ids))


# Pickling data ######################################################
import os; import pickle
directory = "../out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "parse_data.pickle"), "wb") as fle:
    pickle.dump((job_dirs, completed_ids), fle)
# #####################################################################
