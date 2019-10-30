#!/usr/bin/env python

import os
import time

start_time = time.time()

python_ex = "/global/homes/f/flores12/anaconda2/envs/py36/bin/python"
scratch_dir = "/global/cscratch1/sd/flores12"
bulk_calc_dir = "/IrOx_Project_temp_190510/ml_bulk_calculations"

num_loops = 18
loop_list = list(range(num_loops))
for i in loop_list:
    print("loop iteration number: ", i)

    time_i = time.time()

    bash_comm = python_ex + \
        " " + scratch_dir + \
        bulk_calc_dir + "/iro3_calcs/rerun_jobs.py"
    os.system(bash_comm)


    time_j = time.time()
    loop_time = time_j - time_i
    print("loop time (sec): ", loop_time)
    print("loop time (min): ", loop_time / 60)

    if i == loop_list[-1]:
        print("ALL DONE WITH LOOP")
    else:
        time.sleep(4000)
        print("Finished waiting")

end_time = time.time()

time_elapsed_sec = end_time - start_time
time_elapsed_min = time_elapsed_sec / 60

print("time_elapsed (sec): ", time_elapsed_sec)
print("time_elapsed (min): ", time_elapsed_min)
