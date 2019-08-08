import os
import time

#| - Job ID List
job_id_list = [
    "10977313",
    "10977314",
    "10977315",
    "10977332",
    "10977333",
    "10977335",
    "10977456",
    "10977458",
    "10977459",
    "10978865",
    "10978890",
    "10978905",
    "10980003",
    "10980004",
    "10980005",
    "10980064",
    "10980097",
    "10980098",
    "10980100",
    "10980111",
    "10980113",
    "10980122",
    "10980123",
    "10980124",
    "10980126",
    "10980128",
    "10980129",
    "10980136",
    "10980137",
    "10980138",
    "10980139",
    "10980140",
    "10980141",
    "10992748",
    "10992771",
    "10992811",
    "10992849",
    "10992864",
    "10992867",
    "10992902",
    "10992911",
    "10993046",
    "10993065",
    "10993071",
    "10993081",
    ]
#__|

# Command
# scontrol update jobid=10977086 partition=regular

for job_i in job_id_list:
    print(job_i)

    command_i = "scontrol update jobid=" + job_i + " partition=regular"
    os.system(command_i)
    time.sleep(4)
