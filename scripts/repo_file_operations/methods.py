"""
"""

# | - Import Modules
import os
import sys

from json import dump, load
from shutil import copyfile

import subprocess
#__|



def clean_ipynb(ipynb_file_path, overwrite):
    """
    """
    # | - clean_ipynb
    if not overwrite:
        ipynb_file_path = copyfile(
            ipynb_file_path,
            ipynb_file_path.replace(".ipynb", ".cleaned.ipynb")
            )

    # with open(ipynb_file_path) as io:
    with open(ipynb_file_path, "rb") as io:
        ipynb_dict = load(io)


    # language = ipynb_dict["metadata"]["language_info"]["name"]
    # if language == "python":
    #     clean_code = clean_python_code
    # elif language == "julia" and can_clean_julia_code():
    #     clean_code = clean_julia_code
    # else:
    #     return

    cells = []
    for cell_dict in ipynb_dict["cells"]:

        cell_dict["execution_count"] = None
        cell_dict["outputs"] = []

        if (
            "metadata" in cell_dict
            and "jupyter" in cell_dict["metadata"]
            and "source_hidden" in cell_dict["metadata"]["jupyter"]
            ):
            cell_dict["metadata"]["jupyter"].pop("source_hidden")

        if cell_dict["cell_type"] == "code":

            tmp = 42

            # source_join_clean_split = clean_code(
            #     "".join(cell_dict["source"])
            #     ).split(sep="\n")
            #
            # if len(source_join_clean_split) == 1 and source_join_clean_split[0] == "":
            #     continue
            #
            # source_join_clean_split[:-1] = [
            #     "{}\n".format(line) for line in source_join_clean_split[:-1]
            #     ]
            #
            # cell_dict["source"] = source_join_clean_split

        cells.append(cell_dict)

    ipynb_dict["cells"] = cells

    with open(ipynb_file_path, mode="w") as io:
        dump(ipynb_dict, io, indent=1)
        io.write("\n")
    #__|


def get_ipynb_notebook_paths(
    PROJ_irox_path=None,
    relative_path=False,
    ):
    """
    """
    #| - get_ipynb_notebook_paths
    dirs_to_ignore = [
        ".ipynb_checkpoints",
        "__old__",
        "__temp__",
        "__test__",
        "__bkp__",
        ]

    if PROJ_irox_path is None:
        root_dir = os.path.join(
            os.environ["PROJ_irox"])
    else:
        root_dir = PROJ_irox_path


    res = subprocess.check_output(
        "git rev-parse --show-toplevel".split(" ")
        )
    root_git_path = res.decode("UTF-8").strip()

    dirs_list = []
    for subdir, dirs, files in os.walk(root_dir):

        # | - Pass over ignored directories
        ignore_subdir = False
        for ignore_dir in dirs_to_ignore:
            if ignore_dir in subdir:
                ignore_subdir = True
        if ignore_subdir:
            continue
        #__|

        for file in files:
            if ".swp" in file:
                continue

            if file == "convert_jup_to_pyth.ipynb":
                continue

            if ".ipynb" in file:
                file_path = os.path.join(subdir, file)
                full_dir_i = os.path.join(subdir, file)

                dir_i = full_dir_i
                if relative_path:
                    file_i = full_dir_i
                    find_ind = file_i.find(root_git_path)
                    if find_ind == 0:
                        relative_path_i = file_i[len(root_git_path) + 1:]
                    else:
                        print("Problem!! sdijfs")
                        relative_path_i = "Yikes!"

                    dir_i = relative_path_i



                dirs_list.append(dir_i)

    return(dirs_list)
    #__|



#| - __old__
# ipynb_file_path = "ml_plots__v3.ipynb"
# clean_ipynb(ipynb_file_path, True)
#__|
