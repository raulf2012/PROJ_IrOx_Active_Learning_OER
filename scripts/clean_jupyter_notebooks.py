"""
"""

#| - Import Modules
import os
import sys

from json import dump, load
from shutil import copyfile

# from .clean_julia_code import clean_julia_code
# from .clean_python_code import clean_python_code
#__|


def clean_ipynb(ipynb_file_path, overwrite):
    """
    """
    #| - clean_ipynb
    if not overwrite:
        ipynb_file_path = copyfile(
            ipynb_file_path,
            ipynb_file_path.replace(".ipynb", ".cleaned.ipynb")
            )

    with open(ipynb_file_path) as io:
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


# ipynb_file_path = "ml_plots__v3.ipynb"
# clean_ipynb(ipynb_file_path, True)




dirs_to_ignore = [
    ".ipynb_checkpoints",
    "__old__",
    "__temp__",
    "__test__",
    "__bkp__",
    ]

root_dir = os.path.join(
    os.environ["PROJ_irox"],
    # "workflow/tmp_sandbox",
    )
for subdir, dirs, files in os.walk(root_dir):

    #| - Pass over ignored directories
    ignore_subdir = False
    for ignore_dir in dirs_to_ignore:
        if ignore_dir in subdir:
            ignore_subdir = True
    if ignore_subdir:
        continue
    #__|

    for file in files:
        if ".ipynb" in file:
            file_path = os.path.join(subdir, file)
            print(file_path)

            clean_ipynb(file_path, True)