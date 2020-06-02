# PROJ_IrOx_Active_Learning_OER

#TODO

TODO Add labels to OER volcano plot for exp. overpot. lines

TODO Add Columbite IrO2 data to everything #imp


























# Instructions to run all code within this project

Just a note that setting up the conda environment as specified is critical, otherwise you will run into a lot of cross-compatability errors. Recreating the exact environment I used (with conda) is the most robust way of making sure everything runs properly.

## Clone this repo into your computer

git clone https://github.com/raulf2012/PROJ_IrOx_Active_Learning_OER.git

Copy the files `config.yml.dist` and `setup.sh.dist` in `PROJ_IrOx_Active_Learning_OER/config/setup.sh.dist` to files called `setup.sh`  and `config.yml` in the same folder.

## Download the PROJ_DATA data folder

Data which houses about 7.5 GB of raw data that is needed to run scripts. I can try to decrease the size of the download a bit if needed.

Go to the following link and download it into zip file. Put the zip file into a folder called `04_IrOx_surfaces_OER` and unzip it.

https://www.dropbox.com/sh/tbe9u21igj43l4q/AACnwx655hFcJwaA8dwv52bTa?dl=0

Note: Don't use Microsoft Edge to attempt the download, it stalls for some reason.

## Install conda package manager
Conda is used to install all required python packages in an independent virtual environment

A minimal conda installation can be installed as follows on Linux:

`curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda_install.sh`
But replacing the correct URL with your system (Mac, Windows, Linux)

`bash miniconda_install.sh`


## Download required git repos to local directory

### Download my personal PythonModules repo
This has methods that are needed for this Project repo
Clone the release version 2.0 of 00_PythonModules into a location of your choice

git clone --branch 2.0 https://github.com/raulf2012/00_PythonModules.git

### Download CatLearn package
The Gaussian Process regression module is from CatLearn.

git clone --branch 0.6.0 https://github.com/SUNCAT-Center/CatLearn.git

### Download Structure Prototype Analysis Package
Used for structure similrity analysis

git clone https://github.com/chuanxun/StructurePrototypeAnalysisPackage.git --branch v1.0.1

### Download Protosearch codebase

https://github.com/SUNCAT-Center/Protosearch

git clone --branch 0.0.1 https://github.com/SUNCAT-Center/Protosearch.git


# Setting up an isolated conda environment:



## Create a new blank environment

conda create --name PROJ_irox python=3.6 --no-default-packages



## conda install commands
    conda install \
    -c plotly \
    -c conda-forge \
    -c anaconda \
    nodejs jupyterlab jupytext \
    scikit-learn matplotlib scipy pandas \
    plotly chart-studio plotly-orca psutil colormap colorlover \
    ase=3.17.0 pymatgen gpflow \
    nodejs nb_conda_kernels tensorflow==1.15 \
    yaml nbclean

## Install plotly jupyter lab extension

This is needed so that figures appear in notebooks

jupyter labextension install jupyterlab-plotly






# Notebooks that can be run immedietely

workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/create_atoms_df.ipynb































## `pinned` file in the following dir:

$HOME/anaconda3/envs/PROJ_IrOx_Active_Learning_OER/conda-meta/pinned

Add this following line

tensorflow ==1.15.*
