# PROJ_IrOx_Active_Learning_OER

#TODO

TODO Add labels to OER volcano plot for exp. overpot. lines

TODO Add Columbite IrO2 data to everything #imp


























# Instructions to run all code within this project

## Download the PROJ_DATA data folder

Data which houses about TEMP GB of raw data that is needed to run scripts




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

https://github.com/chuanxun/StructurePrototypeAnalysisPackage

### Download Protosearch codebase

https://github.com/SUNCAT-Center/Protosearch

git clone --branch 0.0.1 https://github.com/SUNCAT-Center/Protosearch.git
git clone https://github.com/SUNCAT-Center/Protosearch.git


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
