# PROJ_IrOx_Active_Learning_OER

























# Instructions to run all code within this project

Just a note that setting up the conda environment as specified is critical, otherwise you will run into a lot of cross-compatibility errors. Recreating the exact environment I used (with conda) is the most robust way of making sure everything runs properly.

## Clone this repo into your computer

git clone https://github.com/raulf2012/PROJ_IrOx_Active_Learning_OER.git

Copy the files `config.yml.dist` and `setup.sh.dist` in `PROJ_IrOx_Active_Learning_OER/config/setup.sh.dist` to files called `setup.sh`  and `config.yml` in the same folder.

Fill in the relevent info into `config.yml` (currently only the MP API key is needed) and also start filling in the `setup.sh` file with info (Define the PROJ_irox ENV variable which points to the `PROJ_IrOx_Active_Learning_OER` repo path, PROJ_DATA will also need to be filled in later).

## Download the PROJ_DATA data folder

Data which houses about 7.5 GB of raw data that is needed to run scripts. I can try to decrease the size of the download a bit if needed.

Go to the following link and download it into zip file. Put the zip file into a folder called `04_IrOx_surfaces_OER` and unzip it.


This is a link to the zipped PROJ_DATA files, should be quicker to download this one:

https://www.dropbox.com/s/rs99pfsch7orpiy/04_IrOx_surfaces_OER.zip?dl=0

https://www.dropbox.com/sh/tbe9u21igj43l4q/AACnwx655hFcJwaA8dwv52bTa?dl=0

Note: Don't use Microsoft Edge to attempt the download, it stalls for some reason.

You will then need to go back to `setup.sh` and define the PROJ_DATA variable to the location of the `04_IrOx_surface_OER` directory (not including `04_IrOx_surface_OER`).

## Install conda package manager
Conda is used to install all required python packages in an independent virtual environment

A minimal conda installation can be installed as follows on Linux:

`curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda_install.sh`

But replacing the correct URL with your system (Mac, Windows, Linux). Go to Conda's website for more info on installation:

https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

`bash miniconda_install.sh`


## Run the installer bash script

The file `config/clone_needed_repos.sh` will go through the procedure of cloning needed repo's and installing a working conda environment. The repos will be default be installed into $HOME/PROJ_irox_repos but that can be changed by specifying the root_folder variable in the installation script. This will have to be changed in the setup.sh file as well.


`bash clone_needed_repos.sh`

This will take a while, installing conda and creating the necessary environment (with all the required packages) is time intensive, maybe 20-30 minutes.



## Source the setup.sh file

The `setup.sh` script will activate the conda env, set needed environment variables, and append to PYTHONPATH.

Run the file by sourcing it:

`source setup.sh`

I would just put this into your `.bashrc` to save time.


## Run all notebook files at once

Start a jupyterlab instance by first sourcing `setup.sh`, navigating to the root repo dir, and running jupyter lab. Then copy and paste the resulting URL into a browswer, they usually look lik this:

http://localhost:8888/?token=48ecf18a1544c28cb1b12a9a6670689e2f3bea88200bf27a


From here navigate to the jupyter notebook at 

`PROJ_IrOx_Active_Learning_OER/scripts/run_all_notebooks_1.ipynb`

which will run through all the notebooks. Importantly, this has to be done in the correct order, otherwise data dependencies will not be met and everything will crash. The notebook will take roughly 20 minutes to run completely. As notebooks are run their output will be redirected to the running notebook, if they finish correctly they will terminate by printing a "All done!" line, watch for these as the notebook runs, if there is an error go investigate the file in question, open it in Jupyter Lab and try running it. Hopefully I've ironed out any possible points of failure and the failure is due to a mistyped path variable or something of the like.

The `run_all_notebooks_1.ipynb` notebook will display a dataframe (`df`) that will display whether each notebook completed succesfully (see `completed` column). The desired outcome here is that every single notebook has a `True` value here.

Once this is done, you can than dig more deeply into the notebooks, running them individually, proding, modifying, etc. The main point of the run all notebooks notebook was to make sure that all data is generated at least once (data dependencies are satisfied).

To save time I have coded in a shortcut, by which time intensive data objects are simply read from pre-created files instead of being created in real-time. This behavior is controlled by the variable `read_from_PROJ_DATA` and is located in `data/proj_data_irox.py`. Turn this joff to run everything more manually.











# Notebooks that can be run immedietely

workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/create_atoms_df.ipynb































## `pinned` file in the following dir:

$HOME/anaconda3/envs/PROJ_IrOx_Active_Learning_OER/conda-meta/pinned

Add this following line

tensorflow ==1.15.*
