How to import the yaml config file in a python script:

    import yaml
    path_i = os.path.join(os.environ["PROJ_irox"], "config", "config.yml")
    with open(path_i) as file:
        config_dict = yaml.load(file, Loader=yaml.FullLoader)

    api_key = config_dict["mpcontrib"]["api_key"]


# Required repos info

## Download required git repos to local directory

### Download my personal PythonModules repo
This has methods that are needed for this Project repo
Clone the release version 2.0 of 00_PythonModules into a location of your choice

git clone --branch 2.0 https://github.com/raulf2012/00_PythonModules.git

### Download CatLearn package
The Gaussian Process regression module is from CatLearn.

git clone --branch v0.6.0 https://github.com/SUNCAT-Center/CatLearn.git

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
