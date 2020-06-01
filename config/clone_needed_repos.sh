#!/bin/bash

root_folder=$HOME/PROJ_irox_repos

mkdir $root_folder

#| - Clone Repos

#| - PythonModules
mkdir $root_folder/00_PythonModules

git clone --branch 2.0 https://github.com/raulf2012/00_PythonModules.git $root_folder/00_PythonModules
#__|

echo "****************************************"

#| - CatLearn
mkdir $root_folder/CatLearn

git clone --branch 0.6.0 https://github.com/SUNCAT-Center/CatLearn.git $root_folder/CatLearn
#__|

echo "****************************************"

#| - StructurePrototypeAnalysisPackage
mkdir $root_folder/StructurePrototypeAnalysisPackage

git clone --branch v1.0.1 https://github.com/chuanxun/StructurePrototypeAnalysisPackage.git $root_folder/StructurePrototypeAnalysisPackage
#__|

echo "****************************************"

#| - Protosearch
mkdir $root_folder/Protosearch
git clone https://github.com/SUNCAT-Center/Protosearch.git $root_folder/Protosearch

cd $root_folder/Protosearch
git checkout b7827449d85ed074c623c422c9b3f814715694a5
#__|

echo "****************************************"

#__|


#| - Conda Environment

conda create --name PROJ_irox python=3.6 --no-default-packages

echo "****************************************"

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate PROJ_irox

echo "****************************************"

conda info

echo "****************************************"

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


jupyter labextension install jupyterlab-plotly
#__|

