@echo off

REM set TEMP_VAR=F:\Dropbox
set DROPBOX_DIR=C:\Users\raulf2012\Dropbox

setx PROJ_irox %DROPBOX_DIR%"\01_norskov\00_git_repos\PROJ_IrOx_Active_Learning_OER"
setx PROJ_DATA %DROPBOX_DIR%"\01_norskov\PROJECT_DATA"

REM setx PROJ_irox "F:\Dropbox\01_norskov\00_git_repos\PROJ_IrOx_Active_Learning_OER"
REM setx PROJ_DATA "F:\Dropbox\01_norskov\PROJECT_DATA"


REM setx PYTHONPATH "F:\Dropbox\01_norskov\00_PythonModules;F:\Dropbox\01_norskov\00_git_repos\CatLearn"
REM setx PYTHONPATH "%PYTHONPATH%;F:\Dropbox\01_norskov\00_PythonModules"
REM setx PYTHONPATH "%PYTHONPATH%;F:\Dropbox\01_norskov\00_git_repos\CatLearn"
REM setx PYTHONPATH "%PYTHONPATH%;F:\Dropbox\01_norskov\00_git_repos\PROJ_IrOx_Active_Learning_OER\python_classes"

setx PYTHONPATH %DROPBOX_DIR%"F:\Dropbox\01_norskov\00_PythonModules;F:\Dropbox\01_norskov\00_git_repos\CatLearn;F:\Dropbox\01_norskov\00_git_repos\PROJ_IrOx_Active_Learning_OER\python_classes"

