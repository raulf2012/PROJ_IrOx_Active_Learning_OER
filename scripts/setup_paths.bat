@echo off

REM set DROPBOX_DIR=C:\Users\raulf2012\Dropbox
set DROPBOX_DIR=F:\Dropbox

setx PROJ_irox %DROPBOX_DIR%"\01_norskov\00_git_repos\PROJ_IrOx_Active_Learning_OER"
setx PROJ_irox_2 %DROPBOX_DIR%"\01_norskov\01_projects\04_irox_oer_orr"
setx PROJ_DATA %DROPBOX_DIR%"\01_norskov\PROJECT_DATA"



set TEMP0=%DROPBOX_DIR%\01_norskov\00_PythonModules
set TEMP1=%DROPBOX_DIR%\01_norskov\00_git_repos\CatLearn
set TEMP2=%DROPBOX_DIR%\01_norskov\00_git_repos\PROJ_IrOx_Active_Learning_OER\python_classes
set TEMP3=%DROPBOX_DIR%\01_norskov\00_git_repos\ase

setx PYTHONPATH %TEMP0%;%TEMP1%;%TEMP2%;%TEMP3%







REM set TEMP_VAR=F:\Dropbox
REM setx PROJ_irox "F:\Dropbox\01_norskov\00_git_repos\PROJ_IrOx_Active_Learning_OER"
REM setx PROJ_DATA "F:\Dropbox\01_norskov\PROJECT_DATA"
REM setx PYTHONPATH "F:\Dropbox\01_norskov\00_PythonModules;F:\Dropbox\01_norskov\00_git_repos\CatLearn"
REM setx PYTHONPATH "%PYTHONPATH%;F:\Dropbox\01_norskov\00_PythonModules"
REM setx PYTHONPATH "%PYTHONPATH%;F:\Dropbox\01_norskov\00_git_repos\CatLearn"
REM setx PYTHONPATH "%PYTHONPATH%;F:\Dropbox\01_norskov\00_git_repos\PROJ_IrOx_Active_Learning_OER\python_classes"
REM setx PYTHONPATH %DROPBOX_DIR%"\01_norskov\00_PythonModules;F:\Dropbox\01_norskov\00_git_repos\CatLearn;F:\Dropbox\01_norskov\00_git_repos\PROJ_IrOx_Active_Learning_OER\python_classes"
