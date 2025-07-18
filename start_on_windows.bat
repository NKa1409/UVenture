@echo off
SETLOCAL ENABLEDELAYEDEXPANSION
cd /d "%~dp0"


:: Check for admin rights
net session >nul 2>&1
IF %ERRORLEVEL% NEQ 0 (
    echo.
    echo [WARNING] You are not running this script with administrator privileges.
    echo If Python installation fails or the app does not start correctly,
    echo try running this script again as administrator.
    echo.
)

:: === Find or install Python ===
SET "PYTHON_EXEC="

:: 1. Try local folder python312 in script directory
echo Python not found in default user installation path. Checking in UVenture folder...
SET "LOCAL_PYTHON=%~dp0python312\python.exe"
IF EXIST "%LOCAL_PYTHON%" (
    SET "PYTHON_EXEC=%LOCAL_PYTHON%"
    goto :found_python
)

:: 2. Try global PATH
python --version >nul 2>nul
IF %ERRORLEVEL% EQU 0 (
    for /f "delims=" %%P in ('where python') do (
        SET "PYTHON_EXEC=%%P"
        goto :found_python
    )
)

:: 3. Try default user install location
echo Python not found in PATH. Looking for default user installation...
FOR %%P IN ("%LocalAppData%\Programs\Python\Python3*\python.exe") DO (
    IF EXIST %%P (
        SET "PYTHON_EXEC=%%P"
        goto :found_python
    )
)

:: 4. Not found — install Python now
echo Could not locate any version of Python on your system. Downloading the installer for Python 3.12.0.
echo This will not affect system-wide Python installations.
echo Please wait while the file downloads and installs Python...
set "PYTHON_INSTALLER=python-installer.exe"
powershell -Command "Invoke-WebRequest -Uri https://www.python.org/ftp/python/3.12.0/python-3.12.0-amd64.exe -OutFile '%PYTHON_INSTALLER%'"
:: Get path to current script directory
set "SCRIPT_DIR=%~dp0"
:: Set Python install directory to a subfolder "python312" inside the script's directory
set "PYTHON_TARGET=%SCRIPT_DIR%python312"
:: Ensure the folder exists
mkdir "%PYTHON_TARGET%"
start /wait "" "%PYTHON_INSTALLER%" /quiet InstallAllUsers=0 PrependPath=0 Include_test=0 TargetDir="%PYTHON_TARGET%"
del "%PYTHON_INSTALLER%"


SET "PYTHON_EXEC=%PYTHON_TARGET%\python.exe"

:found_python
echo Using Python at: %PYTHON_EXEC%

:: Check if Python is version 3.12.2 or higher
FOR /F "tokens=2 delims=." %%A IN ('"%PYTHON_EXEC%" --version 2^>nul') DO (
    IF %%A LSS 12 (
        echo [ERROR] Python version 3.12.2 or higher is required.
        echo Please install the correct version of Python and try again.
        exit /b 1
    ) ELSE (
        echo Python version is sufficient: %%A
    )
)

:: Print the python version
"%PYTHON_EXEC%" --version


REM Now use "%PYTHON_EXEC%" instead of "python" for venv, pip, running scripts, etc.


REM === Validate or recreate venv ===
SET "VENV_DIR=venv"
SET "VENV_PY=%VENV_DIR%\Scripts\python.exe"

IF EXIST "%VENV_DIR%\Scripts\activate.bat" (
    echo Virtual environment found. Validating...

    "%VENV_PY%" --version >nul 2>nul
    IF %ERRORLEVEL% NEQ 0 (
        echo Virtual environment is broken. Recreating...
        rmdir /s /q "%VENV_DIR%"
    ) ELSE (
        echo Virtual environment is valid.
    )
)

REM === Create venv if it doesn't exist ===
IF NOT EXIST "%VENV_DIR%\Scripts\activate.bat" (
    echo Creating virtual environment...
    "%PYTHON_EXEC%" -m venv "%VENV_DIR%"
)

REM === Activate venv ===
call "%VENV_DIR%\Scripts\activate.bat"

REM === Install packages if needed ===
IF EXIST requirements.txt (
    echo Checking for required packages...
    pip install --dry-run -r requirements.txt > tmp_check.txt 2>nul
    findstr /I "would install" tmp_check.txt >nul
    IF %ERRORLEVEL%==0 (
        echo Some packages are missing. Installing...
        pip install --upgrade pip
        pip install -r requirements.txt
    ) ELSE (
        echo All required packages are already installed.
    )
    del tmp_check.txt
)

echo Running UVenture_web.py...
start "" /B python UVenture_web.py

:: === Get the configuration file from static/webserver_settings.txt and extract the host and port. ===
set "webserver_settings_filepath=static\webserver_settings.txt"
IF NOT EXIST "%webserver_settings_filepath%" (
    echo Webserver settings file not found. Default settings will be used.
    echo host=127.0.0.1
    echo port=5000
) ELSE (
    echo Webserver settings file found.
)
:: Read the settings from the file
set "host="
set "port="

for /f "usebackq tokens=1,2 delims==" %%A in ("%webserver_settings_filepath%") do (
    if "%%A"=="host" set "host=%%B"
    if "%%A"=="port" set "port=%%B"
)
IF NOT DEFINED host (
    set "host=127.0.0.1"
    set "port=5000"
    echo Using default host: %host%
) ELSE (
    echo Using host from settings: %host%
)
:: If host is 0.0.0.0, change it to 127.0.0.1
IF "%host%"=="0.0.0.0" (
    set "host=127.0.0.1"
)

echo Waiting for server to start on http://%host%:%port% ...

set /a wait_seconds=120
set /a waited=0

:wait_for_server
powershell -Command ^
  "$r=0; try { $r=(Invoke-WebRequest -UseBasicParsing -Uri 'http://%host%:%port%' -TimeoutSec 1).StatusCode } catch {}; if ($r -eq 200) { exit 0 } else { exit 1 }"
IF %ERRORLEVEL%==0 (
    echo Server is up!
    goto :launch_browser
)

:: Wait 1 second and increment counter
timeout /t 1 >nul
set /a waited+=1

IF %waited% LSS %wait_seconds% (
    goto :wait_for_server
) ELSE (
    echo Timeout reached. Proceeding to open browser anyway...
)

:launch_browser
start "" http://%host%:%port%

pause
