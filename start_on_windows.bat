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

:: 1. Try global PATH
python --version >nul 2>nul
IF %ERRORLEVEL% EQU 0 (
    for /f "delims=" %%P in ('where python') do (
        SET "PYTHON_EXEC=%%P"
        goto :found_python
    )
)

:: 2. Try default user install location
echo Python not found in PATH. Looking for default user installation...
FOR %%P IN ("%LocalAppData%\Programs\Python\Python3*\python.exe") DO (
    IF EXIST %%P (
        SET "PYTHON_EXEC=%%P"
        goto :found_python
    )
)

:: 3. Not found — install Python now
echo Could not locate Python. Installing Python for current user...
set "PYTHON_INSTALLER=python-installer.exe"
powershell -Command "Invoke-WebRequest -Uri https://www.python.org/ftp/python/3.12.2/python-3.12.2-amd64.exe -OutFile '%PYTHON_INSTALLER%'"
start /wait "" "%PYTHON_INSTALLER%" /quiet InstallAllUsers=0 PrependPath=1 Include_test=0
del "%PYTHON_INSTALLER%"

:: Try to find Python again after install
FOR %%P IN ("%LocalAppData%\Programs\Python\Python3*\python.exe") DO (
    IF EXIST %%P (
        SET "PYTHON_EXEC=%%P"
        goto :found_python
    )
)

:found_python
echo Using Python at: %PYTHON_EXEC%

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
    %PYTHON_EXEC% -m venv "%VENV_DIR%"
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

echo Waiting for server to start on http://localhost:5000 ...

set /a wait_seconds=30
set /a waited=0

:wait_for_server
powershell -Command ^
  "$r=0; try { $r=(Invoke-WebRequest -UseBasicParsing -Uri 'http://127.0.0.1:5000' -TimeoutSec 1).StatusCode } catch {}; if ($r -eq 200) { exit 0 } else { exit 1 }"
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
start "" http://127.0.0.1:5000

pause
