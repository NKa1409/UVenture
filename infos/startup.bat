@echo off
SETLOCAL
cd /d "%~dp0"
cd ..

REM Check if Python is installed
where python >nul 2>nul
IF %ERRORLEVEL% NEQ 0 (
    echo Python is not installed. Downloading and installing Python...

    REM Define temporary installer path
    set "PYTHON_INSTALLER=python-installer.exe"

    REM Download the Python installer (you can change version if needed)
    powershell -Command "Invoke-WebRequest -Uri https://www.python.org/ftp/python/3.12.2/python-3.12.2-amd64.exe -OutFile '%PYTHON_INSTALLER%'"

    REM Install Python silently for all users and add to PATH
    start /wait "" "%PYTHON_INSTALLER%" /quiet InstallAllUsers=1 PrependPath=1 Include_test=0

    REM Clean up installer
    del "%PYTHON_INSTALLER%"

    REM Refresh environment variables in current session
    powershell -Command "[System.Environment]::SetEnvironmentVariable('PATH', [System.Environment]::GetEnvironmentVariable('PATH','Machine'), 'Process')"

    REM Verify installation
    where python >nul 2>nul
    IF %ERRORLEVEL% NEQ 0 (
        echo Python installation failed. Exiting...
        exit /b 1
    )
)

REM Check if venv exists
IF NOT EXIST "venv\Scripts\activate.bat" (
    echo Virtual environment not found. Creating...
    python -m venv venv
    call venv\Scripts\activate.bat
    echo Installing dependencies...
    pip install --upgrade pip
    IF EXIST requirements.txt (
        pip install -r requirements.txt
    )
) ELSE (
    echo Virtual environment already exists.
    call venv\Scripts\activate.bat

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
)

echo Running UVenture_web.py...
python UVenture_web.py

pause