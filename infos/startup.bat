@echo off
SETLOCAL
cd /d "%~dp0"
cd ..

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