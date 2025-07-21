#!/bin/bash
set -e
cd "$(dirname "$0")"

echo
echo "===================================================================="
echo " Searching for existing Python installation"
echo "===================================================================="

PYTHON_EXEC=""

# 1. Try local folder ./python312/bin/python3
LOCAL_PYTHON="./python312/bin/python3"
if [ -x "$LOCAL_PYTHON" ]; then
    PYTHON_EXEC="$LOCAL_PYTHON"
    echo "Found local Python: $PYTHON_EXEC"
else
    # 2. Try system Python
    if command -v python3 &>/dev/null; then
        PYTHON_EXEC=$(command -v python3)
        echo "Found system Python: $PYTHON_EXEC"
    else
        echo
        echo "===================================================================="
        echo " Python not found. Installing Python 3.12.2"
        echo "===================================================================="

        echo "You need Python 3.12.2 or newer. Installing locally..."
        sudo apt update
        sudo apt install -y build-essential wget libssl-dev zlib1g-dev \
            libncurses5-dev libncursesw5-dev libreadline-dev libsqlite3-dev \
            libgdbm-dev libdb5.3-dev libbz2-dev libexpat1-dev liblzma-dev tk-dev

        mkdir -p python_src
        cd python_src
        wget https://www.python.org/ftp/python/3.12.2/Python-3.12.2.tgz
        tar -xzf Python-3.12.2.tgz
        cd Python-3.12.2
        ./configure --prefix="$(dirname "$(pwd)")/python312" --enable-optimizations
        make -j$(nproc)
        make install
        cd ../..
        PYTHON_EXEC="./python312/bin/python3"
    fi
fi

# Ensure it's Python 3.12.2 or newer
PYTHON_VERSION=$("$PYTHON_EXEC" -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
REQUIRED_VERSION="3.12.2"
if [[ "$(printf '%s\n' "$REQUIRED_VERSION" "$PYTHON_VERSION" | sort -V | head -n1)" != "$REQUIRED_VERSION" ]]; then
    echo "[ERROR] Python version $REQUIRED_VERSION or higher is required."
    echo "Found version: $PYTHON_VERSION"
    exit 1
else
    echo "Using Python: $PYTHON_EXEC (version $PYTHON_VERSION)"
fi

echo
echo "===================================================================="
echo " Create and check virtual environment"
echo "===================================================================="

VENV_DIR="venv"
VENV_PY="$VENV_DIR/bin/python"

if [ -f "$VENV_DIR/bin/activate" ]; then
    echo "Virtual environment found. Validating..."
    if ! "$VENV_PY" --version &>/dev/null; then
        echo "Virtual environment is broken. Recreating..."
        rm -rf "$VENV_DIR"
    else
        echo "Virtual environment is valid."
    fi
fi

if [ ! -f "$VENV_DIR/bin/activate" ]; then
    echo "Creating virtual environment..."
    "$PYTHON_EXEC" -m venv "$VENV_DIR"
fi

source "$VENV_DIR/bin/activate"

echo
echo "===================================================================="
echo " Install packages"
echo "===================================================================="

if [ -f requirements.txt ]; then
    echo "Checking for required packages..."
    pip install --upgrade pip
    pip install -r requirements.txt
else
    echo "No requirements.txt found."
fi

echo
echo "===================================================================="
echo " Start UVenture"
echo "===================================================================="

echo "Running UVenture_web.py..."

# Parse host/port from settings
webserver_settings_filepath="static/webserver_settings.txt"
host="127.0.0.1"
port="5000"

if [ -f "$webserver_settings_filepath" ]; then
    echo "Webserver settings file found."
    while IFS='=' read -r key value; do
        case "$key" in
        host) host="$value" ;;
        port) port="$value" ;;
        esac
    done <"$webserver_settings_filepath"
else
    echo "Webserver settings file not found. Using default host: $host"
fi

# Replace 0.0.0.0 with localhost for browser
[ "$host" = "0.0.0.0" ] && host="127.0.0.1"

# Run server in background
nohup python UVenture_web.py &

echo "Waiting for server to start on http://$host:$port ..."

wait_seconds=120
waited=0

while ! curl -s --max-time 1 "http://$host:$port" >/dev/null; do
    sleep 1
    waited=$((waited + 1))
    if [ "$waited" -ge "$wait_seconds" ]; then
        echo "Timeout reached. Proceeding to open browser anyway..."
        break
    fi
done

echo "Server is up! Launching browser..."
xdg-open "http://$host:$port" || echo "Please open your browser and visit: http://$host:$port"
