#!/bin/bash

set -e

cd "$(dirname "$0")"
cd ..

function pause_on_error {
    echo
    echo "An error occurred. Press Enter to exit."
    read
    exit 1
}

trap pause_on_error ERR

check_requirements_installed() {
    if [ ! -f "requirements.txt" ]; then
        return 0
    fi

    echo "Checking for required packages..."
    missing=$(pip install --dry-run -r requirements.txt 2>/dev/null | grep -i "would install")

    if [ -n "$missing" ]; then
        echo "Some packages are missing. Installing..."
        pip install --upgrade pip
        pip install -r requirements.txt
    else
        echo "All required packages are already installed."
    fi
}

if [ ! -d "venv" ]; then
    echo "Virtual environment not found. Creating..."
    python3 -m venv venv
    source venv/bin/activate
    echo "Installing dependencies..."
    pip install --upgrade pip
    [ -f "requirements.txt" ] && pip install -r requirements.txt
else
    echo "Virtual environment already exists."
    source venv/bin/activate
    check_requirements_installed
fi

echo "Running main.py..."
python main.py || pause_on_error

echo
read -p "Script finished. Press Enter to close."