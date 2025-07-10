# Python Script Setup & Execution Guide

This guide provides detailed instructions to set up and run a Python script from a ZIP archive. It includes:

- Installing Python
- Extracting the ZIP archive
- Creating and activating a virtual environment
- Installing dependencies
- Running the Python script

---

## 1. Install Python

### 1.1 Download Python

Go to the official Python download page: https://www.python.org/downloads/
Download a supported version of Python 3 (3.12)

#### Windows:
1. Run the downloaded `.exe` file.
2. **Check** the box that says: `Add Python 3.x to PATH`
3. Click **Customize installation**, enable all options, and proceed.
4. Select **Install for all users** on the next screen.
5. Click **Install**.

#### macOS:
1. Open the downloaded `.pkg` installer.
2. Follow the installation instructions.
3. Python will be installed in `/usr/local/bin/python3`.

#### Linux (Debian/Ubuntu-based):
Open a terminal and run (you do not have to download the file from python.org):

```bash
sudo apt update
sudo apt install python3 python3-venv python3-pip -y


## 2. Extract the ZIP Archive

1) Locate and download the ZIP archive containing the project files from this repo.
2) Right-click the ZIP file and select Extract All or use an extraction tool.
3) Note the path of the extracted folder (e.g., C:\Users\YourName\Downloads\Project).


## 3. Create a Virtual Environment
### 3.1 Open Terminal / Command Prompt
Windows: Press Win + R, type cmd, and press Enter.

macOS/Linux: Open Terminal.

### 3.2 Navigate to the Project Folder
Replace with the actual path to your extracted project folder:
cd path/to/project-folder

Examples:
Windows:
```bash
cd "C:\Users\YourName\Downloads\Project"

macOS/Linux:
```bash
cd ~/Downloads/Project

