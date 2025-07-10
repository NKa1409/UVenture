# Python Script Setup & Execution Guide

This guide provides detailed instructions to set up and run a Python script from a ZIP archive. It includes:

- Installing Python
- Extracting the ZIP archive
- Creating and activating a virtual environment
- Installing dependencies
- Running the Python script


## If the script has already run once, repeat steps 3.1, 3.2, 4 and 6 only.

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
```


## 2. Extract the ZIP Archive

1) Locate and download the ZIP archive containing the project files from this repo.
2) Right-click the ZIP file and select Extract All or use an extraction tool.
3) Note the path of the extracted folder (e.g., `C:\Users\YourName\Downloads\UVenture-main`).


## 3. Create a Virtual Environment
### 3.1 Open Terminal / Command Prompt
Windows: Press `Win + R`, type `cmd`, and press Enter.

macOS/Linux: Open Terminal.

### 3.2 Navigate to the Project Folder
Replace with the actual path to your extracted project folder:
```bash
cd path/to/project-folder
```

Examples:
Windows:
```bash
cd "C:\Users\YourName\Downloads\UVenture-main"
```

macOS/Linux:
```bash
cd ~/Downloads/UVenture-main
```

## 3.3 Create the Virtual Environment
```bash
python -m venv venv
```
This will create a folder named venv inside your project directory.

## 4. Activate the Virtual Environment
Windows:
```bash
venv\Scripts\activate
```
macOS/Linux:
```bash
source venv/bin/activate
```
You will see the terminal prompt change to indicate the environment is active, e.g., (venv).

## 5. Install Dependencies
Install the required requirements (in requirements.txt):

```bash
pip install -r requirements.txt
```

## 6. Run the Script
```bash
python UVenture_web.py
```
If everything is set up correctly, the script should execute without errors.

## 7. Deactivate the Virtual Environment (Optional)
After you're finished and do not need the program anymore, you can deactivate the environment by running:
```bash
deactivate
```


## Troubleshooting

| Issue                        | Solution                                       |
| ---------------------------- | ---------------------------------------------- |
| `'python' is not recognized` | Use `py` instead of `python`, or fix PATH      |
| `ModuleNotFoundError`        | Run `pip install package-name`                 |
| Permissions errors           | Run terminal as Administrator / use `sudo`     |
| Wrong Python version         | Check with `python --version` or use `python3` |

