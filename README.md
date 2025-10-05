<p align="center">
  <img src="infos/images/uventure_icon.png" alt="UVenture" width="150" height="150">
</p>

## Table of Contents

- **Introduction**
- **Installation & Usage**
- **Screenshots**
- **Notes & Ideas for upcoming Versions**


# UVenture

 **Analyze Full MS / AIF experiments and reconstruct quasi-isolated MS² spectra.**
 
This repository contains the code for a tool developed to process high-resolution Full MS and All-Ion Fragmentation (AIF) mass spectrometry data.  
The program allows researchers to reconstruct quasi-isolated MS² spectra from complex AIF datasets. By nature of the data-independent acquisition (DIA) that AIF is, more information are contained in the resulting measurement file compared to more traditional data-dependent acquisition (DDA) experiments. Comparison with isolated dd-MS2 experiments revealed that the tool delivers comparable results to traditional/manual analysis. 
In contrast to other programs, UVenture creates a local webserver that is reachable by all members of a team via a normal internet browser. The website provides a clean GUI for all interaction with the program.

Future versions with enhanced functionality and optimizations are planned and will be published in this repository.

### Features
-  Automatic installer script for Windows OS
-  Easy-to-use (simple user interface on a local website)
-  Flexible parameter settings for different instruments and experimental setups.
-  Identification of sum formulas of precursor ions and related fragment ions.
-  Reconstruction of quasi-isolated MS² spectra from Full MS / AIF data. 
-  Spectra, chromatograms and logfiles for every compound.
-  Export of a summary file for easy downstream analysis.

---

<div style="page-break-after: always;"></div>

# Installation & Usage

You can either follow the detailled manual installation instructions as shown in infos/installation_guide.md, or use the automatic installer script that is provided in `start_on_windows.bat` (for Windows). For new users it is strongly recommended to use the automatic installer script.

<br>

## Automatic installation (**recommended**)

**Step 1)** Download and extract the program folder on your computer.  After clicking the button shown below, choose `Download ZIP`. The download will start automatically. After the download has finished, unzip the folder.
<p align="center">
  <img src="infos/images/DownloadzipImage.png" alt="Download ZIP Folder" width="300">
</p>

**Step 2)** Double-click the `start_on_windows.bat` (for Windows) script. The installation will be handled automatically - please be patient, installation of the required packages might take a **long** time. Once finished, a new browser tab with the tool will open.

* To restart the program (e.g. after a reboot), just click on the same `start_on_windows.bat` that was used in step 1.

> If the automatic installer script gets flagged by Windows Defender you can click 'More Info' > 'Run anyway'.
> Otherwise you can also install the script manually.

<br><br>

## Manual installation (see infos/installation_guide.md)

**Step 1)** Install Python3 (3.12.0).  
**Step 2)** Download and extract this repository on your computer.  
**Step 3)** Create the virtual environment with the command prompt.  
**Step 4)** Activate the virtual environment.  
**Step 5)** Install the dependencies from requirements.txt.  
**Step 6)** Start UVenture_web.py  

* To restart the program (e.g. after a reboot), redo steps 3.1, 3.2, 4 and 6 from infos/installation_guide.md.

<br><br>

## Usage

1. Open the browser and go to *IP:5000*. Where IP stands for the IP address of the PC/Server where the script is executed. E.g. *127.0.0.1:5000* or *192.168.178.2:5000*.
   This allows remote access to the program from PCs connected to the same local network as the server.

2. Remember to eventually update the provided formula cache (Formula Predictions) to support a larger m/z range. 

---

<div style="page-break-after: always;"></div>

# Screenshots

![Raw Spectra Screenshot](infos/images/WebsiteDescription.png)

Image of the Homepage as can be seen when accessing *IP:5000*. To start a new analysis, a *.mzML file must be uploaded to the server and can then be selected in the dropdown menu. Select the desired type of analysis and click "Upload to server".


| Isolated MS² spectrum of Acetylsalicylic acid as captured with a data-dependent Measurement | Reconstructed MS² spectrum of Acetylsalicylic acid as captured by a data-independent (AIF) measurement and processed with the tool |
|:------------------:|:--------------------------:|
| ![Raw Spectra Screenshot](infos/images/isolatedMSMS.png) | ![Reconstructed Spectra Screenshot](infos/images/quasi_isolatedMSMS.png) |


---

<div style="page-break-after: always;"></div>

# Notes & Ideas for upcoming Versions
Future updates are planned to include:

- [x] Add API endpoints
- [x] Improve batch processing capabilities. Allow to schedule multiple tasks/analyses.
- [x] Bug fix: File browser
- [x] Add estimation for total runtime
- [x] Check RAM usage and prevent from using over 90% of the RAM
- [ ] Add a timeout for processes. If a process takes longer than x minutes (end it / start an additional one / stop it and add peak to peaklist again). 
- [ ] Add file converter to .mzML format (msconvert.exe).
- [ ] Automatic peak detection.
  - [x] Smoothing / Background subtraction of XIC
  - [x] Detection of peaks (RT, FWHM, exact mass, mass deviation)  
  - [ ] Detection of the type of ion. Either precursor ion or fragment ion.
  - [x] Deconvoluion of overlapping peaks.
- [ ] Create a database search tool to handle large numbers of individual raw files.
  - [ ] Detect peaks with a given m/z in every raw file (create XIC). E.g. see if PFOA can be found in the samples.
  - [ ] Calculate sum formula and fragments of every peak for the given m/z value, for every raw file. Are it just isomers or completely different molecules.
  - [ ] Compare the different fragmentation patterns and see if similar patterns can be found across the raw data files in the database. Those might then be the same molecules, although the LC method might have been different.
  - [ ] Give a bar chart for every individual identified compound to directly compare the individual raw files.
- [ ] Create *local* user accounts to allow for settings to be saved per user
  - [ ] Allow user specific file uploads
  - [ ] Allocate x cores / RAM for one user
  - [ ] Allow user specific settings
  - [ ] Access controll for result files
- [ ] Allow for GC/EI-HRMS data analysis with fragment annotation and precursor ion identification.


<div style="page-break-after: always;"></div>

## Contributing
Contributions and recommendations for upcoming versions are welcome!

## Citation
If you use this software for your research, please consider citing:
> **Niklas Karbach, Thorsten Hoffmann**: github.com/NKa1409/UVenture


## License
This project is licensed under the MIT License.
See the LICENSE file for more details.

## Contact
Niklas Karbach: n.karbach@uni-mainz.de
Thorsten Hoffmann: t.hoffmann@uni-mainz.de
