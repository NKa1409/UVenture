# UVenture

 **Analyze Full MS / AIF experiments and reconstruct quasi-isolated MS² spectra.**
 
This repository contains the code for a tool developed to process high-resolution Full MS and All-Ion Fragmentation (AIF) mass spectrometry data.  
The program allows researchers to reconstruct quasi-isolated MS² spectra from complex AIF datasets, facilitating deeper structural elucidation and compound identification.

Future versions with enhanced functionality and optimizations are planned and will be published in this repository.

### Features
-  Easy to use (simple user interface and installation process)
-  Identification of sum formulas of precursor ions and related fragment ions.
-  Reconstruction of quasi-isolated MS² spectra from Full MS / AIF data.
-  Flexible parameter settings for different instruments and experimental setups.
-  Spectra, chromatograms and logfiles for every compound.
-  Export of a summary file for easy downstream analysis.

### Screenshots

| Isolated MS² spectrum of Acetylsalicylic acid as captured with a data-dependent Measurement | Reconstructed MS² spectrum of Acetylsalicylic acid as captured by a data-independent (AIF) measurement and processed with the tool |
|:------------------:|:--------------------------:|
| ![Raw Spectra Screenshot](images/isolatedMSMS.png) | ![Reconstructed Spectra Screenshot](images/quasi_isolatedMSMS.png) |


![Raw Spectra Screenshot](images/WebsiteDescription.png)

Image of the Homepage as can be seen when accessing *IP:5000*. To start a new analysis, a *.mzML file must be uploaded to the server and can then be selected in the dropdown menu. Select the desired type of analysis and click "Upload to server".

---

## Installation
0. Make sure to have python installed on your system / on the server (tested with version: 3.12.0)
   
2. Clone or download & unzip the repository and then go (*cd*) to the new folder:
   ```bash
   cd path/of/the/new/folder # navigate to the newly downloaded folder
   
3. Create and activate a virtual environment
   ```bash
   # Create virtual environment (venv)
   python -m venv path/to/venv

   # Activate venv
   source path/to/venv/bin/activate # on Linux
   path\to\venv\Scripts\activate # on Windows

4. Install the required Python packages:
   ```bash
   pip install -r requirements.txt

5. Copy the formula cache to the appropriate location. This can either be done with the command line or via the desktop GUI. The destination folder where the formula cache should be saved can be checked and changed in the settings.
   

## Usage


To run the program:
1. Activate the venv.
2. Start the *UVenture_web.py* python script.
   
   ```bash
   source path/to/venv/bin/activate # Linux
   path\to\venv\Scripts\activate # Windows
   
   python UVenture_web.py

3. Open the browser and go to *IP:5000*. Where IP stands for the IP address of the PC/Server where the script is executed. E.g. *127.0.0.1:5000* or *192.168.178.2:5000*.
   This allows remote access to the program from PCs connected to the same local network as the server.

---

## Notes & Ideas for upcoming Versions
Future updates are planned to include:

- [ ] Convert proprietary raw data formats in .mzML files.
- [ ] Improve batch processing capabilities. Allow to schedule multiple tasks/analyses.
- [ ] Automatic peak detection.
- [ ]    Smoothing / Background subtraction of XIC
- [ ]    Detection of peaks (RT, FWHM, exact mass, mass deviation)  
- [ ]    Detection of the type of ion. Either precursor ion or fragment ion.
- [ ]    Deconvoluion of overlapping peaks.
- [ ] Allow for GC/EI-HRMS data analysis with fragment annotation and precursor ion identification.

## Contributing
Contributions and recommendations for upcoming versions are welcome!

## Citation
If you use this software for your research, please cite:
> **Niklas Karbach, Thorsten Hoffmann**, *TITLE OF THE PUBLICATION*, *Journal Name*, Year 2025. DOI: `XXXXX`


## License
This project is licensed under the MIT License.
See the LICENSE file for more details.

## Contact
Niklas Karbach: n.karbach@uni-mainz.de
Thorsten Hoffmann: t.hoffmann@uni-mainz.de
