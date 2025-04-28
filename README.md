# UVenture
 **Analyze Full MS / AIF experiments and reconstruct quasi-isolated MS² spectra.**

---

## Introduction

This repository contains the code for a tool developed to process Full MS and All-Ion Fragmentation (AIF) mass spectrometry data.  
The program allows researchers to reconstruct quasi-isolated MS² spectra from complex AIF datasets, facilitating deeper structural elucidation and compound identification.

Future versions with enhanced functionality and optimizations are planned and will be published in this repository.

---

## Features

- ✅ Reconstruction of quasi-isolated MS² spectra from Full MS / AIF data
- ✅ Automated peak picking and deconvolution
- ✅ Flexible parameter settings for different instruments and experimental setups
- ✅ User-friendly graphical outputs
- ✅ Export of reconstructed spectra for downstream analysis

---

## Installation

1. Clone or download & unzip the repository and then go *cd* to the new folder:
   ```bash
   cd path\of\the\new\folder # navigate to the newly downloaded folder
   
2. Create and activate a virtual environment
   ```bash
   # Create virtual environment (venv)
   python -m venv path/to/venv

   # Activate venv
   source path/to/venv/bin/activate # on Linux
   path\to\venv\Scripts\activate # on Windows

2. Install the required Python packages:
   ```bash
   pip install -r requirements.txt

## Usage


To run the program:
1. Activate the venv.
2. Start the *UVenture_web.py* python script.
   
   ```bash
   *(Linux)*
   source path/to/venv/bin/activate
   *(Windows)*
   path\to\venv\Scripts\activate
   
   python UVenture_web.py


You can specify your own input data and parameter files by modifying the configuration options inside the program or via CLI arguments (depending on the version).


## Screenshots

| Raw Data Analysis | Reconstructed MS² Spectra |
|:------------------:|:--------------------------:|
| ![Raw Spectra Screenshot](images/raw_spectra.png) | ![Reconstructed Spectra Screenshot](images/reconstructed_spectra.png) |


## Notes for Upcoming Versions
Future updates are planned to include:

🚀 Full GUI integration (PyQt5 / Tkinter)
🚀 Support for additional mass spectrometry file formats (e.g., mzML, mzXML)
🚀 Advanced noise filtering algorithms
🚀 Improved deconvolution algorithms for overlapping fragments
🚀 Batch processing capabilities
🚀 Direct export to popular libraries like GNPS or METLIN
🚀 Extensive unit tests and CI/CD integration

## Contributing
Contributions are welcome! If you would like to:

Report a bug

Request a feature

Submit a pull request

Please follow the Contribution Guidelines.

## Citation
If you use this software for your research, please cite:
> **Niklas Karbach, Thorsten Hoffmann**, *TITLE OF THE PUBLICATION*, *Journal Name*, Year 2025. DOI: `XXXXX`


## License
This project is licensed under the MIT License.
See the LICENSE file for more details.

## Contact
