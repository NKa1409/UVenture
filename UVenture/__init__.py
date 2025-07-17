# src/__init__.py

"""
This package contains custom modules for the UVenture project.
"""

try:
    import pyteomics
except ImportError:
    print("pyteomics is not installed. Please install it to use UVenture functionalities."
          " You can install it using 'pip install pyteomics'.")
    print("Or use the provided requirements.txt file to install all dependencies needed for the UVenture Package.")
    print("Restart the application after installation.")
    input("Press Enter to continue...")

try:
    import matplotlib
    matplotlib.use('Agg')  # Use a non-interactive backend for matplotlib
except ImportError:
    print("matplotlib is not installed. Please install it to use UVenture functionalities."
          " You can install it using 'pip install matplotlib'.")
    print("Or use the provided requirements.txt file to install all dependencies needed for the UVenture Package.")
    print("Restart the application after installation.")
    input("Press Enter to continue...")

try:
    import similaritymeasures
except ImportError:
    print("similaritymeasures is not installed. Please install it to use UVenture functionalities."
          " You can install it using 'pip install similaritymeasures'.")
    print("Or use the provided requirements.txt file to install all dependencies needed for the UVenture Package.")
    print("Restart the application after installation.")
    input("Press Enter to continue...")

try:
    import scipy
except ImportError:
    print("scipy is not installed. Please install it to use UVenture functionalities."
          " You can install it using 'pip install scipy'.")
    print("Or use the provided requirements.txt file to install all dependencies needed for the UVenture Package.")
    print("Restart the application after installation.")
    input("Press Enter to continue...")

try:
    import pandas
except ImportError:
    print("pandas is not installed. Please install it to use UVenture functionalities."
          " You can install it using 'pip install pandas'.")
    print("Or use the provided requirements.txt file to install all dependencies needed for the UVenture Package.")
    print("Restart the application after installation.")
    input("Press Enter to continue...")

try:
    import numpy
except ImportError:
    print("numpy is not installed. Please install it to use UVenture functionalities."
          " You can install it using 'pip install numpy'.")
    print("Or use the provided requirements.txt file to install all dependencies needed for the UVenture Package.")
    print("Restart the application after installation.")
    input("Press Enter to continue...")





print("UVenture package loaded successfully.")





