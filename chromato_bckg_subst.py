import ast
import copy
import datetime
import math
import pathlib
try:
    import cPickle as pickle
except ModuleNotFoundError:
    print("cPickle not found. Using pickle instead.")
    import pickle
import sys
import traceback
import PIL
import fontTools.t1Lib
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
import matplotlib
import UVenture
import matplotlib.gridspec
import numpy as np
import psutil
import scipy
import MS_functions
import pyteomics
import flask
import time
import threading
import werkzeug
import os
from pyteomics import mzml
import similaritymeasures

from plotly.subplots import make_subplots
import plotly.graph_objects as go

import shutil


def get_mode_of_spec(filter_string):
    if " d " in filter_string and "@hcd" in filter_string:
        ms_ms_masses = []
        filter_parsed = filter_string.split(" ")
        filter_parsed = [x for x in filter_parsed if "hcd" in x]
        for element in filter_parsed:
            ms_ms_masses.append(round(float(element.split("@")[0]), 2))
        return "MS/MS"
    elif " d " not in filter_string and "hcd" not in filter_string:
        return "Full scan"
    elif " d " not in filter_string and "hcd" in filter_string:
        return "AIF"







if __name__ == "__main__":
    print("Start")
    blank_filename = "V://005 Mitarbeiter-aktuell//Karbach//Backups//MessdatenFrankfurtPamChamber//05Dec2024_C18//05Dec2024_PAMPain3_neg_20241008BLANK.mzML"
    data_filename = "V://005 Mitarbeiter-aktuell//Karbach//Backups//MessdatenFrankfurtPamChamber//05Dec2024_C18//05Dec2024_PAMPain2_neg_20241008_ISP_OH_O3_NO.mzML"

    blank_ms_file = UVenture.MS_File(blank_filename)
    #data_ms_file = UVenture.MS_File(data_filename)
    with open('datamsfile.pkl', 'rb') as inp:
        data_ms_file = pickle.load(inp)
    #bg_subst_data_ms_file = copy.deepcopy(data_ms_file)

    print(data_ms_file.rawdata[0])
    max_deviation = 15
    print("start2")
    print(len(blank_ms_file.rawdata))
    print(len(data_ms_file.rawdata))
    for i in range(len(data_ms_file.rawdata)):
        filter_name = get_mode_of_spec(data_ms_file.rawdata[i]["scanList"]["scan"][0]["filter string"])
        print(filter_name)
        spec_data = UVenture.Spec(ms_file=data_ms_file, index=i, spec_requested_filter_mode=filter_name)
        rt = spec_data.rt
        print(rt)
        best_rt_at_blank = min(blank_ms_file.rt_list, key=lambda x: abs(rt - x))
        print(best_rt_at_blank)
        best_index_at_blank = blank_ms_file.rt_list.index(best_rt_at_blank)
        print(best_index_at_blank)
        spec_blank = UVenture.Spec(ms_file=blank_ms_file, index=best_index_at_blank, spec_requested_filter_mode=filter_name)
        print("specdatafinished")
        bsubst_mass_int_dict = copy.deepcopy(spec_data.summarized_mass_intensity_dict)
        print(len(bsubst_mass_int_dict))
        for k in range(len(spec_blank.summarized_masses)):
            b_mass = spec_blank.summarized_masses[k]
            b_int = spec_blank.summarized_intensities[k]
            best_mass = min(list(bsubst_mass_int_dict.keys()), key=lambda x: abs(b_mass - x))
            #print("blank mass: " + str(b_mass) + "  data mass: " + str(best_mass) + "  deviation: " + str(( (abs(best_mass - b_mass)/b_mass) * 1000000 )))
            if ( (abs(best_mass - b_mass)/b_mass) * 1000000 ) <= max_deviation/2:
                bsubst_mass_int_dict[best_mass] = bsubst_mass_int_dict[best_mass] - b_int
                if bsubst_mass_int_dict[best_mass] <= 0:
                    del bsubst_mass_int_dict[best_mass]
        print(len(bsubst_mass_int_dict))
        data_ms_file.rawdata[i]["m/z array"] = np.array(list(bsubst_mass_int_dict.keys()))
        data_ms_file.rawdata[i]["intensity array"] = np.array(list(bsubst_mass_int_dict.values()))
        del spec_blank
        del spec_data
        print("progress: " + str(i) + "/" + str(len(data_ms_file.rawdata)))


    with open('datamsfile.pkl', 'wb') as outp:
        pickle.dump(data_ms_file, outp, pickle.HIGHEST_PROTOCOL)




