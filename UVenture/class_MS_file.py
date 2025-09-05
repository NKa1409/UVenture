import copy
import datetime
import os
import re
import traceback
import numpy as np
from pyteomics import mzml

import UVenture.MS_functions as MS_functions
import UVenture.plotting as plotting


class MS_File:
    def __init__(self, filename=None, **kwargs):
        self.debug_output = True
        self.do_bckg_subtraction = False
        if filename == None:
            print("No filename provided. Cannot read MS file. Returning...")
            return
        else:
            if "parentfolder_msfile" in kwargs: 
                parentfolder = kwargs["parentfolder_msfile"]
            else: 
                parentfolder = os.path.normpath( str(".".join(filename.split(".")[:-1])))
            default_kwargs = {"parentfolder_msfile": parentfolder,
                              "logfile_filepath": os.path.join(parentfolder, "MSfile_logfile.txt"),
                              "msfile_raw_file_retention_time_unit": "sec",
                              "create_2d_spec_of_ms_file": False}
            self.kwargs = {**default_kwargs, **kwargs}
            os.makedirs(self.kwargs["parentfolder_msfile"], exist_ok=True)

            self.get_2d_spec = self.kwargs.get("create_2d_spec_of_ms_file", False)

            self.parentfolder = self.kwargs["parentfolder_msfile"]
            self.filename = filename
            self.file = mzml.read(self.filename)
            self.rawdata = list(self.file)


            self.rename_rawdata_keys()
            print("Renamed rawdata keys: " + str(list(self.rawdata[-1].keys())))

            self.rename_scanlist_keys()
            print("Renamed scanList keys: " + str(list(self.rawdata[-1]["scanList"].keys())))

            # Calculate some basic properties of the MS file
            self.rt_list = [element["scanList"]["scan"][0]["scan time"] for element in self.rawdata]
            self.method_duration = self.rt_list[-1] - self.rt_list[0]
            self.rt_range = [min(self.rt_list), max(self.rt_list)]
            # Now extract the filter strings from the rawdata
            self.all_filters = []
            for index in range(len(self.rawdata)):
                self.all_filters.append(self.rawdata[index]["scanList"]["scan"][0]["filter string"])
            self.mz_range = [ float( self.all_filters[0].split("[")[1].split("-")[0] ), float(self.all_filters[0].split("[")[1].split("-")[1].replace("]", "")) ]

            # Check if retention time unit is seconds or minutes and convert to seconds if necessary
            if self.kwargs["msfile_raw_file_retention_time_unit"] in ["min", "minutes", "minute", "mins", "m"]:
                print("Converting retention time from minutes to seconds...")
                for i in range(len(self.rawdata)):
                    self.rawdata[i]["scanList"]["scan"][0]["scan time"] = self.rawdata[i]["scanList"]["scan"][0]["scan time"] * 60
                self.rt_list = [rt * 60 for rt in self.rt_list]
                self.method_duration = self.method_duration * 60
            
            self.save_ms_file_log_entry("INFO:\t" + "Reading MS file: " + str(self.filename))
            self.save_ms_file_log_entry("INFO:\t" + "Method duration: " + str(self.method_duration))
            self.save_ms_file_log_entry("INFO:\t" + "Number of spectra: " + str(len(self.rawdata)))
            self.save_ms_file_log_entry("INFO:\t" + "RT range: " + str(self.rt_range))
            self.save_ms_file_log_entry("INFO:\t" + "m/z range: " + str(self.mz_range))

            
            # Get the tic from each scan
            self.tic = []
            for index in range(len(self.rawdata)):
                self.tic.append(self.rawdata[index]["total ion current"])
            self.save_ms_file_log_entry("INFO:\t" + "Summed TIC: " + str(sum(self.tic)))
            
            # Check if the last scan has a filter string if not, try to construct it
            try:
                filter_string = self.rawdata[-1]["scanList"]["scan"][0]["filter string"]
            except Exception as e:
                print("No filter string found in the last scan. Trying to construct it...")
                for index in range(len(self.rawdata)):
                    self.rawdata[index]["scanList"]["scan"][0]["filter string"] = self.construct_filter_string(index)
            
            
            # Determine available modes based on filter strings
            self.available_modes = []
            for filter_string in self.all_filters:
                if " d " in filter_string and "@hcd" in filter_string:
                    self.available_modes.append("MS/MS")
                elif " d " not in filter_string and "hcd" not in filter_string:
                    self.available_modes.append("Full scan")
                elif " d " not in filter_string and "hcd" in filter_string:
                    self.available_modes.append("AIF")
            self.all_modes = copy.deepcopy(self.available_modes)
            self.available_modes = list(set(self.available_modes))
            self.save_ms_file_log_entry("INFO:\t" + "Available modes: " + str(self.available_modes))

            # Calculate a background spectrum for AIF and MS1 spectra 
            self.aif_background_spectrum = None
            self.ms1_background_spectrum = None
            self.get_background_spectra(background_range=(3, 10), bckg_m_dev=0.0001, blank_multiplicator=3)
            self.save_ms_file_log_entry("INFO:\t" + "Summed MS1 background signals: " + str(sum(list(self.ms1_background_spectrum.values()))))
            self.save_ms_file_log_entry("INFO:\t" + "Summed AIF background signals: " + str(sum(list(self.aif_background_spectrum.values()))))
            
            if self.debug_output: print("Starting to bg substract the data...")
            starttime = datetime.datetime.now()
            if self.get_2d_spec:
                save_filename = os.path.join(self.parentfolder, "2Dspec.tiff")
                os.makedirs(os.path.dirname(save_filename), exist_ok=True)
                if not os.path.exists(save_filename):
                    plotting.create_2d_massspec_plot(ms_file_object=self, filter_mode="Full scan", save=save_filename, max_dim=5000)
                    self.save_ms_file_log_entry("INFO:\t" + "2D spectrum saved to: " + str(save_filename))
            if self.do_bckg_subtraction:
                if self.debug_output: print("Doing background subtraction...")
                self.save_ms_file_log_entry("INFO:\t" + "Doing background subtraction...")
                self.do_background_subtraction(background_range=(3, 10), multiplicator=3, bckg_m_dev=0.00001)
                self.save_ms_file_log_entry("INFO:\t" + "Finished background subtraction...")
                if self.debug_output: print("Finished background subtraction in: " + str(datetime.datetime.now() - starttime)) 
    
    def get_xic(self, mass, mass_deviation, requested_filter_mode="Full scan"):
        # The mass deviation is defined as the requested mass +1x the mass deviation and -1x the mass deviation.
        # The mass deviation is given in Da (Dalton), and NOT in PPM!!
        # If the requested mass is 1000 and the mass deviation is 0.005, the range is 999.995 to 1000.005.
        valid_modes = {"Full scan", "AIF", "MS/MS"}
        if requested_filter_mode not in valid_modes:
            requested_filter_mode = "Full scan"

        rt_list = []
        original_index_list = []
        intensity_list = []
        for idx, entry in enumerate(self.rawdata):
            
            filter = entry["scanList"]["scan"][0]["filter string"]
            # Determine filter_mode quickly
            if " d " in filter and "@hcd" in filter:
                filter_mode = "MS/MS"
            elif " d " not in filter and "hcd" in filter:
                filter_mode = "AIF"
            elif " d " not in filter and "@hcd" not in filter:
                filter_mode = "Full scan"
            else:
                filter_mode = "Unknown"
            
            # Skip unwanted modes early
            if requested_filter_mode != "all" and filter_mode != requested_filter_mode:
                continue

            rt_list.append(entry["scanList"]["scan"][0]["scan time"])

            # Convert arrays to NumPy for vectorized masking
            mz_array = np.asarray(entry["m/z array"])
            int_array = np.asarray(entry["intensity array"])

            mask = (mz_array >= mass - mass_deviation) & (mz_array <= mass + mass_deviation)
            curr_sum_int = np.sum(int_array[mask]) if np.any(mask) else 0

            intensity_list.append(curr_sum_int)
            original_index_list.append(idx)
        return [rt_list, intensity_list, original_index_list]

    def extract_key_value_pairs(self, d, parent_key=''):
        items = []
        if isinstance(d, dict):
            for k, v in d.items():
                new_key = f"{parent_key}.{k}" if parent_key else k
                if isinstance(v, dict):
                    items.extend(self.extract_key_value_pairs(v, new_key))
                elif isinstance(v, list):
                    for i, item in enumerate(v):
                        indexed_key = f"{new_key}[{i}]"
                        if isinstance(item, dict):
                            items.extend(self.extract_key_value_pairs(item, indexed_key))
                        else:
                            items.append((indexed_key, item))
                else:
                    items.append((new_key, v))
        else:
            items.append((parent_key, d))
        return items

    def get_nested_value(self, data, path, delete=False):
        dict_value = None
        try:
            keys = re.split(r'\.(?![^\[]*\])', path)
            for i, key in enumerate(keys):
                list_match = re.match(r'([^\[]+)\[(\d+)\]', key)
                if list_match:
                    dict_key = list_match.group(1)
                    index = int(list_match.group(2))
                    if i == len(keys) - 1:
                        dict_value = data[dict_key][index]
                        if delete:
                            del data[dict_key][index]
                    else:
                        data = data[dict_key][index]
                else:
                    if i == len(keys) - 1:
                        dict_value = data[key]
                        if delete:
                            del data[key]
                    else:
                        data = data[key]
        except Exception as e:
            print(f"Error accessing path '{path}': {e}")
            dict_value = None
        return dict_value

    def rename_rawdata_keys(self):
        all_kv_pairs = self.extract_key_value_pairs(self.rawdata[-1])
        example_rawdata = self.rawdata[-2]
        rd_keys = list(example_rawdata.keys())
        if not "m/z array" in rd_keys:
            print("No m/z array found in the rawdata. Trying to find the appropriate key...")
            for key in rd_keys:
                if ( ("m/z" in key.lower()) or ("mz" in key.lower()) or ("mass" in key.lower()) ) and \
                    ( ("array" in key.lower()) or ("list" in key.lower()) or ("values" in key.lower()) ) and \
                        ( isinstance(self.rawdata[-2][key], tuple) or isinstance(self.rawdata[-2][key], list) or isinstance(self.rawdata[-2][key], np.ndarray) ):
                    print("Found m/z array key: " + str(key))
                    for index in range(len(self.rawdata)):
                        self.rawdata[index]["m/z array"] = list(self.rawdata[index][key])
                        del self.rawdata[index][key]
        all_kv_pairs = self.extract_key_value_pairs(self.rawdata[-1])
        if not "intensity array" in rd_keys:
            print("No intensity array found in the rawdata. Trying to find the appropriate key...")
            for key in rd_keys:
                if ( ("intensity" in key.lower()) or ("intensities" in key.lower()) or ("int" in key.lower()) ) and \
                      ( ("array" in key.lower()) or ("list" in key.lower()) or ("values" in key.lower()) ) and \
                        ( isinstance(self.rawdata[-2][key], tuple) or isinstance(self.rawdata[-2][key], list) or isinstance(self.rawdata[-2][key], np.ndarray) ):
                    print("Found intensity array key: " + str(key))
                    for index in range(len(self.rawdata)):
                        self.rawdata[index]["intensity array"] = list(self.rawdata[index][key])
                        del self.rawdata[index][key]
        all_kv_pairs = self.extract_key_value_pairs(self.rawdata[-1])
        if not "scanList" in rd_keys:
            print("No scanList found in the rawdata. Trying to find the appropriate key...")
            for key in rd_keys:
                if ( ("scan list" in key.lower()) or ("scanlist" in key.lower()) or ("list of scans" in key.lower()) or \
                          ("scan dict" in key.lower()) or ("scandict" in key.lower()) or ("description" in key.lower()) ) and \
                                ( isinstance(self.rawdata[-2][key], dict) ):
                    print("Found scanList key: " + str(key))
                    for index in range(len(self.rawdata)):
                        self.rawdata[index]["scanList"] = self.rawdata[index][key]
                        del self.rawdata[index][key]
        all_kv_pairs = self.extract_key_value_pairs(self.rawdata[-1])
        if not "total ion current" in rd_keys:
            print("No total ion current found in the rawdata. Trying to find the appropriate key...")
            for key, value in all_kv_pairs:
                if ( ("total ion current" in key.lower()) or ("tic" in key.lower()) or ("total ion current signal" in key.lower()) or ("tic int" in key.lower()) or ("ticint" in key.lower()) ):
                      example_tic = value
                      if ( isinstance(example_tic, float) or isinstance(example_tic, int) ):
                        print("Found total ion current key: " + str(key))
                        for index in range(len(self.rawdata)):
                            curr_tic_val = self.get_nested_value(self.rawdata[index], key, delete=True)
                            self.rawdata[index]["total ion current"] = curr_tic_val

    def rename_scanlist_keys(self):
        scanlist_kv_pairs = self.extract_key_value_pairs(self.rawdata[-1]["scanList"])
        scanlist_keys = list(self.rawdata[-1]["scanList"].keys())
        print(scanlist_keys)
        if not "scan" in scanlist_keys:
            print("No scan found in the scanList. Trying to find the appropriate key...")
            for key in scanlist_keys:
                if ("scan" in key.lower() or "scans" in key.lower()) and \
                    ( isinstance(self.rawdata[-1]["scanList"][key], list) or isinstance(self.rawdata[-1]["scanList"][key], tuple) ):
                    print("Found scanList -> scan key: " + str(key))
                    for index in range(len(self.rawdata)):
                        self.rawdata[index]["scanList"]["scan"] = list(self.rawdata[index]["scanList"][key])
                        del self.rawdata[index]["scanList"][key]
        
        scanlist_scan_keys = list(self.rawdata[-1]["scanList"]["scan"][0].keys())
        if not "scan time" in scanlist_scan_keys:
            print("No scan time found in the scanList -> scan. Trying to find the appropriate key...")
            for key in scanlist_scan_keys:
                if ("rt" in key.lower() or "scanstart" in key.lower().replace(" ", "") or "scantime" in key.lower().replace(" ", "")) and \
                    ( isinstance(self.rawdata[-1]["scanList"]["scan"][0][key], float) or isinstance(self.rawdata[-1]["scanList"]["scan"][0][key], int) ):
                    print("Found scanList -> scan -> scan time key: " + str(key))
                    for index in range(len(self.rawdata)):
                        self.rawdata[index]["scanList"]["scan"][0]["scan time"] = self.rawdata[index]["scanList"]["scan"][0][key]
                        del self.rawdata[index]["scanList"]["scan"][0][key]
        if not "filter string" in scanlist_scan_keys:
            print("No filter string found in the scanList -> scan. Trying to find the appropriate key...")
            for key in scanlist_scan_keys:
                if ("filter" in key.lower() or "filter string" in key.lower() or "filt" in key.lower()) and \
                    ( isinstance(self.rawdata[-1]["scanList"]["scan"][0][key], str) ):
                    print("Found scanList -> scan -> filter string key: " + str(key))
                    for index in range(len(self.rawdata)):
                        self.rawdata[index]["scanList"]["scan"][0]["filter string"] = self.rawdata[index]["scanList"]["scan"][0][key]
                        del self.rawdata[index]["scanList"]["scan"][0][key]

    def construct_filter_string(self, rawdata_index):
        # normal filter string looks like: FTMS - p ESI Full ms [50.0000-750.0000]
        upper_mass = max(self.rawdata[rawdata_index]["m/z array"])
        lower_mass = min(self.rawdata[rawdata_index]["m/z array"])
        upper_mass = round(upper_mass, 4)
        lower_mass = round(lower_mass, 4)
        if upper_mass < lower_mass:
            print("Upper mass is lower than lower mass. Setting to default values.")
            upper_mass = 50000
            lower_mass = 0
        if str(upper_mass).lower() == "nan" or str(upper_mass).lower() == "inf":
            print("Upper mass is NaN or Inf. Setting upper mass to 50000 and lower mass to 0.")
            upper_mass = 50000
        if str(lower_mass).lower() == "nan" or str(lower_mass).lower() == "-inf":
            print("Lower mass is NaN or -Inf. Setting lower mass to 0.")
            lower_mass = 0
        msn = 0
        polarity = "NA"
        all_data_kv_pairs = self.extract_key_value_pairs(self.rawdata[rawdata_index])
        for key, val in all_data_kv_pairs:
            if ("mslevel" in key.lower().replace(" ", "") or "msn" in key.lower()) and isinstance(value, int):
                msn = val
                break
            elif ("mslevel" in key.lower().replace(" ", "") or "msn" in key.lower()) and isinstance(value, str):
                try:
                    msn = int(val)
                    break
                except:
                    try:
                        msn = float(val)
                        break
                    except:
                        print("Could not get the ms level by key value pair. Trying only keys or values now...")
            try:
                value = str(val).lower().replace(" ", "")
                key = str(key).lower().replace(" ", "")
                if value == "ms1" or value == "fullscan":
                    msn = 1
                    break
                elif value == "ms2" or value == "msn" or value == "aif" or value.replace("/", "") == "msms":
                    msn = 2
                    break
                elif value == "ms3":
                    msn = 3
                    break
                elif key == "fullscan" or key == "ms1":
                    msn = 1
                    break
                elif key == "ms2" or key == "msn" or key == "aif" or key.replace("/", "") == "msms":
                    msn = 2
                    break
                else:
                    print("No ms level found in key value pairs. Setting ms level to 1 as default.")
                    msn = 1
            except Exception as e:
                print("Error while trying to convert ms level: " + str(e))
                print(traceback.format_exc())
                msn = 1
        filter_string = "FTMS - p ESI "
        if msn == 1:
            filter_string += "Full ms [" + str(round(lower_mass, 4)) + "-" + str(round(upper_mass, 4)) + "]"
        elif msn == 2:
            filter_string += "hcd [" + str(round(lower_mass, 4)) + "-" + str(round(upper_mass, 4)) + "]"
        print(filter_string)
        return filter_string
                                    
    def get_background_spectra(self, background_range=(3, 10), bckg_m_dev=0.0001, blank_multiplicator=3):
        if "AIF" in self.available_modes and "Full scan" in self.available_modes:
            print("AIF and MS1 spectra available. Proceeding with background calculation...")
        else:
            print("AIF and MS1 spectra not available. Cannot calculate background spectrum. Returning...")
            self.save_ms_file_log_entry("ERROR:\t" + "AIF and MS1 spectra not available. Cannot calculate background spectrum. Returning...")
            self.aif_background_spectrum = None
            self.ms1_background_spectrum = None
            return False
        averaged_aif_spectrum = {}
        number_aif_spectra = 0
        averaged_ms1_spectrum = {}
        number_ms1_spectra = 0
        for index in range(background_range[0], background_range[1]):
            if self.all_modes[index] == "AIF":
                number_aif_spectra += 1
                aif_masses = self.rawdata[index]["m/z array"]
                aif_intensity = self.rawdata[index]["intensity array"]
                curr_mass_int_dict = dict(zip(aif_masses, aif_intensity))
                curr_mass_int_dict = MS_functions.summarize_mass_intensity_dict(curr_mass_int_dict, deviation=11, debug_output=True)
                
                for m, i in curr_mass_int_dict.items():
                    for existing_mass in list(averaged_aif_spectrum.keys()):
                        if abs(existing_mass - m) <= bckg_m_dev:
                            averaged_aif_spectrum[existing_mass] += i
                            break
                    averaged_aif_spectrum[m] = i
            elif self.all_modes[index] == "Full scan":
                number_ms1_spectra += 1
                ms1_masses = self.rawdata[index]["m/z array"]
                ms1_intensity = self.rawdata[index]["intensity array"]
                curr_mass_int_dict = dict(zip(ms1_masses, ms1_intensity))
                curr_mass_int_dict = MS_functions.summarize_mass_intensity_dict(curr_mass_int_dict, deviation=11, debug_output=True)
                for m, i in curr_mass_int_dict.items():
                    for existing_mass in list(averaged_ms1_spectrum.keys()):
                        if abs(existing_mass - m) <= bckg_m_dev:
                            averaged_ms1_spectrum[existing_mass] += i
                            break
                    averaged_ms1_spectrum[m] = i
        averaged_aif_spectrum = MS_functions.summarize_mass_intensity_dict(averaged_aif_spectrum, deviation=11, debug_output=True)
        averaged_ms1_spectrum = MS_functions.summarize_mass_intensity_dict(averaged_ms1_spectrum, deviation=11, debug_output=True)
        averaged_aif_spectrum = {k: ((v / number_aif_spectra)*blank_multiplicator) for k, v in averaged_aif_spectrum.items()}
        averaged_ms1_spectrum = {k: ((v / number_ms1_spectra)*blank_multiplicator) for k, v in averaged_ms1_spectrum.items()}
        self.aif_background_spectrum = averaged_aif_spectrum
        self.ms1_background_spectrum = averaged_ms1_spectrum
        return averaged_aif_spectrum, averaged_ms1_spectrum

    def do_background_subtraction(self, background_range=(3,10), bckg_m_dev=0.0001):
        if self.aif_background_spectrum is None or self.ms1_background_spectrum is None:
            self.get_background_spectra(background_range=background_range, bckg_m_dev=bckg_m_dev)
        
        for index in range(len(self.rawdata)):
            print("Doing background subtraction for index: " + str(index) + " / " + str(len(self.rawdata)))
            curr_masses = self.rawdata[index]["m/z array"]
            curr_intensity = self.rawdata[index]["intensity array"]
            curr_dict = dict(zip(curr_masses, curr_intensity))
            curr_dict = MS_functions.summarize_mass_intensity_dict(curr_dict, deviation=11, debug_output=False)
            if self.all_modes[index] == "AIF":
                for curr_mass, curr_int in curr_dict.items():
                    for mass_bckg, intensity_bckg in self.aif_background_spectrum.items():
                        if abs(mass_bckg - curr_mass) <= bckg_m_dev:
                            curr_dict[curr_mass] -= intensity_bckg
                curr_masses = list(curr_dict.keys())
                curr_intensity = list(curr_dict.values())
                self.rawdata[index]["m/z array"] = curr_masses
                self.rawdata[index]["intensity array"] = curr_intensity
            elif self.all_modes[index] == "Full scan":
                for curr_mass, curr_int in curr_dict.items():
                    for mass_bckg, intensity_bckg in self.ms1_background_spectrum.items():
                        if abs(mass_bckg - curr_mass) <= bckg_m_dev:
                            curr_dict[curr_mass] -= intensity_bckg
                curr_masses = list(curr_dict.keys())
                curr_intensity = list(curr_dict.values())
                self.rawdata[index]["m/z array"] = curr_masses
                self.rawdata[index]["intensity array"] = curr_intensity

    def save_ms_file_log_entry(self, log_entry):
        #get the dirname of the logfile_filepath
        directory_logfile = os.path.dirname(self.kwargs["logfile_filepath"])
        #create the directory if it does not exist
        os.makedirs(directory_logfile, exist_ok=True)
        log_f = open(self.kwargs["logfile_filepath"], "a")
        log_f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\t" + log_entry + "\n")
        log_f.close()
        return True