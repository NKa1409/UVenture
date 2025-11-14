import copy
import os
import sys
import numpy as np
import UVenture.class_MS_file as class_MS_file
import UVenture.mzml_functions as mzml_functions
import UVenture.plotting as plotting
import UVenture.class_Spec as class_Spec

if __name__ == "__main__":

    
    '''
    single_ms_file = class_MS_file.MS_File(mzml_filenames_2ppb[0])
    single_xic = single_ms_file.get_xic(138.0196, 0.001)
    plotting.create_xic(single_xic[0], single_xic[1], "Single, 138_0196 - Nitrophenole (most intense)", os.path.join("res", "Nitrophenole_single.png"))
    single_xic = single_ms_file.get_xic(139.0229, 0.0005)
    plotting.create_xic(single_xic[0], single_xic[1], "Single, 139_0229 - 13C-Nitrophenole", os.path.join("res", "Nitrophenole_13c_single.png"))
    single_xic = single_ms_file.get_xic(139.0238, 0.0004)
    plotting.create_xic(single_xic[0], single_xic[1], "Single, 139_0238 - 17O-Nitrophenole", os.path.join("res", "Nitrophenole_17o_single.png"))
    
    
    for i in range(1, len(mzml_filenames_2ppb), 1):
        summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(mzml_filenames_2ppb[:i+1])
        multiple_xic = summed_ms_obj.get_xic(138.0196, 0.001)
        plotting.create_xic(multiple_xic[0], multiple_xic[1], str(i) + "multiple_138.0196 - Nitrophenole (most intense)", os.path.join("res", str(i) + "multiple_Nitrophenole.png"))
        multiple_xic = summed_ms_obj.get_xic(139.0229, 0.0005)
        plotting.create_xic(multiple_xic[0], multiple_xic[1], str(i) + "multiple_139.0229 - 13C-Nitrophenole", os.path.join("res", str(i) + "multiple_Nitrophenole_13c.png"))
        multiple_xic = summed_ms_obj.get_xic(139.0238, 0.0004)
        plotting.create_xic(multiple_xic[0], multiple_xic[1], str(i) + "multiple_139.0238 - 17O-Nitrophenole", os.path.join("res", str(i) + "multiple_Nitrophenole_17o.png"))

    '''





    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("coffeinsynthetisch") and "Verd0_25" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_SynKoffein_Verd0_25_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("coffeinsynthetisch") and "Verd0_5" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_SynKoffein_Verd0_5_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("coffeinsynthetisch") and "Verd1" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_SynKoffein_Verd1_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("coffeinsynthetisch") and "Verd2" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_SynKoffein_Verd2_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("coffeinsynthetisch") and "Verd3" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_SynKoffein_Verd3_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("coffeinsynthetisch") and "Verd4" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_SynKoffein_Verd4_file.mzML")



    ###########################


    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("kaffeenatuerlich") and "Verd0_25" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_natuerlichKaffee_Verd0_25_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("kaffeenatuerlich") and "Verd0_5" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_natuerlichKaffee_Verd0_5_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("kaffeenatuerlich") and "Verd1" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_natuerlichKaffee_Verd1_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("kaffeenatuerlich") and "Verd2" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_natuerlichKaffee_Verd2_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("kaffeenatuerlich") and "Verd3" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_natuerlichKaffee_Verd3_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("kaffeenatuerlich") and "Verd4" in mzmlfile and not "high" in mzmlfile:
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_natuerlichKaffee_Verd4_file.mzML")


    #################################
    # All concentrations together
    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("kaffeenatuerlich"):
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_natuerlichKaffee_AllConcentrationsTogether_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("KoffeinStandards"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("coffeinsynthetisch"):
            ms_files_cutout.append( os.path.join("KoffeinStandards", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_SynKoffein_AllConcentrationsTogether_file.mzML")


    ####################################

    ms_files_cutout = []
    for mzmlfile in os.listdir("Koffein"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("nat_kaffee"):
            ms_files_cutout.append( os.path.join("Koffein", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_NatKaffee_file.mzML")

    ms_files_cutout = []
    for mzmlfile in os.listdir("Koffein"):
        if mzmlfile.lower().endswith(".mzml") and mzmlfile.lower().startswith("syn_coffein"):
            ms_files_cutout.append( os.path.join("Koffein", str(mzmlfile)) )
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_SynKoffein_file.mzML")



    sys.exit()

    ms_files_cutout = ms_filename_list[:30]
    print("Starting with 30 files...")
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    print("Saving MS FIle Object")
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_file.mzML")
    print("Object saved...")

    ms_files_cutout = ms_filename_list[:100]
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_file.mzML")

    ms_files_cutout = ms_filename_list[:10]
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_file.mzML")

    ms_files_cutout = ms_filename_list[:80]
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_file.mzML")

    ms_files_cutout = ms_filename_list[:50]
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_file.mzML")

    ms_files_cutout = ms_filename_list[:90]
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_file.mzML")

    ms_files_cutout = ms_filename_list[:70]
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_file.mzML")

    ms_files_cutout = ms_filename_list[:60]
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_file.mzML")

    ms_files_cutout = ms_filename_list[:40]
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_file.mzML")
        
    ms_files_cutout = ms_filename_list[:20]
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_summed_file.mzML")
