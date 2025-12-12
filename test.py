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


    #mass, ints = mzml_functions.calculate_averaged_spectrum("22_summed_NatKaffee_file.mzML", [1, 3, 5, 7])
    #import matplotlib.pyplot as plt
    #fig, ax = plt.subplots(nrows=1, ncols=1)
    #ax.bar(mass, ints, width=0.5)
    #fig.savefig("test.png")
    #sys.exit()    

    #mzml_functions.remove_background("22_summed_NatKaffee_file.mzML", remove_completely=True, background_signals=None, save_mzml_filepath="dynamic_bgremoved.mzML", max_signals=20)
    #mzml_functions.remove_background("22_summed_NatKaffee_file.mzML", background_signals=[100.0000], save_mzml_filepath="static_bgremoved.mzML", remove_completely=True, death_time=0.1)
    #sys.exit()


    ms_files_cutout = [os.path.join("1000ppb", str("1000ppb_20.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_21.mzML")),
                       os.path.join("1000ppb", str("1000ppb_22.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_23.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_24.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_25.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_26.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_27.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_28.mzML"))]
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_20_28_summed_Nitrophenole1000ppb.mzML")

    ms_files_cutout = [os.path.join("1000ppb", str("1000ppb_50.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_51.mzML")),
                       os.path.join("1000ppb", str("1000ppb_52.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_53.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_54.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_55.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_56.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_57.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_58.mzML"))]
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_50_58_summed_Nitrophenole1000ppb.mzML")

    ms_files_cutout = [os.path.join("1000ppb", str("1000ppb_90.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_91.mzML")),
                       os.path.join("1000ppb", str("1000ppb_92.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_93.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_94.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_95.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_96.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_97.mzML")), 
                       os.path.join("1000ppb", str("1000ppb_98.mzML"))]
    ms_files_cutout.sort()
    print(ms_files_cutout)
    summed_ms_obj = mzml_functions.sum_multiple_mzmlfiles(ms_files_cutout)
    summed_ms_obj.save_to_mzml_file(str(len(ms_files_cutout)) + "_90_98_summed_Nitrophenole1000ppb.mzML")

    sys.exit()





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
