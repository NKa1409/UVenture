import UVenture



ms_file = UVenture.MS_File("webserver_save/mzml_files/Shark_10_07_24_Spritzenfilter_1.mzML")
print(ms_file)

myspec = UVenture.Spec(ms_file, 426)
print(myspec)

my_prediction = UVenture.OneAnalysis(ms_file, 139.0069, 270.35)
