import ast
import base64
import datetime
import hashlib
import io
import os
import shutil
import signal
import sys
import threading
import time
import traceback
import zipfile
import flask
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
import werkzeug
import UVenture
import MS_functions
import multiprocessing
import UVenture_peakdetection


def caller_func(ms_filepath, mz, rt, settings_dict, parentfolder_msfile):
    try:
        print("Starting analysis for m/z: " + str(mz) + ", rt: " + str(rt) + ", ms_filepath: " + ms_filepath)
        ms_file = UVenture.MS_File(ms_filepath, parentfolder_msfile=parentfolder_msfile, **settings_dict)
        print("MS file loaded")
        myanalysis = UVenture.OneAnalysis(ms_file, mz, rt, **settings_dict)
        return
    except Exception as e:
        print(e)
        print(traceback.format_exc())
        return

def start_pd_process(ms_filepath, parentfolder, peaklist_filename, mzrt_filename, settings_dict):
    ms_file = UVenture.MS_File(ms_filepath, parentfolder_msfile=parentfolder, **settings_dict)
    UVenture_peakdetection.get_all_possible_peaks_2(ms_file, settings_dict, mass_range=1, 
                                        threshold_area=150000, threshold_intensity=50000, 
                                        rt_bins=400, mass_deviation_isotopo=settings_dict["mass_deviation"]/5, height_deviation_isotopo=0.5,
                                        min_peak_width=4, max_peak_width=40,
                                        peaklist_filename=peaklist_filename, mzrt_filename=mzrt_filename)
    return


def resource_path(relative_path):
    """ Get the absolute path to a resource, works for dev and PyInstaller .exe """
    try:
        # PyInstaller creates a temp folder and stores the path in _MEIPASS
        base_path = sys._MEIPASS
        return os.path.join(base_path, relative_path)
    except AttributeError:
        base_path = os.path.abspath(".")
        return os.path.join("./", relative_path)


class Webpage:
    def __init__(self) -> None:
        self.version = "1.0.1"
        self.python_version = "Python " + str(sys.version_info.major) + "." + str(sys.version_info.minor) + "." + str(sys.version_info.micro)
        self.path_of_this_file = os.path.abspath(__file__)
        self.file_modification_time = os.path.getmtime(self.path_of_this_file)
        with open(self.path_of_this_file, "rb") as f:
            bytes = f.read()
            self.readable_hash = hashlib.sha256(bytes).hexdigest()
        self.start_time = datetime.datetime.now()
        self.running_processes = []
        self.process_times = {} # Will be populated with the starttime and the endtime of the process {proc_uuid: (starttime, endtime, file, )}
        self.best_time_approx_per_analysis = 0

        self.num_cores_to_use = os.cpu_count() - 6
        if self.num_cores_to_use < 1:
            self.num_cores_to_use = 1

        self.parentfolder = resource_path("webserver_save/")
        print("Parent folder: " + self.parentfolder)
        os.makedirs(self.parentfolder, exist_ok=True)
        self.mzml_folder = resource_path(self.parentfolder + "mzml_files/")
        print("mzML folder: " + self.mzml_folder)
        os.makedirs(self.mzml_folder, exist_ok=True)
        self.results_folder = resource_path(self.parentfolder + "results/")
        print("Results folder: " + self.results_folder)
        os.makedirs(self.results_folder, exist_ok=True)
        self.app = flask.Flask(__name__)
        self.server = None
        self.available_files = os.listdir(self.mzml_folder)
        self.available_files = [f for f in self.available_files if f.endswith(".mzML")]
        self.curr_ms_file = None
        self.curr_xic_encoded_plot = {}
        self.curr_spec_encoded_plot = {}
        self.curr_mass_deviation = 0
        os.makedirs("static", exist_ok=True)
        self.settings_file_filepath = resource_path("static/settings.txt")
        self.settings_default_file_filepath = resource_path("static/settings_default.txt")
        self.help_page_contents_filepath = resource_path("static/help_page_contents.txt")
        self.taskstorage_filepath = resource_path("static/taskstorage.txt")
        #Clear the task storage file
        with open(self.taskstorage_filepath, "w") as f:
            f.write("")

        @self.app.route("/", methods=["GET", "POST"])
        def index():
            self.start_background_task_checking()
            if len(self.available_files) == 0:
                return flask.redirect("/upload_mzml_file")
            else:
                return flask.redirect("/queue_new_analysis")
        
        @self.app.route("/queue_new_analysis", methods=["GET", "POST"])
        def queue_new_analysis():
            self.available_files = os.listdir(self.mzml_folder)
            self.available_files = [f for f in self.available_files if f.endswith(".mzML")]
            self.start_background_task_checking()
            settings_dict = self.get_settings_dict()
            if flask.request.method == "POST":
                form_data = flask.request.form.to_dict()
                if not "peak_analysis_cb" in form_data:
                    form_data["peak_analysis_cb"] = "false"
                if not "mass_analysis_cb" in form_data:
                    form_data["mass_analysis_cb"] = "false"
                if not "peak_analysis_list_cb" in form_data:
                    form_data["peak_analysis_list_cb"] = "false"
                if form_data["peak_analysis_cb"] == "false" and form_data["mass_analysis_cb"] == "false" and form_data["peak_analysis_list_cb"] == "false":
                    return flask.render_template_string("No analysis mode chosen. \n Please tick the box of the analysis that you want to perform.")
                print(form_data)
                keys_to_check = ["mz_peak_analysis", "mz_mass_analysis", "retention_time"]
                for k in keys_to_check:
                    try:
                        form_data[k] = float(form_data[k])
                    except:
                        pass
                try:
                    form_data["spec_index"] = int(form_data["spec_index"])
                except:
                    pass
                
                print("MS file loaded")
                if form_data["peak_analysis_cb"] == "true":
                    if form_data["mz_peak_analysis"] == "":
                        return flask.render_template_string("No m/z given! Cannot analyze peak without a mass. \n Please enter a peak to analyse.")
                    if form_data["retention_time"] == "":
                        return flask.render_template_string("No retention time given! Cannot analyze peak without retention time. \nPlease enter a retention time.")
                    try:
                        with open(self.taskstorage_filepath, "a") as f:
                            f.write(form_data["fileselection"] + "\t" + str(form_data["mz_peak_analysis"]) + "\t" + str(form_data["retention_time"]) + "\t" + str(settings_dict) + "\n")
                        print("Task added to task storage file")
                    except Exception as e:
                        print(e)
                        print(traceback.format_exc())
                        return flask.render_template_string("An error occured during the analysis: " + str(e) + "\n \n \n" + str(traceback.format_exc()))
                    print("Peak analysis started")

                if form_data["mass_analysis_cb"] == "true":
                    if form_data["mz_mass_analysis"] == "":
                        return flask.render_template_string("No m/z given! Cannot analyze peak without a mass. \n Please enter a peak to analyse")
                    if form_data["spec_index"] == "":
                        return flask.render_template_string("No retention time given! Cannot analyze peak without retention time. \nPlease enter a retention time")
                    def run_prediction_analysis():
                        try:
                            ms_filepath = self.mzml_folder + form_data["fileselection"]
                            ms_file = UVenture.MS_File(ms_filepath, parentfolder_msfile=self.results_folder + str(".".join(form_data["fileselection"].split(".")[:-1])) + "/")
                            settings_dict["spec_requested_filter_mode"] = settings_dict["spec_requested_filter_mode"]
                            myspec = UVenture.Spec(ms_file, form_data["spec_index"], **settings_dict)
                            print("Spec created")
                            print("Starting mass prediction")
                            myanalysis = UVenture.Prediction(ms_file, form_data["mz_mass_analysis"], myspec, **settings_dict)
                        except Exception as e:
                            print(e)
                            print(traceback.format_exc())
                    thread_name = "UVenture_Prediction_starttime_" + str(datetime.datetime.now().strftime("%Y%m%d:%H%M%S")) + "_file_" + str(form_data["fileselection"]) + "_mass_" + str(form_data["mz_mass_analysis"]) + "_specindex_" + str(form_data["spec_index"])
                    thread = threading.Thread(target=run_prediction_analysis, name=thread_name)
                    thread.start()
                    print("Mass analysis started")
                
                if form_data["peak_analysis_list_cb"] == "true":
                    if "peak_list_file" in flask.request.files:
                        peak_list_file = flask.request.files["peak_list_file"]
                        peak_list_file = io.StringIO(peak_list_file.read().decode("utf-8"))
                        print(peak_list_file)
                    else:
                        return flask.render_template_string("No peak list given! Cannot analyze peaks without a list. \n Please enter a peak list to analyse.")

                    lines = peak_list_file.readlines()
                    with open(self.taskstorage_filepath, "a") as f:
                        for l in lines:
                            l = l.strip()
                            f.write(form_data["fileselection"])
                            f.write("\t")
                            f.write(l)
                            f.write("\t")
                            f.write(str(settings_dict))
                            f.write("\n")
                
                print("Analysis started")
                return flask.redirect("/show_currently_running")
            return flask.render_template("queue_new_analysis.html", files=self.available_files)


        @self.app.route("/main_settings", methods=["GET", "POST"])
        def main_settings():
            main_settings = ["mass_deviation", 
                             "charge_of_measured_mass", 
                             "msfile_raw_file_retention_time_unit", 
                             "pred_atoms_to_keep_in_prediction", 
                             "oa_fragments_do_peak_computation_if_area_higher_than", 
                             "oa_fragments_absolute_max_number_of_fragment_masses", 
                             "pred_minimum_assumed_noise"]
            settings_dict = self.get_settings_dict()
            main_settings_dict = {k:v for k, v in settings_dict.items() if k in main_settings}

            if flask.request.method == "POST":
                form_data = flask.request.form.to_dict()
                print(form_data)
                if flask.request.form.get("button") == "save_button":
                    new_main_settings_dict = {}
                    for key in main_settings:
                        new_main_settings_dict[key] = form_data.get(key, "")
                    print(new_main_settings_dict)
                    write_settings_dict = {**settings_dict, **new_main_settings_dict}
                    with open(self.settings_file_filepath, "w") as f:
                        for key, value in write_settings_dict.items():
                            f.write(key + "=" + str(value) + "\n")
                    print("New settings saved")
                    return flask.redirect("/main_settings")
            return flask.render_template("main_settings.html", settings_dict=main_settings_dict)

        @self.app.route("/change_settings", methods=["GET", "POST"])
        def change_settings():
            settings_dict = self.get_settings_dict()
            
            if flask.request.method == "POST":
                if flask.request.form.get("button") == "save_button":
                    form_data = flask.request.form.to_dict()
                    print(form_data)
                    del form_data["button"]
                    with open(self.settings_file_filepath, "w") as f:
                        for key in form_data:
                            f.write(key + "=" + str(form_data[key]) + "\n")
                    print("New settings saved")
                    print(form_data)
                if flask.request.form.get("button") == "restore_default_values":
                    os.remove(self.settings_file_filepath)
                    shutil.copyfile(self.settings_default_file_filepath, self.settings_file_filepath)
                return flask.redirect("/change_settings")
                
            return flask.render_template("change_settings.html", settings_dict=settings_dict)

        @self.app.route("/resultsdownload", methods=["GET", "POST"])
        def resultsdownload():
            self.start_background_task_checking()
            if flask.request.method == "POST":
                if "delete_once_downloaded" in flask.request.form.to_dict():
                    @flask.after_this_request
                    def remove_folder(response):
                        folderpath = self.results_folder + ".".join(flask.request.form["fileselection"].split(".")[:-1]) + "/"
                        try:
                            shutil.rmtree(folderpath)  # Delete the file after sending
                        except Exception as e:
                            print(f"Error deleting file: {e}")
                        return response
                if "download_summary_only" in flask.request.form.to_dict():
                    filepath = self.results_folder + ".".join(flask.request.form["fileselection"].split(".")[:-1]) + "/" + "SUMMARY.txt"
                    if not os.path.exists(filepath):
                        return flask.render_template_string("Summary file does not exist!")
                    else:
                        return flask.send_file(filepath, as_attachment=True, download_name=flask.request.form["fileselection"] + "SUMMARY.txt")
                folderpath = self.results_folder + ".".join(flask.request.form["fileselection"].split(".")[:-1]) + "/"
                print(folderpath)
                # Create a temporary zip file in memory
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                    for root, dirs, files in os.walk(folderpath):
                        for file in files:
                            zip_file.write(os.path.join(root, file), 
                                           os.path.relpath(os.path.join(root, file), 
                                           os.path.join(folderpath, '..')))
                zip_buffer.seek(0)
                return flask.send_file(zip_buffer, as_attachment=True, download_name="results.zip")
            return flask.render_template("resultsdownload.html", files=self.available_files, basefolder=self.results_folder)

        @self.app.route("/show_currently_running", methods=["GET", "POST"])
        def show_currently_running():
            self.start_background_task_checking()

            # Get the list of scheduled tasks from the task storage file
            with open(self.taskstorage_filepath, "r") as f:
                lines = f.readlines()
            lines = [line for line in lines if line.strip() != "" and line[0] != "#" and "\t" in line]
            pending_analyses = len(lines)
            try:
                unique_analysis_files = [f.split("\t")[0] for f in lines]
                unique_analysis_files = list(set(unique_analysis_files))
                unique_analysis_files = [f.split("/")[-1] for f in unique_analysis_files]
            except Exception as e:
                print(e)
                print(traceback.format_exc())
                unique_analysis_files = []

            all_threads = threading.enumerate()
            running_threads = [t.name for t in all_threads if t.is_alive()]
            running_processes = [p.name for p in multiprocessing.active_children()]
            running_analyses = running_threads + running_processes
            running_analyses = [a for a in running_analyses if "UVenture_" in a]
            #Extract the start times of the individual threads
            try:
                self.update_best_time_approx_per_analysis()
                best_esimation_for_process_time = self.best_time_approx_per_analysis
            except Exception as e:
                best_esimation_for_process_time = 0
                print(e)
                print(traceback.format_exc())
                

            tasks_information = []
            for task in running_analyses:
                tasks_information.append({"name": task, "is_alive": True})
            tasks_information.sort(key=lambda x: x["name"])
            tasks_information.append({"name": "█████████████████████████████████████████████████████████████████████", "is_alive": True})
            for task in all_threads:
                if not task.name.startswith("UVenture_"):
                    tasks_information.append({"name": task.name, "is_alive": task.is_alive()})
            if not len(tasks_information) <= 3:
                time_till_finished = ((pending_analyses * best_esimation_for_process_time) / (len(tasks_information)-1)) + best_esimation_for_process_time
                time_till_finished = int(time_till_finished)
                time_till_finished_hours = time_till_finished / 3600
                time_till_finished = datetime.timedelta(seconds=time_till_finished)
                finish_time = datetime.datetime.now() + time_till_finished
                finish_time = finish_time.strftime("%H:%M:%S  %Y-%m-%d")
            else:
                time_till_finished = 0
                time_till_finished_hours = 0
                finish_time = "N/A"
            return flask.render_template("show_currently_running.html", tasks_information=tasks_information, 
                                         pending_analyses=pending_analyses, 
                                         unique_analysis_files=unique_analysis_files, 
                                         best_esimation_for_process_time=best_esimation_for_process_time, 
                                         time_till_finished_hours=time_till_finished_hours, 
                                         finish_time=finish_time)

        @self.app.route("/mzml_file_viewer", methods=["GET", "POST"])
        def mzml_file_viewer():
            return flask.render_template("mzml_file_viewer.html", files=self.available_files)

        @self.app.route("/update_mzml_plot_viewer_plots", methods=["GET", "POST"])
        def update_mzml_plot_viewer_plots():
            if flask.request.method == "POST":
                print("POST request received")
                # Parse the JSON request
                data = flask.request.get_json()
                print(data)
                file_select = ""
                try:
                    xic_mass = float(data['xic_mass'])
                except:
                    xic_mass = 150
                try:
                    mass_deviation = float(data['mass_deviation'])
                except:
                    mass_deviation = 0.001
                try:
                    spec_index = int(data['spec_index'])
                except:
                    spec_index = 1
                file_select = data['fileSelect']
                if file_select == "":
                    print("No file selected!")
                    return "No file selected!"
                print(xic_mass, spec_index, file_select)
                calc_new_file = True
                try:
                    print(self.curr_ms_file.filename, file_select)
                    if self.curr_ms_file.filename == self.mzml_folder + file_select:
                        print("Same file selected")
                        calc_new_file = False
                    else:
                        print("Different file selected")
                        calc_new_file = True
                except:
                    print("First file selected")
                    calc_new_file = True
                if calc_new_file:
                    ms_filepath = self.mzml_folder + file_select
                    self.curr_ms_file = UVenture.MS_File(ms_filepath, parentfolder_msfile=self.results_folder + str(".".join(file_select.split(".")[:-1])) + "/")
                print("MS file loaded")
                if calc_new_file == False and (xic_mass in list(self.curr_xic_encoded_plot.keys())) and (mass_deviation == self.curr_mass_deviation):
                    print("XIC plot already calculated")
                    base64_image_1 = self.curr_xic_encoded_plot[xic_mass]
                else:
                    self.curr_xic_encoded_plot = {}
                    xic = MS_functions.get_xic_fast(self.curr_ms_file.rawdata, xic_mass, mass_deviation)
                    self.curr_mass_deviation = mass_deviation
                    print("XIC calculated")
                    fig, ax = plt.subplots(layout="tight")
                    ax.plot(xic[0], xic[1])
                    ax.set_title("XIC of mass " + str(xic_mass) + " in spectrum " + str(file_select))
                    ax.set_xlabel("Retention time (s)")
                    ax.set_ylabel("Intensity")
                    png_image_1 = io.BytesIO()
                    fig.savefig(png_image_1, format="png")
                    png_image_1.seek(0)
                    base64_image_1 = base64.b64encode(png_image_1.read()).decode()
                    self.curr_xic_encoded_plot[xic_mass] = base64_image_1
                print("XIC plot saved")

                if calc_new_file == False and (spec_index in list(self.curr_spec_encoded_plot.keys())):
                    print("Spectrum plot already calculated")
                    base64_image_2 = self.curr_spec_encoded_plot[spec_index]
                else:
                    self.curr_spec_encoded_plot = {}
                    fig2, ax2 = plt.subplots(layout="tight")
                    masses = list(self.curr_ms_file.rawdata[spec_index]["m/z array"])
                    intensities = list(self.curr_ms_file.rawdata[spec_index]["intensity array"])
                    ax2.bar(masses, intensities)
                    ax2.set_title("Spectrum " + str(spec_index) + " in file " + str(file_select))
                    ax2.set_xlabel("m/z")
                    ax2.set_ylabel("Intensity")
                    png_image_2 = io.BytesIO()
                    fig2.savefig(png_image_2, format="png")
                    png_image_2.seek(0)
                    base64_image_2 = base64.b64encode(png_image_2.read()).decode()
                    self.curr_spec_encoded_plot[spec_index] = base64_image_2
                print("Spectrum plot saved")
                return flask.jsonify({'plot1': base64_image_1, 'plot2': base64_image_2})
            return "test"

        @self.app.route("/all_routes")
        def all_routes():
            """Show links to all the available websites within this server."""
            urls = [str(rule) for rule in self.app.url_map.iter_rules() if rule.endpoint != 'static']
            return flask.render_template("all_routes.html", urls=urls)
        
        @self.app.route("/help_page", methods=["GET", "POST"])
        def help_page():
            if flask.request.method == "POST":
                question = flask.request.form.get("question")
                answer = flask.request.form.get("answer")
                if question != "" and answer != "":
                    with open(self.help_page_contents_filepath, "a") as f:
                        f.write("\nQ: " + question + "\nA: " + answer + "\n")
            
            if not os.path.exists(self.help_page_contents_filepath):
                with open(self.help_page_contents_filepath, "w") as f:
                    f.write("")
            with open(self.help_page_contents_filepath, "r") as f:
                lines = f.read().strip().split("\n")
            contents_list = []
            question, answer = None, None
            for line in lines:
                if line.startswith("Q: "):
                    question = line[3:].strip()
                elif line.startswith("A: "):
                    answer = line[3:].strip()
                    if question is not None and answer is not None:
                        contents_list.append((question, answer))
                        question, answer = None, None

            return flask.render_template("help_page.html", contents_list=contents_list)
        
        @self.app.route("/detect_peaks", methods=["GET", "POST"])
        def detect_peaks():
            # Take the mzml file that the user has uploaded and perform a complete peak detection on the file. Store the results in a folder on the webserver_save/results/mzml file folder.
            # The user can select the mzml file from a dropdown menu.
            available_mzml_files = os.listdir(self.mzml_folder)
            available_mzml_files = [f for f in available_mzml_files if f.endswith(".mzML")]
            if flask.request.method == "POST":
                settings_dict = self.get_settings_dict()
                mass_deviation = settings_dict["mass_deviation"]

                form_data = flask.request.form.to_dict()
                if not "fileselection" in form_data:
                    return flask.render_template_string("No file selected!")
                start_analysis_once_detected = form_data.get("start_analysis_once_detected", "false")

                print(form_data)
                if "fileselection" in form_data:
                    ms_filepath = self.mzml_folder + form_data["fileselection"]
                    parentfolder = self.parentfolder + "results/" + ".".join(form_data["fileselection"].split(".")[:-1]) + "/"
                    #UVenture_OA_starttime_20231001:12:00:00_MZ_
                    proc_name = "UVenture_peak_detection_starttime_" + str(datetime.datetime.now().strftime("%Y%m%d:%H%M%S")) + "_MZ__FILE_" + str(form_data["fileselection"])
                    print("Starting process: " + proc_name)
                    peaklist_filename = parentfolder + "peaklist.txt"
                    taskstorage = ""
                    if start_analysis_once_detected == "true" or start_analysis_once_detected == True:
                        #ms_filepath, parentfolder, peaklist_filename, mzrt_filename
                        taskstorage = self.taskstorage_filepath
                        p = multiprocessing.Process(target=start_pd_process, args=(ms_filepath, parentfolder, peaklist_filename, taskstorage, settings_dict), name=proc_name)
                    else:
                        p = multiprocessing.Process(target=start_pd_process, args=(ms_filepath, parentfolder, peaklist_filename, taskstorage, settings_dict), name=proc_name)
                    p.start()
                    self.running_processes.append(p)
                    return flask.redirect("/show_currently_running")

            return flask.render_template("detect_peaks.html", files=available_mzml_files)


        @self.app.route("/upload_mzml_file", methods=["POST", "GET"])
        def upload_mzml_file():
            """If the user has uploaded a file, save it to the parent folder"""
            if flask.request.method == "POST":
                # Check if the file is in the request
                if "file" in flask.request.files:
                    file = flask.request.files["file"]
                    filename = werkzeug.utils.secure_filename(file.filename)
                    form_data = flask.request.form.to_dict()
                    print(form_data)
                    uploader_info = {"client_ip": flask.request.remote_addr, "user_agent": flask.request.user_agent, "headers": flask.request.headers}
                    uploader_info["cookies"] = flask.request.cookies
                    uploader_info["url_params"] = flask.request.args
                    uploader_info["fullpath"] = flask.request.full_path
                    uploader_info["method"] = flask.request.method
                    uploader_info["referrer"] = flask.request.referrer
                    uploader_info["remote_user"] = flask.request.remote_user
                    uploader_info["scheme"] = flask.request.scheme
                    uploader_info["url"] = flask.request.url
                    uploader_info["url_root"] = flask.request.url_root
                    uploader_info["is_secure"] = flask.request.is_secure
                    uploader_info["host"] = flask.request.host
                    uploader_info["host_url"] = flask.request.host_url
                    uploader_info["base_url"] = flask.request.base_url
                    uploader_info["path"] = flask.request.path
                    uploader_info["mime_type"] = flask.request.mimetype
                
                    # Check if the file is a mzML file
                    if not filename.endswith(".mzML"):
                        return flask.render_template_string("File must be in mzML format!")

                    with open(self.mzml_folder + "file_metadata.txt", "a") as f:
                        f.write(str(datetime.datetime.now().strftime("%D/%m/%Y, %H:%M:%S")) + "\t")
                        f.write(str(filename) + "\t")
                        f.write(str(form_data) + "\t")
                        f.write(str(uploader_info) + "\n")
                        
                    #check if parentfolder exists. If not, create it
                    os.makedirs(self.mzml_folder, exist_ok=True)
                    #check if a file with the same name already exists
                    if filename in self.available_files:
                        return flask.render_template_string("File with the same name already exists!")
                    file.save(self.mzml_folder + filename)

                    self.available_files.append(filename)
                    return flask.redirect("/queue_new_analysis"), 200
            return flask.render_template("upload_mzml_file.html")

        @self.app.route("/file_browser", methods=["GET", "POST"])
        def file_browser():
            self.start_background_task_checking()
            try:
                available_space_left = shutil.disk_usage(self.parentfolder).free / (1024 * 1024 * 1024)  # Convert to GB
                available_space_left = str(round(available_space_left, 3)) + " GB"
            except:
                available_space_left = "NA"
            folder = resource_path("webserver_save/results/")
            #folder = os.path.normpath(folder)
            contents = []
            subfolders = []
            for item in os.listdir(folder):
                full_path = os.path.join(folder, item)
                full_path = full_path.replace(folder, "")
                full_path = "webserver_save/results/" + full_path
                if os.path.isdir(full_path):
                    subfolders.append(full_path.replace("/", "->"))
                else:
                    full_path = full_path.replace("webserver_save/results/", "")
                    full_path = "webserver_save->results->" + full_path
                    contents.append(full_path)
                
            return flask.render_template("file_browser.html", contents=contents, subfolders=subfolders, available_space_left=available_space_left)
        
        @self.app.route("/file_browser/delete/<folder>", methods=["GET", "POST"])
        def file_browser_delete(folder):
            folder = folder.replace("->", "/")
            folder = resource_path(folder)
            folder = os.path.normpath(folder)
            if not str(folder).startswith(str(os.path.normpath(self.results_folder))):
                print("Folder is not in the results folder. Redirecting to file browser.")
                print("Folder: " + str(folder))
                print("Results folder: " + str(self.results_folder))
                return flask.redirect("/file_browser")
            # Get the last part of the folder path
            folder = folder.replace(resource_path(""), "")
            upper_folder = "->".join(folder.replace("\\\\", "->").replace("\\", "->").replace("/", "->").split("->")[:-1])
            print(upper_folder)
            # Check if folder is a file or a directory
            folder = resource_path(folder)
            if os.path.isfile(folder):
                #Delete the file
                os.remove(folder)
                return flask.redirect("/file_browser/" + str(upper_folder))
            else:
                if os.path.exists(folder):
                    shutil.rmtree(folder)
                    return flask.redirect("/file_browser/" + str(upper_folder))
                else:
                    return flask.render_template_string("Folder does not exist!")

        @self.app.route("/file_browser/<folder>", methods=["GET", "POST"])
        def file_browser_folder(folder):
            try:
                available_space_left = shutil.disk_usage(self.parentfolder).free / (1024 * 1024 * 1024)  # Convert to GB
                available_space_left = str(round(available_space_left, 3)) + " GB"
            except:
                available_space_left = "NA"
            image_folder = folder.replace("->", "/") + "/"
            print(image_folder)
            folder = image_folder
            folder = resource_path(folder)
            folder = os.path.normpath(folder)
            print(folder)
            if not str(folder).startswith(str(os.path.normpath(self.results_folder))):
                print("Folder is not in the results folder. Redirecting to file browser.")
                print("Folder: " + str(folder))
                print("Results folder: " + str(self.results_folder))
                return flask.redirect("/file_browser")
            contents = []
            subfolders = []
            # Get the last part of the folder path
            folder = folder.replace(resource_path(""), "")
            upper_folder = "->".join(folder.replace("\\\\", "->").replace("\\", "->").replace("/", "->").split("->")[:-1])
            print(upper_folder)
            # Check if folder is a file or a directory
            folder = resource_path(folder)
            folder = os.path.normpath(folder)
            if os.path.isfile(folder):
                # Serve the file directly as a download
                filename = os.path.basename(folder)
                #Check if the file is an image or a txt file
                if filename.endswith(".png") or filename.endswith(".jpg") or filename.endswith(".jpeg") or filename.endswith(".tiff"):
                    return flask.send_file(folder)
                elif filename.endswith(".txt") or filename.endswith(".csv"):
                    # Show the contents of the text file in the browser
                    with open(folder, "r") as f:
                        contents = f.read()
                    return flask.render_template_string("<pre>" + contents + "</pre>")
                else:
                    # Serve the file as a download
                    return flask.send_file(folder, as_attachment=True, download_name=filename)
            else: #Folder is a directory
                # Check if the folder exists
                if not os.path.exists(image_folder):
                    return flask.render_template_string("Folder does not exist!")
                # Check if the folder is empty
                if len(os.listdir(folder)) == 0:
                    return flask.render_template_string("Folder is empty!")
                # Get the contents of the folder
                contents = []
                subfolders = []
                for item in os.listdir(folder):
                    full_path = os.path.join(folder, item)
                    if os.path.isdir(full_path):
                        subf = full_path.replace("/", "->")
                        subf = subf.replace("\\", "->")
                        subfolders.append(subf)

                    else:
                        cont = full_path.replace("/", "->")
                        cont = cont.replace("\\", "->")
                        contents.append(cont)
                # Sort the contents and subfolders
                contents.sort()
                subfolders.sort()
            print(contents, subfolders)
            return flask.render_template("file_browser.html", contents=contents, subfolders=subfolders, available_space_left=available_space_left)

        @self.app.route("/api/get_taskstorage", methods=["GET"])
        def get_taskstorage():
            with open(self.taskstorage_filepath, "r") as f:
                lines = f.readlines()
            lines = [line for line in lines if line.strip() != "" and line[0] != "#" and "\t" in line]
            return flask.jsonify(lines)
        
        @self.app.route("/api/delete_taskstorage", methods=["GET", "POST"])
        def delete_taskstorage():
            with open(self.taskstorage_filepath, "w") as f:
                f.write("")
            return flask.jsonify({"status": "success", "message": "Task storage file deleted"})
        
        @self.app.route("/api/kill_all_processes", methods=["GET", "POST"])
        def kill_processes():
            with open(self.taskstorage_filepath, "w") as f:
                f.write("")
            # Kill all running processes
            for p in self.running_processes:
                try:
                    print("Killing process: " + str(p))
                    try:
                        os.kill(p.pid, signal.SIGKILL)
                    except:
                        try:
                            os.kill(p.pid, signal.SIGILL)
                        except:
                            try:
                                os.kill(p.pid, signal.SIGINT)
                            except:
                                try:
                                    os.kill(p.pid, signal.SIGQUIT)
                                except:
                                    os.kill(p.pid, signal.SIGTERM)
                except Exception as e:
                    print(e)
                    print(traceback.format_exc())
            self.start_background_task_checking()
            return flask.jsonify({"status": "success", "message": "Processes killed"})

        @self.app.route("/api/get_settings", methods=["GET"])
        def get_settings():
            settings = self.get_settings_dict()
            return flask.jsonify(settings)
        
        @self.app.route("/api/stop_script", methods=["GET", "POST"])
        def stop_script():
            import sys
            sys.exit("Script stopped via API command...")
            return "Script stopped via API command..."

        @self.app.route("/api/get_available_files", methods=["GET"])
        def get_available_files():
            # Get the list of available files in the mzml folder
            self.available_files = os.listdir(self.mzml_folder)
            self.available_files = [f for f in self.available_files if f.endswith(".mzML")]
            return flask.jsonify(self.available_files)

        @self.app.route("/api/queue_analysis", methods=["GET", "POST"])
        def queue_analysis():
            mz = flask.request.args.get("mz", type=float)
            rt = flask.request.args.get("rt", type=float)
            settings_dict = self.get_settings_dict()
            filename = flask.request.args.get("filename", type=str)
            if mz is None or rt is None or filename is None:
                return flask.jsonify({"status": "error", "message": "Missing parameters"})
            with open(self.taskstorage_filepath, "a") as f:
                f.write(filename + "\t" + str(mz) + "\t" + str(rt) + "\t" + str(settings_dict) + "\n")
            print("Task added to task storage file")
            return flask.jsonify({"status": "success", "message": "Task added to task storage file"})
        
        @self.app.route("/api/get_status", methods=["GET"])
        def get_status():
            try:
                import psutil
                server_uptime = datetime.datetime.now() - datetime.datetime.fromtimestamp(psutil.boot_time())
                cpu_utilization = psutil.cpu_percent(interval=0.01)
                available_cpu_count = psutil.cpu_count(logical=False)
                available_cpu_count_logical = psutil.cpu_count(logical=True)
                try:
                    cpu_temperature = psutil.sensors_temperatures()
                except:
                    cpu_temperature = "NA"
                memory_utilization = psutil.virtual_memory().percent
                available_memory = str(round(psutil.virtual_memory().available / (1024 * 1024 * 1024), 3)) + " GB"
                cpu_utilization = psutil.cpu_percent(interval=0.3)
            except ImportError:
                server_uptime = datetime.datetime.now() - self.start_time
                cpu_utilization = "NA"
                available_cpu_count = "NA"
                available_cpu_count_logical = "NA"
                cpu_temperature = "NA"
                memory_utilization = "NA"
                available_memory = "NA"
            
            # Get the available files in the mzml folder
            self.available_files = os.listdir(self.mzml_folder)
            self.available_files = [f for f in self.available_files if f.endswith(".mzML")]
            total, used, free = shutil.disk_usage(self.parentfolder)
            used_percent = round((used / total) * 100, 2)
            # Get the number of pending analyses
            with open(self.taskstorage_filepath, "r") as f:
                lines = f.readlines()
            lines = [line for line in lines if line.strip() != "" and line[0] != "#" and "\t" in line]
            pending_analyses = len(lines)
            status = {
                "sha256_of_this_script": str(self.readable_hash),
                "last_modification_of_this_script": str(self.file_modification_time),
                "abspath_of_this_script": str(self.path_of_this_file),
                "disk_path": self.parentfolder,
                "version": str(self.version),
                "python_version": str(self.python_version),
                "server_uptime": str(round(server_uptime.days, 2)) + " days",
                "script_runtime": str(round((datetime.datetime.now() - self.start_time).days, 2)) + " days",
                "available_cores": str(self.num_cores_to_use),

                "cpu_utilization": str(cpu_utilization) + "%",
                "cpu_count": str(available_cpu_count),
                "cpu_count_logical": str(available_cpu_count_logical),
                "cpu_temperature": str(cpu_temperature),
                "memory_utilization": str(memory_utilization) + "%",
                "available_memory": available_memory,
                "disk_utilization": str(used_percent) + "%",
                "disk_total": str(round(total / (1024 * 1024 * 1024), 3)) + " GB",
                "disk_used": str(round(used / (1024 * 1024 * 1024), 3)) + " GB",
                "disk_free": str(round(free / (1024 * 1024 * 1024), 3)) + " GB",
                
                "available_files": self.available_files,
                "pending_analyses": pending_analyses,
            }
            return flask.jsonify(status)

    def get_settings_dict(self):
        settings_dict = {}
        with open(self.settings_file_filepath, "r") as f:
            file_contents_raw = f.read()
            lines = file_contents_raw.split("\n")
            for line in lines:
                line = line.strip()
                if not "=" in line:
                    continue
                if line[0] == "#":
                    continue
                key, value = line.split("=")
                try:
                    value = float(value)
                except:
                    try:
                        value = int(value)
                    except:
                        try:
                            value = ast.literal_eval(value)
                        except:
                            value = str(value)
                settings_dict[key] = value
        return settings_dict

    def start_background_task_checking(self):
        def task_check_bckg_func():
            while True:
                try:
                    # Check if a new task can be loaded
                    last_check_time = self.check_and_load_new_task()
                    # Check if the running_processes list is still up to date
                    self.update_running_processes()
                    time.sleep(0.1)
                except Exception as e:
                    print(e)
                    print(traceback.format_exc())
        #Check if the thread is already running
        all_threads = threading.enumerate()
        running_threads = [t.name for t in all_threads if t.is_alive()]
        running_processes = [p.name for p in multiprocessing.active_children()]
        running_analyses = running_threads + running_processes
        running_analyses = [a for a in running_analyses if "UVenture_task_check_thread_" in a]
        if len(running_analyses) > 0:
            print("Task checking thread already running. Not starting a new one.")
            return
        task_check_thread_name = "UVenture_task_check_thread_" + str(datetime.datetime.now().strftime("%Y%m%d:%H%M%S"))
        task_check_thread = threading.Thread(target=task_check_bckg_func, name=task_check_thread_name)
        task_check_thread.start()
        task_check_thread.join(0.1)  # Wait for the thread to start before returning
        print("Background task checking started")
        return

    def check_and_load_new_task(self):
        settings_dict = self.get_settings_dict()
        core_count = self.num_cores_to_use
        if core_count <= 1:
            core_count = 1
        all_threads = threading.enumerate()
        running_threads = [t.name for t in all_threads if t.is_alive()]
        running_processes = [p.name for p in multiprocessing.active_children()]
        running_analyses = running_threads + running_processes
        running_analyses = [a for a in running_analyses if "UVenture_" in a]
        #print("Running analyses: " + str(running_analyses))
        if len(running_analyses) <= core_count:
            #print("Checking for new tasks")
            pass
        else:
            #print("Too many tasks running. Not checking for new tasks.")
            return datetime.datetime.now()
        if not os.path.exists(self.taskstorage_filepath):
            print("Task storage file does not exist. Creating it.")
            with open(self.taskstorage_filepath, "w") as f:
                f.write("")
                return datetime.datetime.now()
        try:
            with open(self.taskstorage_filepath, "r") as f:
                lines = f.readlines()
            lines = [line for line in lines if line.strip() != "" and line[0] != "#" and "\t" in line]
            if len(lines) == 0:
                #print("No tasks to load")
                return datetime.datetime.now()
            
            try:
                ms_filepath = self.mzml_folder + lines[0].split("\t")[0]
                parentfolder_msfile = self.results_folder + str(".".join(lines[0].split("\t")[0].split(".")[:-1])) + "/"
                mz = lines[0].split("\t")[1]
                rt = lines[0].split("\t")[2]
                mz = float(mz)
                rt = float(rt)
                settings_dict = ast.literal_eval(lines[0].split("\t")[3])
                print("New task loaded: file: " + str(ms_filepath) + " mz: " + str(mz) + " rt: " + str(rt))
                self.start_one_oa(ms_filepath,
                                    mz,
                                    rt,
                                    settings_dict,
                                    parentfolder_msfile)
                print("Started new process for task: " + str(ms_filepath) + " mz: " + str(mz) + " rt: " + str(rt))
            except Exception as e:
                print(e)
                print(traceback.format_exc())

            # Remove the first line from the file
            print("Removing first line from task storage file")
            lines = lines[1:]
            with open(self.taskstorage_filepath, "w") as f:
                for line in lines:
                    line = line.strip()
                    if line == "" or line[0] == "#" or not "\t" in line or line == "\n":
                        continue
                    f.write(line + "\n")
            print("Removed first line from task storage file")

        except Exception as e:
            print(e)
            print(traceback.format_exc())
            return datetime.datetime.now()
        
        return datetime.datetime.now()
    
    def update_running_processes(self):
        for p in self.running_processes:
            if not p.is_alive():
                print("Process " + str(p.name) + " is not alive. Removing it from the list.")
                try:
                    # process name looks like: UVenture_OA_starttime_20231001:120000_MZ_150.0_RT_1.0_FILE_test.mzML
                    self.process_times[p.name] = [datetime.datetime.strptime(p.name.split("starttime_")[1].split("_MZ_")[0], "%Y%m%d:%H%M%S"), datetime.datetime.now(), p.name.split("_FILE_")[1]]
                    self.update_best_time_approx_per_analysis()
                    print("Process " + str(p.name) + " finished. Time taken: " + str(self.process_times[p.name][1] - self.process_times[p.name][0]))
                except Exception as e:
                    print(e)
                    print(traceback.format_exc())
                self.running_processes.remove(p)
    
    def update_best_time_approx_per_analysis(self):
        if len(self.process_times) == 0:
            self.best_time_approx_per_analysis = 0.001
            print("No process times available. Setting best time to 0.001")
            return
        times = []
        for k, v in self.process_times.items():
            print(k, v)           
            if len(v) == 3:
                times.append((v[1] - v[0]).total_seconds())
        print("Process times: " + str(times))
        average_time_seconds = sum(times)/len(times)
        average_time_div_by_cores = average_time_seconds / self.num_cores_to_use
        self.best_time_approx_per_analysis = average_time_div_by_cores
        print("Best time approximation per analysis: " + str(self.best_time_approx_per_analysis))
        return self.best_time_approx_per_analysis

    def start_one_oa(self, ms_filepath, mz, rt, settings_dict, parentfolder_msfile):
        mz = float(mz)
        rt = float(rt)
        thread_name = ("UVenture_OA_starttime_" + str(datetime.datetime.now().strftime("%Y%m%d:%H%M%S")) + "_MZ_" + str(mz) + "_RT_" + str(rt) + "_FILE_" + str(ms_filepath.split("/")[-1]))
        print("Starting new process: " + thread_name)
        proc = multiprocessing.Process(target=caller_func, args=[ms_filepath, mz, rt, settings_dict, parentfolder_msfile], name=thread_name)
        print("Process created: " + str(proc))
        proc.start()
        print("Process started: " + str(proc))
        self.running_processes.append(proc)
        print("Process added to running processes list")
        self.update_running_processes()
        





if __name__ == "__main__":
    multiprocessing.freeze_support()
    webapp = Webpage()
    #webapp.app.run(host="0.0.0.0", debug=False, use_reloader=False, port=5000)
    webapp.app.run(debug=False, use_reloader=False, port=5000)
    print("Server running")
    while True:
        time.sleep(1)
    
