import ast
import base64
import datetime
import io
import os
import shutil
import subprocess
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
import psutil


def caller_func(ms_filepath, mz, rt, settings_dict, parentfolder_msfile):
    try:
        ms_file = UVenture.MS_File(ms_filepath, parentfolder_msfile=parentfolder_msfile)
        myanalysis = UVenture.OneAnalysis(ms_file, mz, rt, **settings_dict)
        return
    except Exception as e:
        print(e)
        print(traceback.format_exc())
        return


class Webpage:
    def __init__(self) -> None:
        self.parentfolder = "webserver_save/"
        os.makedirs(self.parentfolder, exist_ok=True)
        self.mzml_folder = self.parentfolder + "mzml_files/"
        os.makedirs(self.mzml_folder, exist_ok=True)
        self.results_folder = self.parentfolder + "results/"
        os.makedirs(self.results_folder, exist_ok=True)
        self.app = flask.Flask(__name__)
        self.server = None
        self.available_files = os.listdir(self.mzml_folder)
        self.available_files = [f for f in self.available_files if f.endswith(".mzML")]
        self.curr_ms_file = None
        self.curr_xic_encoded_plot = {}
        self.curr_spec_encoded_plot = {}
        self.curr_mass_deviation = 0
        self.settings_file_filepath = "static/settings.txt"
        self.settings_default_file_filepath = "static/settings_default.txt"
        self.help_page_contents_filepath = "static/help_page_contents.txt"

        @self.app.route("/", methods=["GET", "POST"])
        def index():
            if len(self.available_files) == 0:
                return flask.redirect("/upload_mzml_file")
            else:
                return flask.redirect("/queue_new_analysis")
        
        @self.app.route("/queue_new_analysis", methods=["GET", "POST"])
        def queue_new_analysis():
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
                        self.start_one_oa(self.mzml_folder + form_data["fileselection"],
                                          mz=form_data["mz_peak_analysis"],
                                          rt=form_data["retention_time"],
                                          settings_dict=settings_dict,
                                          parentfolder_msfile=self.results_folder + str(".".join(form_data["fileselection"].split(".")[:-1])) + "/")
                        all_threads = threading.enumerate()
                        running_threads = [t.name for t in all_threads if t.is_alive()]
                        running_processes = [p.name for p in multiprocessing.active_children()]
                        print("Running threads:", running_threads)
                        print("Running processes:", running_processes)
                        ...
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
                    print(lines)

                    def create_task(lines):
                        core_count = os.cpu_count()
                        core_count = 6

                        ms_filepath = self.mzml_folder + form_data["fileselection"]
                        parentfolder_msfile = self.results_folder + str(".".join(form_data["fileselection"].split(".")[:-1])) + "/"

                        for line in lines:
                            line.strip()
                            print(line)
                            mz, rt = line.split("\t")
                            try:
                                mz = float(mz)
                                rt = float(rt)
                            except:
                                continue

                            all_threads = threading.enumerate()
                            running_threads = [t.name for t in all_threads if t.is_alive()]
                            running_processes = [p.name for p in multiprocessing.active_children()]
                            running_analyses = running_threads + running_processes
                            running_analyses = [a for a in running_analyses if "UVenture_" in a]
                            print("Running analyses: " + str(running_analyses))

                            while len(running_analyses) >= core_count:
                                all_threads = threading.enumerate()
                                running_threads = [t.name for t in all_threads if t.is_alive()]
                                running_processes = [p.name for p in multiprocessing.active_children()]
                                running_analyses = running_threads + running_processes
                                running_analyses = [a for a in running_analyses if "UVenture_" in a]
                                time.sleep(1)

                            try:
                                print("scheduling")
                                self.start_one_oa(ms_filepath,
                                                  mz,
                                                  rt,
                                                  settings_dict,
                                                  parentfolder_msfile)
                                print("Scheduled new one analysis")
                            except Exception as e:
                                print(e)
                                print(traceback.format_exc())
                                continue


                    #create_task(lines)
                    thread_name = "UVenture_ListAnalysis_starttime_" + str(datetime.datetime.now().strftime("%Y%m%d:%H%M%S")) + "_AnalysisList_" + str(form_data["fileselection"]) + "WHOLELIST"
                    thread = threading.Thread(target=create_task, args=[lines], name=thread_name)
                    thread.start()
                    print("List analysis started")
                
                print("Analysis started")
                return flask.redirect("/show_currently_running")
            return flask.render_template("queue_new_analysis.html", files=self.available_files)

        @self.app.route("/change_settings", methods=["GET", "POST"])
        def change_settings():
            settings_dict = {}
            file_contents_raw = ""
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
                
            return flask.render_template("change_settings.html", settings_dict=settings_dict, file_contents_raw=file_contents_raw)

        @self.app.route("/resultsdownload", methods=["GET", "POST"])
        def resultsdownload():
            if flask.request.method == "POST":
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
            all_threads = threading.enumerate()
            running_threads = [t.name for t in all_threads if t.is_alive()]
            running_processes = [p.name for p in multiprocessing.active_children()]
            running_analyses = running_threads + running_processes
            running_analyses = [a for a in running_analyses if "UVenture_" in a]

            tasks_information = []
            for task in running_analyses:
                tasks_information.append({"name": task, "is_alive": True})
            tasks_information.sort(key=lambda x: x["name"])
            tasks_information.append({"name": "█████████████████████████████████████████████████████████████████████", "is_alive": True})
            for task in all_threads:
                if not task.name.startswith("UVenture_"):
                    tasks_information.append({"name": task.name, "is_alive": task.is_alive()})                
            return flask.render_template("show_currently_running.html", tasks_information=tasks_information)

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
                    self.curr_ms_file = UVenture.MS_File(ms_filepath, parentfolder=self.results_folder + str(".".join(file_select.split(".")[:-1])) + "/")
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
            contents_dict = {}
            with open(self.help_page_contents_filepath, "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    line = line.replace("\\n", "\n")
                    if line == "":
                        continue
                    if line[0] == "#":
                        continue

                    line = line.split("   esUWFds!&$!=====asfeDVW   ")
                    question = line[0]
                    answer = line[1]
                    contents_dict[question] = answer
            if flask.request.method == "POST":
                print(flask.request.form.to_dict())
                form_data = flask.request.form.to_dict()
                if "add_question" in form_data:
                    question = form_data["question"].replace("\n", "\\n")
                    answer = form_data["answer"].replace("\n", "\\n")
                    with open(self.help_page_contents_filepath, "a") as f:
                        f.write(str(question) + "   esUWFds!&$!=====asfeDVW   " + str(answer) + "\n")
                    contents_dict[question] = answer
                if "delete_question" in form_data:
                    question = form_data["question"].replace("\n", "\\n")
                    with open(self.help_page_contents_filepath, "r") as f:
                        lines = f.readlines()
                    with open(self.help_page_contents_filepath, "w") as f:
                        for line in lines:
                            if question in line:
                                continue
                            f.write(line)
                    contents_dict.pop(question)
            return flask.render_template("help_page.html", contents_dict=contents_dict)
        
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
                
                    with open(self.mzml_folder + "file_metadata.txt", "a") as f:
                        f.write(str(datetime.datetime.now().strftime("%D/%m/%Y, %H:%M:%S")) + "\t")
                        f.write(str(filename) + "\t")
                        f.write(str(form_data) + "\t")
                        f.write(str(uploader_info) + "\n")
                        
                    # check if the file ends with .mzML
                    if not filename.endswith(".mzML"):
                        return flask.render_template_string("File must be in mzML format!")
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
            folder = "webserver_save/results/"
            contents = []
            subfolders = []
            for item in os.listdir(folder):
                full_path = os.path.join(folder, item)
                if os.path.isdir(full_path):
                    subfolders.append(full_path.replace("/", "->"))
                else:
                    contents.append(full_path)
                
            return flask.render_template("file_browser.html", contents=contents, subfolders=subfolders)
        
        @self.app.route("/file_browser/<folder>", methods=["GET", "POST"])
        def file_browser_folder(folder):
            image_folder = folder.replace("->", "/") + "/"
            print(image_folder)
            folder = image_folder
            contents = []
            subfolders = []
            for item in os.listdir(folder):
                full_path = os.path.join(folder, item)
                if os.path.isdir(full_path):
                    subfolders.append(full_path.replace("/", "->"))
                else:
                    contents.append(full_path)
            print(contents, subfolders)
            return flask.render_template("file_browser.html", contents=contents, subfolders=subfolders)

        # Create a new thread for running the Flask application
        #self.server = werkzeug.serving.make_server('127.0.0.1', 5000, self.app)
        #self.thread = threading.Thread(target=self.server.serve_forever)
        #self.thread.start()
        #self.app.run(debug=True, use_reloader=True, port=5000)

    def __del__(self):
        # Stop the Flask application when the object is deleted
        if self.server is not None:
            self.server.shutdown()
        if self.thread is not None:
            self.thread.join()
        print("Webpage object deleted and server stopped")

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


    def start_one_oa_thread(self, ms_filepath, mz, rt, settings_dict, parentfolder_msfile):
        mz = float(mz)
        rt = float(rt)
        def caller_func(ms_filepath, mz, rt, settings_dict, parentfolder_msfile):
            try:
                ms_file = UVenture.MS_File(ms_filepath, parentfolder_msfile=parentfolder_msfile)
                myanalysis = UVenture.OneAnalysis(ms_file, mz, rt, **settings_dict)
                return
            except Exception as e:
                print(e)
                print(traceback.format_exc())
                return
        thread_name = ("UVenture_OA_starttime_" + str(datetime.datetime.now().strftime("%Y%m%d:%H%M%S")) + "_MZ_" + str(mz) + "_RT_" + str(rt))
        thread = threading.Thread(target=caller_func, args=[ms_filepath, mz, rt, settings_dict, parentfolder_msfile], name=thread_name)
        thread.start()

    def start_one_oa(self, ms_filepath, mz, rt, settings_dict, parentfolder_msfile):
        mz = float(mz)
        rt = float(rt)
        thread_name = ("UVenture_OA_starttime_" + str(datetime.datetime.now().strftime("%Y%m%d:%H%M%S")) + "_MZ_" + str(mz) + "_RT_" + str(rt))
        thread = multiprocessing.Process(target=caller_func, args=[ms_filepath, mz, rt, settings_dict, parentfolder_msfile], name=thread_name)
        thread.start()





if __name__ == "__main__":
    webapp = Webpage()
    webapp.app.run(debug=False, use_reloader=False, port=5000)
    print("Server running")
    while True:
        time.sleep(1)
    
