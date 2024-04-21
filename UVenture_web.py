import os
import threading
import time
import traceback
import flask
import werkzeug
import UVenture


class Webpage:
    def __init__(self) -> None:
        self.parentfolder = "webserver_save/"
        self.mzml_folder = self.parentfolder + "mzml_files/"
        self.results_folder = self.parentfolder + "results/"
        self.app = flask.Flask(__name__)
        self.server = None
        self.available_files = os.listdir(self.mzml_folder)

        @self.app.route("/", methods=["GET", "POST"])
        def index():
            return flask.redirect("/queue_new_analysis")
        
        @self.app.route("/queue_new_analysis", methods=["GET", "POST"])
        def queue_new_analysis():
            if flask.request.method == "POST":
                form_data = flask.request.form.to_dict()
                if not "peak_analysis_cb" in form_data:
                    form_data["peak_analysis_cb"] = "false"
                if not "mass_analysis_cb" in form_data:
                    form_data["mass_analysis_cb"] = "false"
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
                ms_filepath = self.mzml_folder + form_data["fileselection"]
                ms_file = UVenture.MS_File(ms_filepath, parentfolder=self.results_folder + str(".".join(form_data["fileselection"].split(".")[:-1])) + "/")
                print("MS file loaded")
                if form_data["peak_analysis_cb"] == "true":
                    if form_data["mz_peak_analysis"] == "":
                        return flask.render_template_string("No m/z given! Cannot analyze peak without a mass. \n Please enter a peak to analyse")
                    if form_data["retention_time"] == "":
                        return flask.render_template_string("No retention time given! Cannot analyze peak without retention time. \nPlease enter a retention time")
                    try:
                        myanalysis = UVenture.OneAnalysis(ms_file, form_data["mz_peak_analysis"], form_data["retention_time"], mass_deviation=11, mass_deviation_xic=0.001, charge_of_measured_mass=-1, formula_cache_folder_path="C://Users//Admin//Desktop//UVenture//Formula_Predictions//Formula_Predictions//")
                    except Exception as e:
                        print(e)
                        print(traceback.format_exc())
                        return flask.render_template_string("An error occured during the analysis" + str(e) + "\n \n \n" + str(traceback.format_exc()))
                    print("Peak analysis started")
                if form_data["mass_analysis_cb"] == "true":
                    myspec = UVenture.Spec(ms_file, form_data["spec_index"], requested_filter_mode="whatever", save_plot=True)
                    myanalysis = UVenture.Prediction(ms_file, form_data["mz_mass_analysis"], myspec, save_plot_of_isotopologues=True, charge_of_measured_mass=-1, formula_cache_folder_path="C://Users//Admin//Desktop//UVenture//Formula_Predictions//Formula_Predictions//")
                    print("Mass analysis started")
                print("Analysis started")
                return flask.redirect("/")
            return flask.render_template("queue_new_analysis.html", files=self.available_files)


        @self.app.route("/all_routes")
        def all_routes():
            """Show links to all the available websites within this server."""
            urls = [str(rule) for rule in self.app.url_map.iter_rules() if rule.endpoint != 'static']
            return flask.render_template("all_routes.html", urls=urls)
            
        @self.app.route("/analysis_started", methods=["POST", "GET"])
        def analysis_started():
            return flask.render_template("analysis_started.html")
        
        @self.app.route("/upload_mzml_file", methods=["POST", "GET"])
        def upload():
            """If the user has uploaded a file, save it to the parent folder"""
            if "file" in flask.request.files:
                file = flask.request.files["file"]
                filename = werkzeug.utils.secure_filename(file.filename)
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
                return flask.redirect("/queue_new_analysis")
            return flask.render_template("upload_mzml_file.html")

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





if __name__ == "__main__":
    webapp = Webpage()
    webapp.app.run(debug=True, use_reloader=True, port=5000)
    print("Server running")
    while True:
        time.sleep(1)
    