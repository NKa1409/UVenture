import MS_functions
import pyteomics
import flask
import time
import threading
import werkzeug
import os

class ms_file:
    def __init__(self, filename) -> None:
        self.file = None
        self.filename = filename
        self.rawdata = None
        self.all_ms_spectra = []
        self.all_filters = []
        self.method_duration = None
        self.available_modes = []
        self.tic = []
        self.open_mzml_file()
        self.xic_spectra = []
        self.ms_spectra = []
    
    def __call__(self):
        return self.rawdata
    
    def __repr__(self):
        return str(self.filename)
    
    def open_mzml_file(self):
        from pyteomics import mzml
        self.file = mzml.read(self.filename)
        self.rawdata = list(self.file)
        for index in range(len(self.rawdata)):
            masses = list(self.rawdata[index]["m/z array"])
            intensities = list(self.rawdata[index]["intensity array"])
            self.all_ms_spectra.append({"masses":masses, "intensities":intensities})
            self.tic.append(self.rawdata[index]["total ion current"])
            self.all_filters.append(self.rawdata[index]["scanList"]["scan"][0]["filter string"])
        for filter_string in self.all_filters:
            if " d " in filter_string and "@hcd" in filter_string:
                self.available_modes.append("MS/MS")
            elif " d " not in filter_string and "hcd" not in filter_string:
                self.available_modes.append("Full scan")
            elif " d " not in filter_string and "hcd" in filter_string:
                self.available_modes.append("AIF")
        self.available_modes = list(set(self.available_modes))
        self.method_duration = self.rawdata[-1]["scanList"]["scan"][0]["scan time"]
        return self
        



class Webpage:
    def __init__(self) -> None:
        self.parentfolder = None
        self.app = flask.Flask(__name__)
        self.server = None

        @self.app.route("/")
        def index():
            """Show links to all the available websites within this server."""
            urls = [str(rule) for rule in self.app.url_map.iter_rules() if rule.endpoint != 'static']
            return flask.render_template("index.html", urls=urls)
            
            
        
        @self.app.route("/upload_mzml_file", methods=["POST", "GET"])
        def upload():
            self.parentfolder = "test"
            """If the user has uploaded a file, save it to the parent folder"""
            if "file" in flask.request.files:
                file = flask.request.files["file"]
                filename = werkzeug.utils.secure_filename(file.filename)
                #check if parentfolder exists. If not, create it
                os.makedirs(self.parentfolder, exist_ok=True)
                file.save(self.parentfolder + "/" + filename)
                return flask.redirect("/upload_mzml_file")
            return flask.render_template("upload_mzml_file.html")

        # Create a new thread for running the Flask application
        self.server = werkzeug.serving.make_server('127.0.0.1', 5000, self.app)
        self.thread = threading.Thread(target=self.server.serve_forever)
        self.thread.start()

    def __del__(self):
        # Stop the Flask application when the object is deleted
        if self.server is not None:
            self.server.shutdown()
        if self.thread is not None:
            self.thread.join()
        print("Webpage object deleted and server stopped")


        
        
            
            



#mzml_filename = "F12_HRAIF4_3.mzML"
#ms_file = ms_file(mzml_filename)

test = Webpage()






