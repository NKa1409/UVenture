import os
import threading
import time
import flask
import werkzeug


class Webpage:
    def __init__(self) -> None:
        self.parentfolder = "webserver_save/"
        self.mzml_folder = self.parentfolder + "mzml_files/"
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
        self.app.run(debug=True, use_reloader=True, port=5000)

    def __del__(self):
        # Stop the Flask application when the object is deleted
        if self.server is not None:
            self.server.shutdown()
        if self.thread is not None:
            self.thread.join()
        print("Webpage object deleted and server stopped")





webapp = Webpage()
print("Server running")
while True:
    time.sleep(1)
    