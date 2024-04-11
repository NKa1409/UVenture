import os
import threading
import time

import flask
import werkzeug


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





webapp = Webpage()
print("Server running")
while True:
    time.sleep(1)
    