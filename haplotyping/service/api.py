import os,logging,configparser
from multiprocessing import Process
from flask import Flask, Blueprint, Response, render_template, current_app, g
from flask_restx import Api, Resource
import json, sqlite3

import haplotyping

from haplotyping.service.api_tools import namespace as ns_api_tools
from haplotyping.service.api_country import namespace as ns_api_country
from haplotyping.service.api_collection import namespace as ns_api_collection
from haplotyping.service.api_variety import namespace as ns_api_variety
from haplotyping.service.api_dataset import namespace as ns_api_dataset
from haplotyping.service.api_kmer import namespace as ns_api_kmer

class API:
    
    def __init__(self, location, doStart=True):
        
        self.location = str(location)
        
        #set logging
        logger_server = logging.getLogger(__name__+".server")
        
        self.config = configparser.ConfigParser()
        self.config.read(os.path.join(self.location,"server.ini")) 
        logger_server.info("read configuration file") 
        if self.config.getboolean("api","debug"):
            logger_server.info("run in debug mode") 
        
        #restart on errors
        while doStart:
            try:
                process_api = Process(target=self.process_api_messages, args=[])
                #start everything
                logger_server.info("start server on port {:s}".format(self.config["api"]["port"]))
                process_api.start()
                #wait until ends  
                process_api.join()
            except Exception as e:  
                logger_server.error("error: "+ str(e))   
                
    def get_db_connection():
        db_connection = getattr(g, "_database", None)
        if db_connection is None:
            config = current_app.config.get("config")
            location = current_app.config.get("location")
            db_connection = g._database = sqlite3.connect(os.path.join(location,config["settings"]["sqlite_db"]))
        return db_connection
                                
    def process_api_messages(self, doStart=True):    
        
        #--- initialize Flask application ---  
        logging.getLogger("werkzeug").disabled = True
        os.environ["WERKZEUG_RUN_MAIN"] = "true"
        app = Flask(__name__, static_url_path="/static", static_folder=os.path.join(self.location,"static"), 
                    template_folder=os.path.join(self.location,"templates"))  
        app.config["config"] = self.config
        app.config["location"] = self.location

        #--- blueprint ---      
        blueprint = Blueprint("api", __name__, url_prefix="/api")
        api = Api(blueprint)

        #namespaces
        api.add_namespace(ns_api_tools)
        api.add_namespace(ns_api_country)
        api.add_namespace(ns_api_collection)
        api.add_namespace(ns_api_variety)
        api.add_namespace(ns_api_dataset)
        api.add_namespace(ns_api_kmer)
        
        app.register_blueprint(blueprint) 
        app.config.SWAGGER_UI_DOC_EXPANSION = "list"

        logger_api = logging.getLogger(__name__)
        parser = api.parser()
        
        #--- database ---    
        @app.teardown_appcontext
        def close_connection(exception):
            db = getattr(g, "_database", None)
            if db is not None:
                logger_api.debug("close database connection")
                db.close()

        #--- site ---
        @app.route("/")
        def index():
            return render_template("index.html")
        
        #--- start webserver ---
        if doStart:
            app.run(host=self.config["api"]["host"], port=self.config["api"]["port"], 
                    debug=self.config.getboolean("api","debug"), 
                    use_reloader=False)   
        else:
            return app
        
        
        