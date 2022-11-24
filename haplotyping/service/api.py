import os,sys,logging,configparser
from multiprocessing import Process,get_start_method
from flask import Flask, Blueprint, Response, render_template, current_app, g, request
from flask_restx import Api, Resource
from werkzeug.middleware.proxy_fix import ProxyFix
import json, sqlite3
from waitress import serve

import haplotyping

from haplotyping.service.api_tools import namespace as ns_api_tools
from haplotyping.service.api_country import namespace as ns_api_country
from haplotyping.service.api_collection import namespace as ns_api_collection
from haplotyping.service.api_variety import namespace as ns_api_variety
from haplotyping.service.api_dataset import namespace as ns_api_dataset
from haplotyping.service.api_kmer import namespace as ns_api_kmer
from haplotyping.service.api_split import namespace as ns_api_split
from haplotyping.service.api_marker import namespace as ns_api_marker

from haplotyping.service.api_kmer import cache as cache_api_kmer
from haplotyping.service.api_split import cache as cache_api_split

class API:
    
    def __init__(self, location):
                                
        self.location = str(location)
        
        #set logging
        logger_server = logging.getLogger(__name__+".server")
        
        self.config = configparser.ConfigParser()
        self.config.read(os.path.join(self.location,"server.ini")) 
        logger_server.info("read configuration file") 
        if self.config.getboolean("api","debug"):
            logger_server.info("run in debug mode") 
            logger_server.setLevel(logging.DEBUG)
        else:
            logger_server.setLevel(logging.INFO)
        
        #solve reload problem when using spawn method (osx/windows)
        if get_start_method()=="spawn":
            frame = sys._getframe()
            while frame:
                if "__name__" in frame.f_locals.keys():
                    if not frame.f_locals["__name__"]=="__main__":
                        return                    
                frame = frame.f_back
                
        #restart on errors
        while True:
            try:
                process_api = Process(target=self.process_api_messages, args=[])
                #start everything
                logger_server.info("start server on port {:s}".format(self.config["api"].get("port","8080")))
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
    
    def get_kmc_query_library():
        try:
            config = current_app.config.get("config")
            return config["settings"]["kmc_query_library"]
        except:
            return None
    
    def get_kmc_query_binary_location():
        try:
            config = current_app.config.get("config")
            return config["settings"]["kmc_query_binary_location"]
        except:
            return None
        
    def get_data_kmer_location():
        try:
            config = current_app.config.get("config")
            return config["settings"]["data_location_kmer"]
        except:
            return None
                                
    def get_data_marker_location():
        try:
            config = current_app.config.get("config")
            return config["settings"]["data_location_marker"]
        except:
            return None
                                
    def process_api_messages(self):    
        
        #--- initialize Flask application ---  
        logging.getLogger("werkzeug").disabled = True
        #os.environ["WERKZEUG_RUN_MAIN"] = "true"
        app = Flask(__name__, static_url_path="/static", 
                    static_folder=os.path.join(self.location,"static"), 
                    template_folder=os.path.join(self.location,"templates")) 
        if self.config.getboolean("api","proxy_prefix"):
            app.wsgi_app = ProxyFix(app.wsgi_app, x_prefix=1)
        app.config["config"] = self.config
        app.config["location"] = self.location
                
        #--- blueprint ---      
        blueprint = Blueprint("api", __name__, url_prefix="/api")
        api = Api(blueprint)
        
        #logger
        logger_api = logging.getLogger(__name__)
        if self.config.getboolean("api","debug"):
            logger_api.setLevel(logging.DEBUG)
        else:
            logger_api.setLevel(logging.INFO)

        #cache
        cache_config = {
            "CACHE_TYPE": "NullCache",
            "CACHE_NO_NULL_WARNING": True
        }
        if ("cache" in app.config["config"]) and ("type" in app.config["config"]["cache"]):
            if app.config["config"]["cache"]["type"]=="SimpleCache":
                logger_api.debug("caching in memory")
                cache_config = {
                    "CACHE_TYPE": "SimpleCache"
                }
                if app.config["config"]["cache"]["timeout"]:
                    cache_config["CACHE_DEFAULT_TIMEOUT"] = int(app.config["config"]["cache"]["timeout"])
                if app.config["config"]["cache"]["threshold"]:
                    cache_config["CACHE_THRESHOLD"] = int(app.config["config"]["cache"]["threshold"])      
                cache_api_kmer.init_app(app, config=cache_config)
                cache_api_split.init_app(app, config=cache_config)
            elif (app.config["config"]["cache"]["type"]=="FileSystemCache") and ("dir" in app.config["config"]["cache"]):       
                logger_api.debug("caching on disk: "+str(app.config["config"]["cache"]["dir"]))
                cache_config = {
                    "CACHE_TYPE": "FileSystemCache",
                    "CACHE_DIR": app.config["config"]["cache"]["dir"]
                }
                if ("timeout" in app.config["config"]["cache"]) and app.config["config"]["cache"]["timeout"]:
                    cache_config["CACHE_DEFAULT_TIMEOUT"] = int(app.config["config"]["cache"]["timeout"])
                if ("threshold" in app.config["config"]["cache"]) and app.config["config"]["cache"]["threshold"]:
                    cache_config["CACHE_THRESHOLD"] = int(app.config["config"]["cache"]["threshold"]) 
                #set specific caches    
                cache_config_kmer = cache_config
                cache_config_kmer["CACHE_DIR"] = os.path.join(app.config["config"]["cache"]["dir"],"kmer")
                cache_api_kmer.init_app(app, config=cache_config_kmer)
                cache_config_split = cache_config
                cache_config_split["CACHE_DIR"] = os.path.join(app.config["config"]["cache"]["dir"],"split")
                cache_api_split.init_app(app, config=cache_config_split)                
            else:
                logger_api.debug("caching disabled")
                cache_api_kmer.init_app(app, config=cache_config)
                cache_api_split.init_app(app, config=cache_config)
        else:
            cache_api_kmer.init_app(app, config=cache_config)
            cache_api_split.init_app(app, config=cache_config)
                            
            
    
        #namespaces
        api.add_namespace(ns_api_tools)
        api.add_namespace(ns_api_country)
        api.add_namespace(ns_api_collection)
        api.add_namespace(ns_api_variety)
        api.add_namespace(ns_api_dataset)
        api.add_namespace(ns_api_kmer)
        api.add_namespace(ns_api_split)
        api.add_namespace(ns_api_marker)
        
        app.register_blueprint(blueprint) 
        app.config.SWAGGER_UI_DOC_EXPANSION = "list"

        @app.before_request
        def before_request_func():
            app.logger.info("Request to "+str(request.endpoint)+" ("+str(request.content_length)+")")
            if self.config.getboolean("api","proxy_prefix"):
                if api._schema:
                    api._schema["basePath"] = api.base_path
                    api.__schema__["basePath"] = api.base_path

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
        #app.run(host=self.config["api"].get("host", "::"), port=self.config["api"].get("port", "8080"), 
        #        debug=self.config.getboolean("api","debug"), 
        #        use_reloader=False)   
        serve(app, 
              host=self.config["api"].get("host", "::"), 
              port=self.config["api"].get("port", "8080"), 
              threads=self.config["api"].get("threads", "10"))                            
        
        
        