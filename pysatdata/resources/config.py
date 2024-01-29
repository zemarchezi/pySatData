import json
import pysatdata
class openConfigFile():
    def __init__(self, satellite):
        self.satellite = satellite
        pathstring = pysatdata.__file__.replace("\\", "/")
        path = "/".join(pathstring.split("/")[:-1])

        path = f"{path}/resources/{self.satellite}/config_file.json"
        print("loading config file from: ", path)
        with open(path, 'r') as f:
            self.config_file_sat = json.load(f)