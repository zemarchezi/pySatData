import json
import pysatdata
class openConfigFile():
    def __init__(self):
        pathstring = pysatdata.__file__.replace("\\", "/")
        path = "/".join(pathstring.split("/")[:-1])

        path = path + "/resources/config_file.json"
        with open(path, 'r') as f:
            self.config_file_sat = json.load(f)