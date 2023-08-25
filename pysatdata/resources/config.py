import json
import pysatdata
class openConfigFile():
    def __init__(self):
        path = "/".join(pysatdata.__file__.split("/")[:-1])

        path = path + "/resources/config_file.json"
        with open(path, 'r') as f:
            self.config_file_sat = json.load(f)