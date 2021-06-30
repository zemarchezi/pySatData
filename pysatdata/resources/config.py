import json

class openConfigFile():
    def __init__(self, path):
        with open(path, 'r') as f:
            self.config_file_sat = json.load(f)