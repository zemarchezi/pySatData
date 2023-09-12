import os

VERSION_PATH = os.path.join(os.path.dirname(__file__), "VERSION.txt")

with open(VERSION_PATH, "r") as version_file:
  VERSION = version_file.read().strip()

__version__ = VERSION
__license__ = "MIT"
__author__ = "José Paulo Marchezi"
__status__ = "Development"