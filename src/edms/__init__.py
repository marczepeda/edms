"""
src/edms                        EDMS package    
├── __init__.py                 Initializer    
├── __main__.py                 Main entry point    
├── _version.py                 Version
├── main.py                     Main module
├── utils/                      Utilities
├── config.py                   Configuration
├── resources/                  Resources
├── gen/                        General module
├── bio/                        Biology module 
├── dat/                        Database module
├── scripts/                    Scripts directory
└── resources/                  Resources directory
"""
# Load Arial font for matplotlib
from importlib.resources import files
from matplotlib import font_manager
import matplotlib as mpl

ttf = files("edms").joinpath("resources/Arial.ttf")
font_manager.fontManager.addfont(str(ttf))

family = font_manager.FontProperties(fname=str(ttf)).get_name()

mpl.rcParams["font.family"] = family
mpl.rcParams["mathtext.fontset"] = "custom"
mpl.rcParams["mathtext.rm"] = family
mpl.rcParams["mathtext.it"] = family
mpl.rcParams["mathtext.bf"] = family