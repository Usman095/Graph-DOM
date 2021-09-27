from configparser import ConfigParser
import os
import ast
from collections import OrderedDict


DEFAULT_PARAM_PATH = './config.ini'
PARAM_PATH = None

config = None


def get_config(section='input', key=None):
    """Read the configuration parameters and return a dictionary."""

    global config

    # If file path is given use it otherwise use default.
    file_path = PARAM_PATH if PARAM_PATH else DEFAULT_PARAM_PATH
    
    # Read config and convert each value to appropriate type.
    # Only for the first time.
    if not config:
        config = dict()
        config_ = ConfigParser()
        assert isinstance(file_path, str)
        config_.read(file_path)
        for section_ in config_.sections():
            config[section_] = dict()
            for key_ in config_[section_]:
                try:
                    config[section_][key_] = ast.literal_eval(config_[section_][key_])
                except (ValueError, SyntaxError):
                    config[section_][key_] = config_[section_][key_]

    if section and section in config:
        if key and key in config[section]:
            return config[section][key]
        return config[section]
    return config

elements = {"C":12000, "H":1008, "O":15995, "N":14003, "S":31972}
use_ns = get_config(section='params', key='use_NS')
if use_ns:
    atom_list = ['C', 'H', 'O', 'N', 'S']
else:
    atom_list = ['C', 'H', 'O']

#the order of neutral losses will be the same in the final output
neutral_losses = OrderedDict([("O", 15995), ("CH4", 16031),
                              ("H2O", 18010), ("CO", 27995), ("CH2O", 30011), ("CH4O", 32026), ("CO2", 43990)])
alt_losses = OrderedDict([("C2H4", 28032), ("CH2O", 30011), ("C2H6", 30046), ("C2H2O", 42011), 
                          ("C3H6", 42047), ("C2H6O", 46042), ("C4H8", 56063), ("C2H4O2", 60021)])


input_file_path = get_config(section='params', key='input_file_path')

multiple = get_config(section='params', key='multiple') #multiple of each neutral loss to consider to find the next peak
tolerance = get_config(section='params', key='tolerance') #tolerance of +-1mDa
nominal_tolerance = get_config(section='params', key='nominal_tolerance') # fragments within +-nominal_tolerance will be considered precursors.
overlap_len = get_config(section='params', key='overlap_len') # Overlap length threshold of pathways for creating families. Should be at least 2.
