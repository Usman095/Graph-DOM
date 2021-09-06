from collections import OrderedDict

elements = {"C":12000, "H":1008, "O":15995, "N":14003, "S":31972}
atom_list = ['C', 'H', 'O']

#the order of neutral losses will be the same in the final output
neutral_losses = OrderedDict([("CH2", 14016), ("CH3", 15023), ("O", 15995), ("CH4", 16031),
                              ("H2O", 18010), ("CO", 27995), ("CH4O", 32026), ("CO2", 43990)])
alt_losses = OrderedDict([("C2H4", 28032), ("CH2O", 30011), ("C2H6", 30046), ("C2H2O", 42011), 
                          ("C3H6", 42047), ("C2H6O", 46042), ("C4H8", 56063), ("C2H4O2", 60021)])


input_file_path = '2D_MSMS_CHO_02_08_2021.xlsx'

multiple = 4 #multiple of each neutral loss to consider to find the next peak
tolerance = 3 #tolerance of +-1mDa
nominal_tolerance = 1 # fragments within +-nominal_tolerance will be considered precursors.
overlap_len = 3 # Overlap length threshold of pathways for creating families. Should be at least 2.