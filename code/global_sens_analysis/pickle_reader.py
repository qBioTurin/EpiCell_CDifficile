
# python3 pickle_reader.py

import pickle
import sys
import pandas

def read_pickle_file(file):
  with open(file, "rb") as f:
    pickle_data = pickle.load(f)
  print(type(pickle_data))
  return pickle_data

# Defining the model inputs.
Problem = read_pickle_file(file = "/home/raucello/EpiCell_CDifficile/code/global_sens_analysis/Problem.pkl")
Problem = pandas.DataFrame.from_dict(Problem)
Problem.to_csv("Problem.csv", sep = ",", index=True, header=True)

# Model outputs
Y = read_pickle_file(file = "/home/raucello/EpiCell_CDifficile/code/global_sens_analysis/Y.pkl")
Y = pandas.DataFrame.from_dict(Y)
Y.to_csv("Y.csv", sep = ",", index=False, header=False)

# Perform Sobol analysis on model outputs.
# Returns a dictionary with keys 'S1', 'S1_conf', 'ST', and 'ST_conf',
# where each entry is a list of size D (the number of parameters)
# containing the indices in the same order as the parameter file.
# If calc_second_order is True,
# the dictionary also contains keys 'S2' and 'S2_conf'.

# dictionary containing the indices
Si = read_pickle_file(file = "/home/raucello/EpiCell_CDifficile/code/global_sens_analysis/Si.pkl")
Si = pandas.DataFrame.from_dict(Si)
# print(Si.head(5))
Si.to_csv("Si.csv", sep = ",", index=True, header=True)
