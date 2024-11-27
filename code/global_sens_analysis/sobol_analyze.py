from SALib.analyze import sobol

import numpy
import scipy
import pickle

problem = None
with open('Problem.pkl', 'rb') as f:
    problem = pickle.load(f)

Y = None
with open('Y.pkl', 'rb') as f:
    Y = pickle.load(f)

# Perform Sobol analysis on model outputs.
# Returns a dictionary with keys 'S1', 'S1_conf', 'ST', and 'ST_conf',
# where each entry is a list of size D (the number of parameters)
# containing the indices in the same order as the parameter file.
# If calc_second_order is True,
# the dictionary also contains keys 'S2' and 'S2_conf'.

Si = sobol.analyze(problem,
                   Y,
                   calc_second_order=False,
                   num_resamples=100,
                   conf_level=0.95,
                   parallel=False)

with open('Si.pkl', 'wb') as f:
    pickle.dump(Si, f, pickle.HIGHEST_PROTOCOL)
