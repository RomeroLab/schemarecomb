""" Choose a library based on SCHEMA energy and mutational levels.

Run syntax: "python step3.py <libraries_fn> <chosen_lib_fn>"
Example: "python step3.py bgl3_libraries.json bgl3_chosen_lib.json"

Command Line Args:
    libraries_fn: name of libraries file generated in step 2
    chosen_lib_fn: name of output chosen library file

Outputs:
    chosen library in <chosen_lib_fn>, used in step 4

You can also replace all instances of "sys.argv" in the code with the input
filenames directly, then run "python step2.py".

Matplotlib is a required package. The Romero lab group server has it installed.
"""

import json
import sys

import matplotlib.pyplot as plt

libraries_fn = sys.argv[1]
chosen_lib_fn = sys.argv[2]

with open(libraries_fn, 'r') as f:
    l = json.load(f)

M_threshold = 1
energy_threshold = 1000

data_points = [(l[lib]['energy'],l[lib]['M'],l[lib]['GG_prob']) 
            for lib in l.keys() if l[lib]['energy'] < energy_threshold and l[lib]['M'] > M_threshold]

plt.figure()
e,m,g = zip(*data_points)
plt.scatter(m,e,c=g,vmax=1,vmin=.95)
plt.colorbar().set_label('GG_prob',rotation=270,labelpad=12)
plt.title('Libraries')
plt.xlabel('M')
plt.ylabel('energy')
plt.show()

accepted_libraries = {k:v for (k,v) in l.items() if v['energy'] < energy_threshold and v['M'] > M_threshold}

print(accepted_libraries)

##set lib = l key of chosen library 
chosen_lib = (lib, l[lib])
with open(chosen_lib_fn, 'w') as f:
    json.dump(chosen_lib, f)
