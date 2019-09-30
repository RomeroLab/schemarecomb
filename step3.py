""" Choose a library based on SCHEMA energy and mutational levels.

Run syntax: "python step3.py <libraries_fn> <chosen_lib_fn>"
Example: "python step3.py bgl3_libraries.json bgl3_chosen_lib.json"

Command Line Args:
    libraries_fn: name of libraries file generated in step 2
    chosen_lib_fn: name of output chosen library file

Outputs:
    chosen library in <chosen_lib_fn>, used in step 4
    libraries.png: graph of SCHEMA-RASPP curve with chosen_library


You can also replace all instances of "sys.argv" in the code with the input
filenames directly, then run "python step2.py".

As written, the script automatically chooses the library with maximum m-e where
m is the average number of mutations and e is the SCHEMA energy for the
library. To manually select a library, comment out the automatic selection,
find the breakpoints key for the desired library, then save the desired library
key-value pair in a tuple for json packaging.

Matplotlib is a required package. The Romero lab group server has it installed.
"""

import json
import sys

import matplotlib.pyplot as plt

libraries_fn = sys.argv[1]
chosen_lib_fn = sys.argv[2]

with open(libraries_fn, 'r') as f:
    libraries = {tuple(k): v for k, v in json.load(f)}


gg_filtered_libs = {k: v for k, v in libraries.items() if v['GG_prob'] >= 0.95}

chosen_lib_bps = max(gg_filtered_libs, key=lambda x: gg_filtered_libs[x]['M']
                     - gg_filtered_libs[x]['energy'])
chosen_lib_attrs = libraries[chosen_lib_bps]

graph_data = [(l['M'], l['energy'], l['GG_prob']) for l in libraries.values()]
graph_chosen_lib = (chosen_lib_attrs['M'], chosen_lib_attrs['energy'])

plt.figure(figsize=(16, 9))
m, e, g = zip(*graph_data)
plt.scatter(m, e, c=g, vmax=1, vmin=.95)
plt.annotate('chosen library', graph_chosen_lib)
plt.colorbar().set_label('GG_prob', rotation=270, labelpad=12)
plt.title('Libraries')
plt.xlabel('M')
plt.ylabel('energy')
plt.savefig('libraries.png')

chosen_lib = (chosen_lib_bps, chosen_lib_attrs)
with open(chosen_lib_fn, 'w') as f:
    json.dump(chosen_lib, f)
