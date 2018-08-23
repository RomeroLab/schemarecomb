import raspp_tools
import pickle

# define the library properties min block length, max block length, number of blocks, and parents
minBL = 40
maxBL = 400
num_bl = 8

alignment_fn = 'alignment.p'
contacts_fn = 'contacts.p'
libraries_fn = 'libraries.p'

alignment = pickle.load(open(alignment_fn, 'rb'))
contacts=pickle.load(open(contacts_fn, 'rb'))

# A breakpoint specifies the first position of a new block. 
breakpoints = raspp_tools.find_GG_breakpoints(alignment)

# The same E, but weighted by the values in the contacts dict. This could be the contact frequency. To get unweighted, just set all contact weight==1""
E_matrix = raspp_tools.generate_weighted_E_matrix(alignment,contacts)

# generate all allowed blocks
blocks = raspp_tools.generate_blocks(breakpoints,minBL,maxBL)

# run RASPP
#libraries = raspp_tools.fast_shortest_path_recombination(num_bl,blocks,E_matrix,False) # fast version
libraries = raspp_tools.shortest_path_recombination(num_bl,blocks,E_matrix,False) # slow, but thorough version

print('\n Updating M', flush=True)

# Add M values to the libraries dictionary
raspp_tools.update_M(libraries,alignment)

pickle.dump(libraries, open(libraries_fn, 'wb'))

#Requires 'normalized_ligation_counts_18h_37C.p' file to be in same dir
libraries = raspp_tools.update_GG_prob(libraries, alignment)

pickle.dump(libraries, open('gg_' + libraries_fn, 'wb'))
