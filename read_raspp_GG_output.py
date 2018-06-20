import pickle
import matplotlib.pyplot as plt

l = pickle.load(open('20180612_libraries_gg.p','rb'))

M_threshold = 177
energy_threshold = 43

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
#chosen_lib = (lib, l[lib])
#pickle.dump(chosen_lib, open('chosen_lib.p','wb'))