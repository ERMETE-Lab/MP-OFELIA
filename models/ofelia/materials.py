import openmc

from .constants import *

################################ OpenMC materials ######################

# Dictionary which contains the used materials
mat_dict = dict()
mat_dict['coolant'] = [] # list of water layers
mat_dict['fuel'] = [] # list of fuel slices
mat_dict['non_updated'] = [] # list of materials whose temperature is constant (i.e., He, Clad)

# water0 and fuel0 are used for cloning the real materials that will be put in the system
idx_mat = 0
mat_dict['non_updated'].append(openmc.Material(name='water0'))
mat_dict['non_updated'][idx_mat].set_density('g/cm3', 0.75)
mat_dict['non_updated'][idx_mat].add_element('H', 2)
mat_dict['non_updated'][idx_mat].add_element('O',1)
mat_dict['non_updated'][idx_mat].add_s_alpha_beta('c_H_in_H2O')

# Cloning water materials
for kk in range(n_div):
    mat_dict['coolant'].append(mat_dict['non_updated'][idx_mat].clone())
    mat_dict['coolant'][kk].name = 'coolant_'+str(kk+1)    

idx_mat += 1
mat_dict['non_updated'].append(openmc.Material(name='fuel0'))
mat_dict['non_updated'][idx_mat].set_density('g/cm3', 10.45)
mat_dict['non_updated'][idx_mat].add_nuclide('U235', 9.3472e-4, 'ao')
mat_dict['non_updated'][idx_mat].add_nuclide('U238', 2.1523e-2, 'ao')
mat_dict['non_updated'][idx_mat].add_nuclide('U234', 9.1361e-6, 'ao')
mat_dict['non_updated'][idx_mat].add_nuclide('O16', 4.4935e-02, 'ao')
mat_dict['non_updated'][idx_mat].temperature = 800

# Cloning fuel materials
for kk in range(n_div):
    mat_dict['fuel'].append(mat_dict['non_updated'][idx_mat].clone()) 
    mat_dict['fuel'][kk].name = 'fuel_'+str(kk+1)
     
idx_mat += 1
mat_dict['non_updated'].append(openmc.Material(name='Helium for gap'))
mat_dict['non_updated'][idx_mat].set_density('g/cm3', 0.001598)
mat_dict['non_updated'][idx_mat].add_element('He', 2.4044e-4)
mat_dict['non_updated'][idx_mat].temperature = 600

idx_mat += 1
mat_dict['non_updated'].append(openmc.Material(name="Zirc4"))
mat_dict['non_updated'][idx_mat].set_density('g/cm3', 6.44)
mat_dict['non_updated'][idx_mat].add_nuclide('O16', 1.192551825E-03, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('O17', 4.82878E-07, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Cr50', 4.16117E-05, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Cr52', 8.34483E-04, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Cr53', 9.64457E-05, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Cr54', 2.446E-05, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Fe54', 1.1257E-04, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Fe56', 1.8325E-03, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Fe57', 4.3077E-05, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Fe58', 5.833E-06, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Zr90', 4.9786E-01, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Zr91', 1.0978E-01, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Zr92', 1.6964E-01, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Zr94', 1.7566E-01, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Zr96', 4.28903E-02, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Sn116', 1.98105E-03, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Sn117', 1.05543E-03, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Sn119', 1.20069E-03, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Sn120', 4.5922E-03, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Sn122', 6.63497E-04, 'wo')
mat_dict['non_updated'][idx_mat].add_nuclide('Sn124', 8.43355E-04, 'wo')
mat_dict['non_updated'][idx_mat].temperature = 600 # [K]

idx_mat += 1
mat_dict['non_updated'].append(openmc.Material(name="Water inlet"))
mat_dict['non_updated'][idx_mat].set_density('g/cm3', 0.996557)
mat_dict['non_updated'][idx_mat].add_element('H', 2)
mat_dict['non_updated'][idx_mat].add_element('O', 1)
mat_dict['non_updated'][idx_mat].add_s_alpha_beta('c_H_in_H2O')
mat_dict['non_updated'][idx_mat].temperature = Tin

idx_mat += 1
mat_dict['non_updated'].append(openmc.Material(name="Water outlet"))
mat_dict['non_updated'][idx_mat].set_density('g/cm3', 0.996557)
mat_dict['non_updated'][idx_mat].add_element('H', 2)
mat_dict['non_updated'][idx_mat].add_element('O', 1)
mat_dict['non_updated'][idx_mat].add_s_alpha_beta('c_H_in_H2O')
mat_dict['non_updated'][idx_mat].temperature = Tout