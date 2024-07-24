from math import log10, pi

import matplotlib.pyplot as plt
import numpy as np
import openmc

###############################################################################
# --- MATERIALS 
from .materials import *

mat_list = list()
for key in list(mat_dict.keys()):
    for mat in mat_dict[key]:
        mat_list.append(mat)
materials = openmc.Materials(mat_list)
materials.export_to_xml(path_to_run+'materials.xml')

###############################################################################
# --- GEOMETRY

# Create cylindrical surfaces
Fuel_or = openmc.ZCylinder(r=fuel_or, name='Fuel OR')
Clad_ir = openmc.ZCylinder(r=clad_ir, name='Clad IR')
Clad_or = openmc.ZCylinder(r=clad_or, name='Clad OR')

# Create planes for the channel
square_side = pitch
east_boundary = openmc.XPlane(x0= square_side/2, boundary_type = 'reflective', name='right boundary')
west_boundary = openmc.XPlane(x0=-square_side/2, boundary_type = 'reflective', name='left  boundary')
north_boundary = openmc.YPlane(y0=square_side/2, boundary_type = 'reflective', name='north boundary')
south_boundary = openmc.YPlane(y0=-square_side/2, boundary_type = 'reflective', name='south boundary')

z0_bottom = -l_active/2
z0_top = l_active/2
dz = (z0_top-z0_bottom)/n_div

# Create boundaries for the 3D pin
top_domain = openmc.ZPlane(z0=end_domain_top, boundary_type = 'vacuum', name ='top domain')
top_clad = openmc.ZPlane(z0=z0_top+plug_length,  name ='top clad')
top_active= openmc.ZPlane(z0=z0_top, name='top boundary')

bot_active= openmc.ZPlane(z0=z0_bottom,name='bot boundary')
bot_clad = openmc.ZPlane(z0=z0_bottom-plug_length, name ='bot clad')
bot_domain = openmc.ZPlane(z0=end_domain_bot, boundary_type = 'vacuum', name ='bot domain')


# Create planes for material subdivision
z_planes = []

for j in range(n_div-1):#7 piani in mezzo
    plane = openmc.ZPlane(z0=z0_bottom + dz*(j+1))
    z_planes.append(plane)


# Create cells and assigning materials to regions in ascendent order
cells_list = []
for div in range(n_div):
    if div==0: #bottom zone
        fuel_cell=openmc.Cell(fill=mat_dict['fuel'][div], region = -Fuel_or & +bot_active & -z_planes[div], name='fuelzone' + str(div))
        water_cell=openmc.Cell(fill=mat_dict['coolant'][div], region=+Clad_or & +west_boundary & -east_boundary  & +south_boundary & -north_boundary  & +bot_active & -z_planes[div], name='waterzone' + str(div))
        
    elif div==(n_div-1):
        fuel_cell=openmc.Cell(fill=mat_dict['fuel'][div], region = -Fuel_or & +z_planes[div-1] & -top_active, name='fuelzone' + str(div))
        water_cell=openmc.Cell(fill=mat_dict['coolant'][div], region=+Clad_or & +west_boundary & -east_boundary  & +south_boundary & -north_boundary  & +z_planes[div-1] & -top_active, name='waterzone' + str(div))

    else:
        fuel_cell=openmc.Cell(fill=mat_dict['fuel'][div], region = -Fuel_or & +z_planes[div-1] & -z_planes[div], name='fuelzone' + str(div))
        water_cell=openmc.Cell(fill=mat_dict['coolant'][div], region=+Clad_or & +west_boundary & -east_boundary  & +south_boundary & -north_boundary  & +z_planes[div-1] & -z_planes[div], name='waterzone' + str(div))

    cells_list.append(fuel_cell)
    cells_list.append(water_cell)



# Adding Clad and Helium which are not updated in temperature
gap = openmc.Cell(fill=mat_dict['non_updated'][2], region=+Fuel_or & -Clad_ir & +bot_active & -top_active)
clad = openmc.Cell(fill=mat_dict['non_updated'][3], region=+Clad_ir & -Clad_or & +bot_active & -top_active)
clad_top = openmc.Cell(fill=mat_dict['non_updated'][3], region= -Clad_or & +top_active & -top_clad)
clad_bot = openmc.Cell(fill=mat_dict['non_updated'][3], region= -Clad_or & +bot_clad & -bot_active)
fill_water_clad_bot = openmc.Cell(fill=mat_dict['non_updated'][4], region=+Clad_or & +west_boundary & -east_boundary  & +south_boundary & -north_boundary  & +bot_clad & -bot_active)
fill_water_clad_top = openmc.Cell(fill=mat_dict['non_updated'][5], region=+Clad_or & +west_boundary & -east_boundary  & +south_boundary & -north_boundary  & +top_active & -top_clad)
fill_water_bot =  openmc.Cell(fill=mat_dict['non_updated'][4], region= +west_boundary & -east_boundary  & +south_boundary & -north_boundary  & +bot_domain & -bot_clad)
fill_water_top =  openmc.Cell(fill=mat_dict['non_updated'][5], region= +west_boundary & -east_boundary  & +south_boundary & -north_boundary  & +top_clad & -top_domain)

cells_list.append(gap)
cells_list.append(clad)
cells_list.append(clad_top)
cells_list.append(clad_bot)
cells_list.append(fill_water_clad_top)
cells_list.append(fill_water_clad_bot)
cells_list.append(fill_water_top)
cells_list.append(fill_water_bot)

# Create a geometry and export to XML
geometry = openmc.Geometry(cells_list)
geometry.export_to_xml(path_to_run+'geometry.xml')

###############################################################################
# --- SETTINGS

# Indicate how many particles to run
settings = openmc.Settings()
settings.batches = batches
settings.inactive = round(batches / 5)
settings.particles = s1_val

if initialUniformSource:
    #Create an initial uniform spatial source distribution over fissionable zones
    lower_left = (-square_side/2, -square_side/2, -35)
    upper_right = (square_side/2, square_side/2, 35)
    uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
    settings.source = openmc.source.Source(space=uniform_dist)

if shannonEntropy:
    # For source convergence checks, add a mesh that can be used to calculate the Shannon entropy
    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.lower_left = (-Fuel_or.r, -Fuel_or.r)
    entropy_mesh.upper_right = (Fuel_or.r, Fuel_or.r)
    entropy_mesh.dimension = (10, 10)
    settings.entropy_mesh = entropy_mesh

settings.temperature = {'method': 'interpolation'}
settings.export_to_xml(path_to_run+'settings.xml')

################################################################################
# --- TALLIES

#### 1 #### SPATIAL TALLY

# Regular mesh on z direction 
mesh_z = openmc.RegularMesh()
mesh_z.dimension = [1, 1, meshSize] #z mesh
mesh_z.lower_left = [-square_side/2, -square_side/2, -(l_active + 2*plug_length)/2]
mesh_z.upper_right = [square_side/2, square_side/2, +(l_active + 2*plug_length)/2]

# Create a mesh filter that can be used in a tally
mesh_z_filter = openmc.MeshFilter(mesh_z)

# Now use the mesh filter in a tally and indicate what scores are desired
mesh_z_tally = openmc.Tally(name=tallyDict['mesh_z'])
mesh_z_tally.filters = [mesh_z_filter]
mesh_z_tally.scores = ['flux', 'fission']

#### 2 #### INTEGRAL TALLY 
integral_tally = openmc.Tally(name = tallyDict['integral'])
integral_tally.scores = ['kappa-fission']

#### 3 #### ENERGY SPECTRUM
# Let's also create a tally to get the flux energy spectrum. We start by creating an energy filter
e_min, e_max = 1e-5, 20.0e6
groups = 100
energies = np.logspace(log10(e_min), log10(e_max), groups + 1)
energy_filter = openmc.EnergyFilter(energies)

spectrum_tally = openmc.Tally(name=tallyDict['spectrum'])
spectrum_tally.filters = [energy_filter]
spectrum_tally.scores = ['flux']

# Instantiate a Tallies collection and export to XML
tallies = openmc.Tallies([mesh_z_tally, integral_tally, spectrum_tally])
tallies.export_to_xml(path_to_run+'tallies.xml')

###############################################################################
# --- PLOT

p_xz = openmc.Plot()
p_xz.basis = 'xz'
#p_xz.origin(0.0 , 2.0, 0.0) #plot xz centrato in y=0
p_xz.filename = './pictures/pin3D_xz'
p_xz.width = (2, 400)
p_xz.pixels = (2000, 2000)
p_xz.color_by = 'material'


p_xy = openmc.Plot()
p_xy.basis = 'xy'
#p_xy.origin(0.0 , 0.0, 2.0) #plot xy centrato in z=0
p_xy.filename = './pictures/pin3D_xy'
#p_xy.width = (-3, 3)
p_xy.pixels = (2000, 2000)
p_xy.color_by = 'material'

plots = openmc.Plots([p_xy, p_xz])
plots.export_to_xml(path_to_run+'plots.xml')

# openmc.plot_geometry()
