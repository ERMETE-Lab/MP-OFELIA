import numpy as np
import os

path_to_run = './build_xml/it0/'
if not os.path.exists(path_to_run):
    os.makedirs(path_to_run)

### MC settings
s1_val = 750 # initial particle population
batches = 500 # number of cycles
n_div = 8 # Number of axial divisions
Power = 65e3 # W
initialUniformSource = False
shannonEntropy = False
meshSize = 100

tallyDict = {'integral': "Integral Quantities",
             'mesh_xy': "Mesh xy tally",
             'mesh_z': "Mesh z tally",
             'spectrum': "Flux Spectrum"}

# Picard settings
maxIter = 25
tol   = 0.05

# Geometrical and physical parameters
l_active=366 #  Active length (cm)
L_clad_o = 1.27
total_length = l_active  + 2 * L_clad_o

pin_r = 0.819/2 # Pin radius (cm)
Tin = 273.15 + 290 # K
Tout = Tin + 40
fuel_or = 0.819/2
clad_ir = 0.836/2
clad_or = 0.95/2
pitch = 1.26 # cm
square_side = 1.26 # [pin pitch]
plug_length = 1.27
end_domain_top = 200
end_domain_bot = -200
J_to_eV = 1.602E-19 # J/eV
pin_length = l_active + 2 * plug_length # total length of the pin

## Parameters to compute the bulk T of the coolant and HTC
pressure = 155 # bar
T_w_average = Tin + 20 # K # this value actually depends on the power imposed !!!!
u_in = 5.3 * 100 #cm/s