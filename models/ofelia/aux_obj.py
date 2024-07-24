import os
import numpy as np

from pyXSteam.XSteam import XSteam
from scipy.stats import linregress as LR

from dolfinx.fem import (assemble_scalar, form)
import ufl
from ufl import inner, grad
import openmc

from mpi4py import MPI

# Update material files
class updateXML():
    def __init__(self, mat_dict: dict, n_div: int, pressure : float = 155., Tmin : float = 280, Tmax : float = 330):
        self.mat_dict = mat_dict
        self.n_div = n_div
        
        steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)

        Temperature = np.linspace(Tmin, Tmax, 50) + 273.15 # Kelvin
        density = np.zeros_like(Temperature)

        for ii in range(len(density)):
            density[ii] = steamTable.rho_pt(pressure, Temperature[ii]-273.15) / 1e3 # g/cm3

        self.fitting = LR(Temperature, density)
            
    def update(self, Tf, Tc, sn, index):

        for k in range(self.n_div):
            # print(Tf[k])
            self.mat_dict['fuel'][k].temperature = Tf[k]

        for k in range(self.n_div):
            # print(Tc[k])
            rho = self.fitting.intercept + self.fitting.slope*Tc[k]
            self.mat_dict['coolant'][k].temperature = Tc[k]
            self.mat_dict['coolant'][k].set_density('g/cm3', rho)

        # Collect the materials together and export to XML
        mat_list = list()
        for key in list(self.mat_dict.keys()):
            for mat in self.mat_dict[key]:
                mat_list.append(mat)
        materials = openmc.Materials(mat_list)
        materials.export_to_xml()

        #Create a new folder, copying the it0 and moving the new 'materials.xml' 
        path = os.getcwd()
        folder_0 = path + '/build_xml/it0'
        new_folder = path + '/build_xml/it' + str(index)
        os.system('cp -r ' + folder_0 + ' ' + new_folder)
        os.system('mv materials.xml ' + new_folder + '/materials.xml')
        UpdateParticles(new_folder, sn)


############################################################################################################

# This class is used in FEniCSx to compute norms, integrals and inner products
class norms():
  def __init__(self, funSpace, domain):
    self.domain = domain
    self.trial = ufl.TrialFunction(funSpace)
    self.test  = ufl.TestFunction(funSpace)

    self.dx = ufl.Measure('dx', domain=self.domain)
    self.L2inner = inner(self.trial, self.test) * self.dx
    self.H1inner = inner(grad(self.trial), grad(self.test)) * self.dx

  def L2norm(self, u):
    repl_form = form(ufl.replace(self.L2inner, {self.trial: u, self.test: u}))
    return np.sqrt( self.domain.comm.allreduce(assemble_scalar(repl_form), op=MPI.SUM) )
    
  def H1norm(self, u):
    repl_form = form(ufl.replace(self.H1inner, {self.trial: u, self.test: u}))
    return np.sqrt( self.domain.comm.allreduce(assemble_scalar(repl_form), op=MPI.SUM) )
    
  def L2innerProd(self, u, v):
    repl_form = form(ufl.replace(self.L2inner, {self.trial: u, self.test: v}))
    return self.domain.comm.allreduce(assemble_scalar(repl_form), op=MPI.SUM)

  def Linftynorm(self, u):
    return self.domain.comm.allreduce(np.max(np.abs(u.x.array)), op = MPI.MAX)
   
# Evaluate the normalization of the quantity 'q', knowing the Power 'P' and the pin-lenght 'l' with radius 'r'
class extract_power():
    def __init__(self, n_div: int, power: float, mesh_size: int, pin_length: float, pin_radius: float,
                      J_to_eV: float, tally_dict: dict):
        self.n_div = n_div
        self.power = power
        self.mesh_size = mesh_size
        
        self.pin_length = pin_length
        self.pin_radius = pin_radius
        
        self.J_to_eV = J_to_eV
        self.tally_dict = tally_dict
        
    def eval(self, sp, i, Ef = 200e6):
 
        # INTEGRAL TALLY
        tally_integral = sp.get_tally(name = self.tally_dict['integral'])

        # Integral fission energy
        heating_integral = tally_integral.get_slice(scores=['kappa-fission']) 
        Qp = float(heating_integral.mean) # measure of the fission energy (eV/src) 

        ### Axial scores
        tally_mesh = sp.get_tally(name = self.tally_dict['mesh_z'])

        # Fission reaction rate (z)
        fissionZ = tally_mesh.get_slice(scores=['fission']) # fission comes in (fissions/source)
        RR_Z, uRR_Z, Area = self.normalisation_z(fissionZ, Qp) # fissions normalized to (fissions / cm3 s)

        # Power density (z)
        dz = self.pin_length / self.mesh_size
        z = np.arange(-self.pin_length/2, self.pin_length/2., dz )
    
        q3 = RR_Z*Ef*self.J_to_eV # (fiss/cm3 s * eV/fiss * J/eV) = (W/cm3)
        q3std = uRR_Z*Ef*self.J_to_eV

        return q3, z, q3std, Area

    # very important function :)
    def normalisation_z(self, qty_to_norm, Qp):
        
        H1 = self.J_to_eV * Qp # (J/source)
        Vol = (self.pin_radius**2) * self.pin_length * np.pi # pin volume (cm3)
        
        f = self.mesh_size * self.power / (H1 * Vol) # Normalization factor ( source/(s cm3) )
        #f = self.power / (H1 * Vol) # Normalization factor ( source/(s cm3) )
        
        # Put quantities in the right shape
        q_mean = qty_to_norm.mean # (fissions / source)
        q_mean.shape = self.mesh_size 
        
        q_std = qty_to_norm.std_dev
        q_std.shape = self.mesh_size 

        q_mean_normalised = q_mean * f # (fissions / (s cm3))
        q_std_normalised = q_std * f 
        
        return q_mean_normalised, q_std_normalised, Vol/self.pin_length

    def getSpectrum(self, sp, phiE_list):
        
        ### Energy score
        tally_energy = sp.get_tally(name = self.tally_dict['spectrum'])

        # Flux in energy
        energy_filter = tally_energy.filters[0]
        energies = energy_filter.bins[:, 0]

        # Get the flux values
        mean = tally_energy.get_values(value='mean').ravel()
        uncertainty = tally_energy.get_values(value='std_dev').ravel()

        EnergyStairs = np.zeros(len(energies)+1)
        EnergyStairs[0:-1] = energies
        EnergyStairs[-1] = energies[-1]

        phiE_list[0].append(np.array(mean))
        phiE_list[1].append(np.array(uncertainty))

        return phiE_list, EnergyStairs



def UpdateParticles(new_folder, sn):
    file = open (new_folder+'/settings.xml', "r")
    list_lines = file.readlines()
    list_lines[3] = ( "  <particles>"+str(int(sn))+"</particles>\n" )
    file = open (new_folder+'/settings.xml', "w")
    file.writelines(list_lines)
    file.close()
    
    
def RemoveFolders(path=None):
    if path is None:
        path = os.getcwd()
    folder_to_check = path + '/build_xml'
    n_folders = len(next(os.walk(folder_to_check))[1])
    
    if n_folders>1:
        for k in range(n_folders-1):
            folder_to_delete = folder_to_check + '/it' + str(k+1)
            os.system('rm -r ' + folder_to_delete)