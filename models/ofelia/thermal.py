import numpy as np
from scipy.interpolate import interp1d

import dolfinx
from dolfinx import fem
from dolfinx.fem import (Function, FunctionSpace, assemble_scalar, form)
import ufl
from ufl import SpatialCoordinate, inner, grad
from petsc4py import PETSc
from mpi4py import MPI
from pyXSteam.XSteam import XSteam

class thermal_solver():
    def __init__(self, domain: dolfinx.mesh.Mesh, ct: dolfinx.cpp.mesh.MeshTags_int32, ft: dolfinx.cpp.mesh.MeshTags_int32,
                 physical_param : dict, regions_markers: list, robin_mark = int, degree : int = 1):

        # Storing domain, cell tags, face tags and the functional space
        self.domain = domain
        self.ct = ct
        self.ft = ft
        self.funSpace = FunctionSpace(domain, ("Lagrange", degree))
        self.Qn = FunctionSpace(domain, ("DG", 0))

        self.gdim = self.domain.geometry.dim
        self.fdim = self.gdim - 1

        # Defining physical parameters and labels for the regions and boundaries        
        self.regions = regions_markers # the first element (in position [0] must be the fuel)
        self.phys_param = physical_param
        self.robin_mark = robin_mark
        
        self.k = Function(self.Qn) # thermal conductivity
        for idx, regionI in enumerate(self.regions):
            cells = self.ct.find(regionI)
            self.k.x.array[cells] = self.phys_param['th_cond'][idx]
        
        self.h = fem.Constant(self.domain, PETSc.ScalarType(self.phys_param['htc'])) # heat transfer coefficient
        
        # Definition of the trial and test space
        self.T = ufl.TrialFunction(self.funSpace)
        self.v = ufl.TestFunction(self.funSpace)
        self.solution = Function(self.funSpace)
        
        # Definition of the bulk temperature of the coolant and the power density
        self.T_b = Function(self.funSpace)
        self.q  = Function(self.funSpace)
        
        # Definition of the surface and volume element for cylindrical coordinate system
        self.ds = ufl.Measure('ds', domain=self.domain, subdomain_data=self.ft)
        self.dx = ufl.Measure('dx', domain=self.domain, subdomain_data=self.ct)
        
        [self.z_, self.r_] = SpatialCoordinate(domain)

    def assemble(self, direct : bool = False):

        # Assembling the linear forms
        self.left_side  = (inner(self.k * grad(self.T), grad(self.v)) * np.abs(self.r_) * self.dx + 
                           inner(self.h * self.T, self.v) * self.ds(self.robin_mark) )
        
        self.right_side = (inner(self.q, self.v) * np.abs(self.r_) * self.dx(self.regions[0]) +
                           inner(self.h * self.T_b, self.v) * self.ds(self.robin_mark) )
        
        self.a = form(self.left_side)
        self.L = form(self.right_side)
        
        # Creating and storing the matrices
        self.A = fem.petsc.assemble_matrix(self.a)
        self.A.assemble()
        self.b = fem.petsc.create_vector(self.L)
        
        # Definition of the solver
        self.solver = PETSc.KSP().create(self.domain.comm)
        self.solver.setOperators(self.A)
        
        if direct:
            self.solver.setType(PETSc.KSP.Type.PREONLY)
            self.solver.getPC().setType(PETSc.PC.Type.LU)
        else:
            self.solver.setType(PETSc.KSP.Type.CG)
            self.solver.getPC().setType(PETSc.PC.Type.SOR)
        
    def solve(self, power_fun, T_b_fun): # power_fun and T_b_fun should be a callable interpolant
        
        # Updating the power density and the bulk temperature of the coolant
        self.q.interpolate(lambda x: power_fun(x[0], x[1]))
        self.T_b.interpolate(lambda x: T_b_fun(x[0], x[1]))
        
        # Updating the rhs vector
        with self.b.localForm() as loc:
            loc.set(0)
        fem.petsc.assemble_vector(self.b, self.L)
        self.b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        
        # Solving the linear system
        self.solver.solve(self.b, self.solution.vector)
        self.solution.x.scatter_forward()
        
        return self.solution.copy()

    def computeSolidAverageT(self, is_region: int, T_sol: dolfinx.fem.Function, slices: np.ndarray):

        # Definition of the function and the form for the integration
        sliceID = Function(self.Qn)
        integral_form = sliceID * T_sol * self.dx(is_region)
        domain_form   = sliceID * self.dx(is_region)
        
        # Computing the average temperature
        aveT = np.zeros([len(slices)-1])
        for ii in range(len(slices)-1):
            bounds = np.array([slices[ii], slices[ii+1]])
            
            sliceID.interpolate(lambda x: np.piecewise(x[0], [x[0]<bounds[0],
                                                            np.logical_and(x[0]>=bounds[0], x[0]<=bounds[1]),
                                                            x[0]>=bounds[1]],
                                                            [0., 1., 0.]))
                                                            
            aveT[ii] = self.domain.comm.allreduce(assemble_scalar(form(integral_form)), op=MPI.SUM) / self.domain.comm.allreduce(assemble_scalar(form(domain_form)), op=MPI.SUM)
            
        return aveT
    
    def extract_2D_data(self, T, L, R, Nx = 400, Ny = 100):
        
        x_grid = np.linspace(-L/2, L/2, Nx)
        y_grid = np.linspace(-R, R, Ny)

        T_matrix = np.zeros((Nx, Ny))

        for ii in range(Nx):
            points = np.zeros((3, Ny))
            points[0, :] = x_grid[ii]
            points[1, :] = y_grid

            bb_tree = dolfinx.geometry.BoundingBoxTree(self.domain, self.domain.topology.dim)
            cells = []
            points_on_proc = []
            cell_candidates = dolfinx.geometry.compute_collisions(bb_tree, points.T)
            colliding_cells = dolfinx.geometry.compute_colliding_cells(self.domain, cell_candidates, points.T)
            for i, point in enumerate(points.T):
                if len(colliding_cells.links(i))>0:
                    points_on_proc.append(point)
                    cells.append(colliding_cells.links(i)[0])
            xPlot = np.array(points_on_proc, dtype=np.float64)

            T_matrix[ii, :] = T.eval(xPlot, cells).flatten()
        
        X, Y = np.meshgrid(x_grid, y_grid)
        res2d = dict()
        res2d['xgrid'] = x_grid
        res2d['ygrid'] = y_grid
        res2d['X'] = X
        res2d['Y'] = Y
        res2d['T'] = T_matrix
        return res2d

class thermal_inputs():
    def __init__(self, coolant_T: float, coolant_p: float):
        
        # Storing coolant properties
        self.coolant_T = coolant_T
        self.coolant_p = coolant_p
        
        # Computing flow thermophysical properties
        steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)
        self.rho = steamTable.rho_pt(self.coolant_p, self.coolant_T - 273.15) / 1e3 # g/cm3 
        self.cp  = steamTable.Cp_pt( self.coolant_p, self.coolant_T - 273.15) # J/g K
        self.k   = steamTable.tc_pt( self.coolant_p, self.coolant_T - 273.15) / 100 # W/cm-K
        self.mu  = steamTable.my_pt( self.coolant_p, self.coolant_T - 273.15) * 10  # g / cm s^2
        self.Pr  = self.cp * self.mu / self.k

    def compute_htc(self, pitch: float, clad_or: float, u_in: float):
        
        
        # Computing hydraulic diameter
        flow_area = pitch**2 - np.pi * clad_or**2
        Dh = 4 * flow_area / (2 * np.pi * clad_or)

        # Computing equivalent Reynolds Number
        Re = self.rho * u_in * Dh / self.mu
        
        # The convective HTC is computed using the Dittus-Boelter correlation
        h = self.k / Dh * (0.023 * Re**0.8 * self.Pr**0.4) # W/cm2 - K
        self.u_in = u_in
        self.flow_area = flow_area
        
        return h

    def mapping_q3_Tb(self, z_omc: np.ndarray, q3_omc: np.ndarray, Tin: float, L: float, L_active: float, fuel_or: float):
        
        # Mapping q3 from OpenMC to scipy.interpolate
        q3_fuel = interp1d(z_omc, q3_omc, kind='linear', fill_value="extrapolate")
        dimz_full = z_omc.shape[0]
        if z_omc.shape[0] < 100:
            dimz_full = 100
        else:
            dimz_full = z_omc.shape[0]+100
        
        z_full = np.linspace(-L/2, L/2, dimz_full)
        q3_full = np.zeros_like(z_full)
        
        for kk in range(len(z_full)):
            if z_full[kk] < -L_active/2:
                value1 = 0.
            elif (( z_full[kk] >= -L_active/2 ) & ( z_full[kk] <= L_active/2 )):
                value1 = q3_fuel(z_full[kk])
            else:
                value1 = 0.
            q3_full[kk] = value1
               
        q3_fuel = interp1d(z_full, q3_full, kind='linear', fill_value="extrapolate")
        
        # Energy balance to compute the bulk temperature of the fuel
        m_cp = (self.rho * self.flow_area * self.u_in) * self.cp # g/s * J/g/K # cp at 300 C and 155 bar

        T_b = np.zeros((len(z_full), ))
        T_b[0] = Tin
        for ii in range(1, len(z_full)):
            T_b[ii] = T_b[0] + 1. / m_cp * np.pi * fuel_or**2 * np.trapz(q3_fuel(z_full[:ii+1]), z_full[:ii+1])
        
        # Creating Tb interpolant
        self.T_b_fun = interp1d(z_full, T_b, kind='linear',fill_value='extrapolate')
        
        # Storing the results
        mapping_res = dict()
        mapping_res['z'] = z_full
        mapping_res['q3'] = q3_fuel
        mapping_res['T_bulk'] = self.T_b_fun

        return lambda x,y: q3_fuel(x) + 0.0 * y, lambda x,y: self.T_b_fun(x) + 0.0 * y, mapping_res
    
    def computeWaterAverageT(self, slices: np.ndarray, dim_z: int):
        
        # Definition of the average and the z for the integration
        aveT = np.zeros([len(slices)-1])
        z = np.linspace(min(slices), max(slices), dim_z)

        # Computing average temperature
        for ii in range(len(slices)-1):
            bounds = np.array([slices[ii], slices[ii+1]])
            sliceID = lambda x: np.piecewise(x, [x<bounds[0],
                                            np.logical_and(x>=bounds[0], x<=bounds[1]),
                                            x>=bounds[1]],
                                            [0., 1., 0.])
                            
            aveT[ii] = np.trapz( self.T_b_fun(z) * sliceID(z), z) / np.trapz( sliceID(z), z) 
        return aveT
    