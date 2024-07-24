import dolfinx
import numpy as np
import ufl
from tqdm import tqdm
from dolfinx import fem
from dolfinx.fem import (Function, FunctionSpace, assemble_scalar, form, 
                        locate_dofs_topological, dirichletbc, locate_dofs_geometrical)
from ufl import grad, inner, div, nabla_grad, dot
from petsc4py import PETSc

class steady_thermal_diffusion():
    r"""
    This class implements the solution of the steady thermal diffusion equation with internal power generation made by neutrons.
    
    Parameters
    ----------
    domain : dolfinx.mesh.Mesh
        Mesh imported.
    ct : dolfinx.cpp.mesh.MeshTags_int32
        Cell tags of the regions in the domain.
    ft : dolfinx.cpp.mesh.MeshTags_int32
        Face tags of the boundaries.
    physical_param : dict
        Dictionary with the main physical parameters at reference condition.
    regions_markers : dict
        List of the region markers inside the domain.
    void_bound_marker : int
        Integer describing the index of the void boundary onto which BC is imposed (Dirichlet or Robin).
    h : float, optional (Default = None)
        If not `None`, it, denoted as :math:`h`, represents the heat transfer coefficients in Newton's cooling law.
    TD : float, optional (Default = None)
        If not `None`, it, denoted as :math:`T_D`, represents either the bulk temperature in Newton's cooling law or the Dirichlet condition to impose.
    coupling : str, optional (Default = None)
        If not `None`, it indicates the type of coupling with the temperature: *log* and *sqrt* are available.
            
    """
    def __init__(self, domain: dolfinx.mesh.Mesh, ct: dolfinx.cpp.mesh.MeshTags_int32, ft: dolfinx.cpp.mesh.MeshTags_int32, 
                       physical_param: dict, regions_markers: list, void_bound_marker: int, 
                       h : float = None, TD: float = 300., coupling=None):

        self.domain = domain
        self.ct = ct 
        self.ft = ft
        self.fdim = domain.geometry.dim - 1
        self.void_bound_marker = void_bound_marker
        
        metadata = {"quadrature_degree": 4}
        self.dx = ufl.Measure('dx', domain=self.domain, subdomain_data=self.ct, metadata=metadata)
        self.ds = ufl.Measure('ds', domain=self.domain, subdomain_data=self.ft, metadata=metadata)

        # Defining the number of energy groups
        self.G = physical_param['Energy Groups']
        
        # Temperature Functional spaces
        self.finiteElement = ufl.FiniteElement("Lagrange", domain.ufl_cell(), 1)
        self.V = FunctionSpace(self.domain, self.finiteElement)
        
        # Material Properties Functional spaces
        self.Qn = FunctionSpace(domain, ("DG", 0))
        if coupling is not None:
            self.Told = Function(self.Qn)
            self.Tref = physical_param['Tref']

        self.coupling = coupling
        self.phys_param = physical_param
        self.regions = regions_markers 

        # Parameters
        self.k = Function(self.Qn) # thermal conductivity
        for idx, regionI in enumerate(self.regions):
            cells = self.ct.find(regionI)
            self.k.x.array[cells] = self.phys_param['th_cond'][idx]
        self.xs_fg = list()
        for g in range(self.G):
            self.xs_fg.append(Function(self.Qn))

        self.Ef = physical_param['Ef']
        self.k_eff = Function(self.Qn)
    
        # Definig trial and test functions
        self.T = ufl.TrialFunction(self.V)
        self.theta = ufl.TestFunction(self.V)
        
        # Setting up boundary conditions
        self.TD = Function(self.V)
        self.TD.x.set(TD)
        self.h = h # HTC - if None, Dirichlet is imposed
        if h is None:
            self.bcs = [dirichletbc(self.TD, locate_dofs_topological(self.V, self.fdim, self.ft.find(self.void_bound_marker)))]
        
        # Create solution function
        self.solution = Function(self.V)
        
        # Creating flux function and power form
        self.phi_g = list()
        for g in range(self.G):
            self.phi_g.append(Function(self.V))
        
        self.q3 = self.Ef * sum([self.xs_fg[g] * self.phi_g[g] for g in range(self.G)]) / self.k_eff
        
    def assembleForm(self, direct = True):
        r"""
        This function assembles the variational forms and the correspondent linear algebra structures.
        
        Parameters
        ----------
        direct : bool, optional (Default=True)
            If `True`, a direct linear solver is used (Gauss Elimination), otherwise an iterative one is adopted (preconditioned GMRES with ILU).
            
        """
        self.left_side   = self.k * dot(grad(self.T), grad(self.theta)) * self.dx
        self.right_side  = dot(self.q3, self.theta) * self.dx

        # Boundary condition (convection)
        if self.h is not None:
            self.left_side  += dot(self.h * self.T, self.theta) * self.ds(self.void_bound_marker)
            self.right_side += dot(self.h * self.TD, self.theta) * self.ds(self.void_bound_marker)
        
        self.bilinear  = form(self.left_side)
        self.linear  = form(self.right_side)
        
        self.A = fem.petsc.create_matrix(self.bilinear)
        self.A.zeroEntries()
        
        if self.h is None:
            fem.petsc.assemble_matrix(self.A, self.bilinear, self.bcs)
        else:
            fem.petsc.assemble_matrix(self.A, self.bilinear)
    
        self.A.assemble()  
        self.b = fem.petsc.create_vector(self.linear)

        self.solver = PETSc.KSP().create(self.domain.comm)
        self.solver.setOperators(self.A)

        if direct:
            self.solver.setType(PETSc.KSP.Type.PREONLY)
            self.solver.getPC().setType(PETSc.PC.Type.LU)  
        else: 
            self.solver.setType(PETSc.KSP.Type.GMRES)
            self.solver.getPC().setType(PETSc.PC.Type.ILU)   
            
    def update_xs(self):
        """
        This function updates the cross sections according to the coupling scheme chosen.
        """
        
        if self.coupling is None: # no thermal feedback
            for idx, regionI in enumerate(self.regions):
                cells = self.ct.find(regionI)
                for g in range(self.G):
                    self.xs_fg[g].x.array[cells]    = self.phys_param['xs_f'][g][idx]
                    
        elif self.coupling == 'log': # log law
            for idx, regionI in enumerate(self.regions):
                cells = self.ct.find(regionI)
                for g in range(self.G):
                    self.xs_fg[g].x.array[cells]    = self.phys_param['xs_f'][g][0][idx]    + self.phys_param['xs_f'][g][1][idx]    * np.log(self.Told.x.array[cells] / self.Tref)
                     
        elif self.coupling == 'sqrt': # sqrt law
            for idx, regionI in enumerate(self.regions):
                cells = self.ct.find(regionI)
                for g in range(self.G):
                    self.xs_fg[g].x.array[cells]    = self.phys_param['xs_f'][g][0][idx]    * ( 1 + self.phys_param['xs_f'][g][1][idx]    * ((np.sqrt(self.Told.x.array[cells]) - np.sqrt(self.Tref))) )
                    
    def solve(self, phi: list[dolfinx.fem.Function], k_eff: float, temperature: dolfinx.fem.Function = None):
        r"""
        This function solves the steady equation.
        
        Parameters
        ----------
        phi : list[dolfinx.fem.Function]
            List of `dolfinx.fem.Function` for the internal power generation, :math:`q''' = E_f\cdot \sum_{g=1}^G \Sigma_{f,g}\phi_g`.
        k_eff : float
            Effective multiplication factor.
        temperature : dolfinx.fem.Function, optional (Default = None)
            Temperature field as a `dolfinx.fem.Function` used to calculate the cross sections.
        
        Returns
        -------
        solution : dolfinx.fem.Function
            Output of the linear system.
        """
        # Update k_eff
        self.k_eff.x.set(k_eff)
        
        # Updating fluxes
        for g in range(self.G):
            if len(self.phi_g[g].x.array[:]) == len(phi[g].x.array[:]):
                self.phi_g[g].x.array[:] = phi[g].x.array[:]
            else:
                self.phi_g[g].interpolate(fem.Expression( phi[g], self.V.element.interpolation_points() ))

        # Update XS
        if self.coupling is not None:
            assert(temperature is not None)
            # Updating temperature
            if len(self.Told.x.array[:]) == len(temperature.x.array[:]):
                self.Told.x.array[:] = temperature.x.array[:]
            else:
                self.Told.interpolate(fem.Expression( temperature, self.Qn.element.interpolation_points() ))
                
        self.update_xs()

        with self.b.localForm() as loc_b:
            loc_b.set(0)
        fem.petsc.assemble_vector(self.b, self.linear)
        if self.h is None:
            # Apply Dirichlet boundary condition to the vector
            fem.petsc.apply_lifting(self.b, [self.bilinear], [self.bcs])
            self.b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
            fem.petsc.set_bc(self.b, self.bcs)

        # Solve linear problem
        self.solver.solve(self.b, self.solution.vector)
        self.solution.x.scatter_forward()  
                
        return self.solution

class transient_thermal_diffusion():
    r"""
    This class implements the solution of the transient thermal diffusion equation with internal power generation made by neutrons.
    
    Parameters
    ----------
    domain : dolfinx.mesh.Mesh
        Mesh imported.
    ct : dolfinx.cpp.mesh.MeshTags_int32
        Cell tags of the regions in the domain.
    ft : dolfinx.cpp.mesh.MeshTags_int32
        Face tags of the boundaries.
    physical_param : dict
        Dictionary with the main physical parameters at reference condition.
    regions_markers : dict
        List of the region markers inside the domain.
    void_bound_marker : int
        Integer describing the index of the void boundary onto which BC is imposed (Dirichlet or Robin).
    h : float, optional (Default = None)
        If not `None`, it, denoted as :math:`h`, represents the heat transfer coefficients in Newton's cooling law.
    TD : float, optional (Default = None)
        If not `None`, it, denoted as :math:`T_D`, represents either the bulk temperature in Newton's cooling law or the Dirichlet condition to impose.
    coupling : str, optional (Default = None)
        If not `None`, it indicates the type of coupling with the temperature: *log* and *sqrt* are available.
            
    """
    def __init__(self, domain: dolfinx.mesh.Mesh, ct: dolfinx.cpp.mesh.MeshTags_int32, ft: dolfinx.cpp.mesh.MeshTags_int32, 
                       physical_param: dict, regions_markers: list, void_bound_marker: int,
                       h : float = None, TD: float = 300., coupling=None):

        self.domain = domain
        self.ft = ft
        self.ct = ct
        self.fdim = domain.geometry.dim - 1
        self.void_bound_marker = void_bound_marker
        
        metadata = {"quadrature_degree": 4}
        self.dx = ufl.Measure('dx', domain=self.domain, subdomain_data=self.ct, metadata=metadata)
        self.ds = ufl.Measure('ds', domain=self.domain, subdomain_data=self.ft, metadata=metadata)

        # Defining the number of energy groups
        self.G = physical_param['Energy Groups']
        
        # Temperature Functional spaces
        self.finiteElement = ufl.FiniteElement("Lagrange", domain.ufl_cell(), 1)
        self.V = FunctionSpace(self.domain, self.finiteElement)
        
        # Material Properties Functional spaces
        self.Qn = FunctionSpace(domain, ("DG", 0))
        if coupling is not None:
            self.Told = Function(self.Qn)
            self.Tref = physical_param['Tref']

        self.coupling = coupling
        self.phys_param = physical_param
        self.regions = regions_markers 

        # Parameters
        self.k = Function(self.Qn) # thermal conductivity
        self.rho_cp = Function(self.Qn) # thermal conductivity
        for idx, regionI in enumerate(self.regions):
            cells = self.ct.find(regionI)
            self.k.x.array[cells]      = self.phys_param['th_cond'][idx]
            self.rho_cp.x.array[cells] = self.phys_param['rho_cp'][idx]
        self.xs_fg = list()
        for g in range(self.G):
            self.xs_fg.append(Function(self.Qn))

        self.Ef = physical_param['Ef']
        self.k_eff = Function(self.Qn)
        self.k_eff.x.set(physical_param['k_eff'])

        # Definig trial and test functions
        self.T = ufl.TrialFunction(self.V)
        self.theta = ufl.TestFunction(self.V)
        
        # Setting up boundary conditions
        self.TD = Function(self.V)
        self.TD.x.set(TD)
        self.h = h # HTC - if None, Dirichlet is imposed
        if h is None:
            self.bcs = [dirichletbc(self.TD, locate_dofs_topological(self.V, self.fdim, self.ft.find(self.void_bound_marker)))]
        
        # Create old function
        self.T_old = Function(self.Qn)
        self.dt = Function(self.V)
        self.T_new = Function(self.V)
        
        # Creating flux function and power form
        self.phi_g = list()
        for g in range(self.G):
            self.phi_g.append(Function(self.V))
        
        self.q3 = self.Ef * sum([self.xs_fg[g] * self.phi_g[g] for g in range(self.G)]) / self.k_eff

    def assembleForm(self, T_ss: dolfinx.fem.Function, dt: float, direct = True):
        r"""
        This function assembles the variational forms and the correspondent linear algebra structures.
        
        Parameters
        ----------
        T_ss : dolfinx.fem.Function
            Steady solution, i.e. the initial condition.
        dt : float
            Time step size.
        direct : bool, optional (Default=True)
            If `True`, a direct linear solver is used (Gauss Elimination), otherwise an iterative one is adopted (preconditioned GMRES with ILU).
            
        """
        self.dt.x.set( dt )
        
        self.left_side   = dot(1. / self.dt * self.T, self.theta) * self.dx
        self.left_side  += dot(self.k / self.rho_cp * grad(self.T), grad(self.theta)) * self.dx
        
        self.right_side  = dot(1. / self.dt * self.T_old, self.theta) * self.dx
        self.right_side += dot(self.q3 / self.rho_cp, self.theta) * self.dx
        
        # Boundary condition (convection)
        if self.h is not None:
            self.left_side  += dot(self.h / self.rho_cp * self.T, self.theta) * self.ds(self.void_bound_marker)
            self.right_side += dot(self.h / self.rho_cp * self.TD, self.theta) * self.ds(self.void_bound_marker)

        self.bilinear  = form(self.left_side)
        self.linear  = form(self.right_side)
        
        self.A = fem.petsc.create_matrix(self.bilinear)
        self.A.zeroEntries()
        if self.h is None:
            fem.petsc.assemble_matrix(self.A, self.bilinear, self.bcs)
        else:
            fem.petsc.assemble_matrix(self.A, self.bilinear)
        self.A.assemble()  
        self.b = fem.petsc.create_vector(self.linear)

        self.solver = PETSc.KSP().create(self.domain.comm)
        self.solver.setOperators(self.A)
        
        if direct:
            self.solver.setType(PETSc.KSP.Type.PREONLY)
            self.solver.getPC().setType(PETSc.PC.Type.LU)  
        else: 
            self.solver.setType(PETSc.KSP.Type.CG)
            self.solver.getPC().setType(PETSc.PC.Type.SOR)   

        # Initialise the old function
        self.T_old.interpolate( fem.Expression(T_ss, self.Qn.element.interpolation_points()))
 
    def update_xs(self):
        """
        This function updates the cross sections according to the coupling scheme chosen.
        """
        
        if self.coupling is None: # no thermal feedback
            for idx, regionI in enumerate(self.regions):
                cells = self.ct.find(regionI)
                for g in range(self.G):
                    self.xs_fg[g].x.array[cells]    = self.phys_param['xs_f'][g][idx]
                    
        elif self.coupling == 'log': # log law
            for idx, regionI in enumerate(self.regions):
                cells = self.ct.find(regionI)
                for g in range(self.G):
                    self.xs_fg[g].x.array[cells]    = self.phys_param['xs_f'][g][0][idx]    + self.phys_param['xs_f'][g][1][idx]    * np.log(self.T_old.x.array[cells] / self.Tref)
                     
        elif self.coupling == 'sqrt': # sqrt law
            for idx, regionI in enumerate(self.regions):
                cells = self.ct.find(regionI)
                for g in range(self.G):
                    self.xs_fg[g].x.array[cells]    = self.phys_param['xs_f'][g][0][idx]    * ( 1 + self.phys_param['xs_f'][g][1][idx]    * ((np.sqrt(self.T_old.x.array[cells]) - np.sqrt(self.Tref))) )
                    
    def advance(self, phi: list):
        r"""
        This function advances in time with a single step.
        
        Parameters
        ----------
        phi : list[dolfinx.fem.Function]
            List of `dolfinx.fem.Function` for the internal power generation, :math:`q''' = E_f\cdot \sum_{g=1}^G \Sigma_{f,g}\phi_g```
        
        Returns
        -------
        solution : dolfinx.fem.Function
            Solution at the new time step..
        """
        # Updating fluxes
        for g in range(self.G):
            if len(self.phi_g[g].x.array[:]) == len(phi[g].x.array[:]):
                self.phi_g[g].x.array[:] = phi[g].x.array[:]
            else:
                self.phi_g[g].interpolate(fem.Expression( phi[g], self.V.element.interpolation_points() ))

        # Update XS
        self.update_xs()
        
        # Assemble RHS with lifting for the Dirichlet BC
        with self.b.localForm() as loc_b:
            loc_b.set(0)
        fem.petsc.assemble_vector(self.b, self.linear)
        if self.h is None:
            # Apply Dirichlet boundary condition to the vector
            fem.petsc.apply_lifting(self.b, [self.bilinear], [self.bcs])
            self.b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
            fem.petsc.set_bc(self.b, self.bcs)

        # Solve linear problem
        self.solver.solve(self.b, self.T_new.vector)
        self.T_new.x.scatter_forward()  
                
        # Update old 
        self.T_old.interpolate( fem.Expression(self.T_new, self.Qn.element.interpolation_points()))
        
        return self.T_new