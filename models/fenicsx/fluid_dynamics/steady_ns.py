import numpy as np

# Mesh generation
from mpi4py import MPI
import dolfinx
from dolfinx import fem
from dolfinx.fem import (Function, FunctionSpace, dirichletbc, locate_dofs_topological, 
                         form, assemble_scalar, locate_dofs_geometrical, apply_lifting)
import ufl
from ufl import grad, div, nabla_grad, dx, inner, dot
from petsc4py import PETSc
from dolfinx.mesh import (CellType, GhostMode, create_rectangle, locate_entities_boundary)


# Function to mark x = 0, x = 1 and y = 0
def noslip_boundary(x):
    return np.logical_or(np.logical_or(np.isclose(x[0], 0.0),
                                       np.isclose(x[0], 1.0)),
                         np.isclose(x[1], 0.0))

# Function to mark the lid (y = 1)
def lid(x):
    return np.isclose(x[1], 1.0)

# Lid velocity
def lid_velocity_expression(x):
    return np.stack((np.ones(x.shape[1]), np.zeros(x.shape[1])))

# Define boundary conditions
class fixedValueVelocity():
    def __init__(self, u_fixed: np.ndarray, t: float, gdim: int):
        
        self.u_fixed = u_fixed
        self.t = t # this may be required for time varying velocity
        self.gdim = gdim

    def __call__(self, x):
        values = np.zeros((self.gdim, x.shape[1]), dtype=PETSc.ScalarType)
        
        for dim in range(self.gdim):
            values[dim] = self.u_fixed[dim]
            
        return values

class ns_steady_nl():
    def __init__(self, domain: dolfinx.mesh.Mesh, ct: dolfinx.cpp.mesh.MeshTags_int32, ft: dolfinx.cpp.mesh.MeshTags_int32,
                       bound_markers: dict, degree_u: int = 2):

        # Domain
        self.domain = domain
        self.ct = ct
        self.ft = ft
        
        self.gdim = self.domain.geometry.dim
        self.fdim = self.gdim - 1

        self.bound_markers = bound_markers

        metadata = {"quadrature_degree": 4}
        self.dx = ufl.Measure('dx', domain=self.domain, subdomain_data=self.ct, metadata=metadata)
        self.ds = ufl.Measure('ds', domain=self.domain, subdomain_data=self.ft, metadata=metadata)

        # Functional Spaces
        self.P2 = ufl.VectorElement("Lagrange", self.domain.ufl_cell(), degree_u)
        self.P1 = ufl.FiniteElement("Lagrange", self.domain.ufl_cell(), degree_u - 1)
        self.mixEl = ufl.MixedElement(self.P2, self.P1)
        self.W = FunctionSpace(self.domain, self.mixEl)

        # Test and trial functions: monolithic
        (self.v, self.q) = ufl.TestFunctions(self.W)
        self.dup = ufl.TrialFunction(self.W)
        self.up  = fem.Function(self.W)
        (self.u, self.p) = ufl.split(self.up)

    def parameters(self, nu: float):

        # Physical Parameters
        self.nu = fem.Constant(self.domain, PETSc.ScalarType(nu))

    def set_cavity_bc(self):
        
        self.lid_velocity = Function(self.W.sub(0).collapse()[0])
        self.lid_velocity.interpolate(lid_velocity_expression)
        self.ft_lid = locate_entities_boundary(self.domain, self.fdim, lid)
        self.dofs_lid = locate_dofs_topological((self.W.sub(0), self.W.sub(0).collapse()[0]), self.fdim, self.ft_lid)
        self.bc_lid = dirichletbc(self.lid_velocity, self.dofs_lid, self.W.sub(0))

        self.no_slip  = Function(self.W.sub(0).collapse()[0])
        self.ft_walls = locate_entities_boundary(self.domain, self.fdim, noslip_boundary)
        self.dofs_walls = locate_dofs_topological((self.W.sub(0), self.W.sub(0).collapse()[0]), self.fdim, self.ft_walls)
        self.bc_w = dirichletbc(self.no_slip, self.dofs_walls, self.W.sub(0))

        self.zero_p = Function(self.W.sub(1).collapse()[0])
        self.zero_p.x.set(0.0)
        self.dofs_p = locate_dofs_geometrical((self.W.sub(1), self.W.sub(1).collapse()[0]), 
                                              lambda x: np.isclose(x.T, [0, 0, 0]).all(axis=1))
        self.bc_p = dirichletbc(self.zero_p, self.dofs_p, self.W.sub(1))

        self.bcs = [self.bc_lid, self.bc_w, self.bc_p]


    def set_bc(self, boundary_type: dict, boundary_value: dict): 
               
        # boundary type can be either 
        #   - 0: inlet, fixed-Value velocity and zero-gradient for pressure
        #   - 1: walls, no slip velocity and zero-gradient for pressure
        #   - 2: outlet, zero-gradient for velocity and null pressure
        
        id_boundaries = list(self.bound_markers.keys())
        
        fixed_velocities = dict()
        dofs_ = dict()
        self.bcs = list()
        for idx, bound in enumerate(id_boundaries):
            
            if np.isclose(boundary_type[bound], 0): # fixedValue velocity
                fixed_velocities[bound] = Function(self.W.sub(0).collapse()[0])
                fixed_velocities[bound].interpolate(fixedValueVelocity(boundary_value[bound], 0., self.gdim))
                dofs_[bound] = locate_dofs_topological((self.W.sub(0), self.W.sub(0).collapse()[0]), self.fdim, 
                                                       self.ft.find(self.bound_markers[bound]))
                self.bcs.append( dirichletbc(fixed_velocities[bound], dofs_[bound], self.W.sub(0)) )
            elif np.isclose(boundary_type[bound], 1): # no slip velocity
                fixed_velocities[bound] = Function(self.W.sub(0).collapse()[0])
                dofs_[bound] = locate_dofs_topological((self.W.sub(0), self.W.sub(0).collapse()[0]), self.fdim, 
                                                       self.ft.find(self.bound_markers[bound]))
                self.bcs.append( dirichletbc(fixed_velocities[bound], dofs_[bound], self.W.sub(0)) )
            else: # outlet
                out_pressure = Function(self.W.sub(1).collapse()[0])
                dofs_[bound] = locate_dofs_topological((self.W.sub(1), self.W.sub(1).collapse()[0]), self.fdim, 
                                                       self.ft.find(self.bound_markers[bound]))
                self.bcs.append( dirichletbc(out_pressure, dofs_[bound], self.W.sub(1) ) )
        
    def create_snes_solution(self) -> PETSc.Vec:  # type: ignore[no-any-unimported]
        """
        Create a petsc4py.PETSc.Vec to be passed to petsc4py.PETSc.SNES.solve.

        The returned vector will be initialized with the initial guess provided in `self._solution`.
        """
        x = self._solution.vector.copy()
        with x.localForm() as _x, self._solution.vector.localForm() as _solution:
            _x[:] = _solution
        return x

    def update_solution(self, x: PETSc.Vec) -> None:  # type: ignore[no-any-unimported]
        """Update `self._solution` with data in `x`."""
        x.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        with x.localForm() as _x, self._solution.vector.localForm() as _solution:
            _solution[:] = _x
            
    def obj_fun(self, snes: PETSc.SNES, x: PETSc.Vec) -> np.float64:
            """Compute the norm of the residual."""
            self.F_fun(snes, x, self._obj_vec)
            return self.b.norm()  # type: ignore[no-any-return]

    def F_fun(self, snes: PETSc.SNES, x: PETSc.Vec, F_vec: PETSc.Vec) -> None:
            """Assemble the residual."""
            self.update_solution(x)
            with F_vec.localForm() as F_vec_local:
                F_vec_local.set(0.0)
            fem.petsc.assemble_vector(F_vec, self._F)
            apply_lifting(F_vec, [self._J], [self.bcs], x0=[x], scale=-1.0)
            F_vec.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
            fem.set_bc(F_vec, self.bcs, x, -1.0)

    def J_fun( self, snes: PETSc.SNES, x: PETSc.Vec, J_mat: PETSc.Mat, P_mat: PETSc.Mat) -> None:
            """Assemble the jacobian."""
            J_mat.zeroEntries()
            fem.petsc.assemble_matrix(J_mat, self._J, self.bcs, diagonal=1.0)
            J_mat.assemble()

    def assemble(self, maxIter = 20, verbose = False):

        # Variational forms
        self.F = (self.nu * inner(grad(self.u), grad(self.v)) * self.dx
                + inner(grad(self.u) * self.u, self.v) * self.dx
                - inner(self.p, ufl.div(self.v)) * self.dx
                + inner(div(self.u), self.q) * self.dx)
        self.J = ufl.derivative(self.F, self.up, self.dup)

        self._F = form(self.F)
        self._J = form(self.J)

        # Create matrix and vector
        self._solution = self.up
        self._obj_vec = fem.petsc.create_vector(self._F)
        self.b = fem.petsc.create_vector(self._F)
        self.A = fem.petsc.create_matrix(self._J)

        # Solver settings
        self.solver = PETSc.SNES().create(self.domain.comm)
        self.solver.setTolerances(max_it=maxIter)
        
        self.solver.getKSP().setType("preonly")
        # self.solver.getKSP().setType("gmres")
        self.solver.getKSP().getPC().setType("lu")
        self.solver.getKSP().getPC().setFactorSolverType("mumps")

        self.solver.setObjective(self.obj_fun)
        self.solver.setFunction(self.F_fun, self.b)
        self.solver.setJacobian(self.J_fun, J=self.A, P=None)
        if verbose:
            self.solver.setMonitor(lambda _, it, residual: print(it, residual))

    def solve(self):
        up_copy = self.create_snes_solution()
        self.solver.solve(None, up_copy)
        self.update_solution(up_copy)
        return self._solution