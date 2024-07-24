# Installation notes

The solvers has been tested on MacOS and Ubuntu 18.04/20.04 machines with **Python3.10**.

## Dependencies
The *OFELIA* solvers requires the following dependencies:

```python
import numpy
import h5py

import pyvista
import gmsh
import dolfinx

import openmc
```

Be sure to install *gmsh* and *gmsh-api* before *dolfinx=0.6.0* (the package has been tested with real mode of the PETSc library). The instructions to install *dolfinx* are available at [https://github.com/FEniCS/dolfinx#binary](https://github.com/FEniCS/dolfinx#binary).

### Set up a conda environment for *ISF*

At first create a new conda environment
```bash
conda create --name <env_name>
```
If not already done, add conda-forge to the channels
```bash
conda config --add channels conda-forge
```
After having activate it, install 
```bash
conda install python=3.10
```
This provides also *pip* which is necessary to install *gmsh* as
```bash
python -m pip install gmsh gmsh-api
```
To install *FEniCSx* (for real numbers)
```bash
conda install fenics-dolfinx=0.6.0 petsc mpich pyvista
```
Add the following packages
```bash
conda install meshio tqdm
```
Downgrade the following (necessary?)
```bash
python -m pip install setuptools==62.0.0
conda install numpy=1.23.5
```
Once this is completed, it may be necessary to re-install *gmsh*
```bash
python -m pip install gmsh gmsh-api
```

In the end, the *OpenMC* code can be installed following the steps reported below [https://docs.openmc.org/en/stable/usersguide/install.html#installing-on-linux-mac-with-mamba-and-conda-forge](https://docs.openmc.org/en/stable/usersguide/install.html#installing-on-linux-mac-with-mamba-and-conda-forge)
```bash
conda install mamba
mamba search openmc
mamba install openmc
```

## How to use the solvers?

Once all the dependencies have been installed, the classes in `models` folder can be imported and used to solve Neutron Diffusion problems.