# Introduction

**Authors**: Lorenzo Loi & Stefano Riva, Carolina Introini, Antonio Cammi

**OFELIA: Openmc-FEnicsx for muLtiphysics tutorIAl**

The aim of the repository consists in collecting tutorials for MultiPhysics problems in Nuclear Reactors Engineering implemented using open-source packages for Python. Specifically, two different aspects are going to be considered: neutronics and thermal-hydraulics, fundamental topics for nuclear reactors.
Two packages are used: [OpenMC](https://docs.openmc.org/en/stable/), a Monte Carlo code for stationary neutron transport, and [FEniCSx (v.0.6.0)](https://fenicsproject.org/), a Finite Element library able to solve both multi-group neutron diffusion and thermal-hydraulics problems (governed by heat diffusion or Navier-Stokes, mainly laminar flow). These packages are going to be used both singularly for different case studies and coupled in a MultiPhysics scheme for a specific case study of fundamental interest for nuclear reactor applications.

This documentation includes a brief introduction to the neutronics adn thermal-hydrualics models, the numerical methods adopted, some basic theory for Finite Elements, the API documentation of some solvers and some basic tutorials on how to use OpenMC and FEniCSx for Nuclear Reactors applications; moreover, the code implementing MP case of {cite:p}`OFELIA_2024` is reported.

This work has been carried out at the [Nuclear Reactors Group - ERMETE Lab](https://github.com/ERMETE-Lab) at [Politecnico di Milano](https://polimi.it), under the supervision of Prof. Antonio Cammi.

---

## How to cite

If you use this code in your research, please cite the following paper

- Loi, L., Riva, S., Introini, C., Giacobbo, F., Wang, X., & Cammi, A. (2024). OFELIA: an OpenMC-FEniCSx Coupling for Neutronic Calculation with Temperature Feedback. accepted at Nuclear Engineering and Design.

Here reported in `.bib` format
```latex
@article{OFELIA_2024,
title = {{OFELIA: an OpenMC-FEniCSx Coupling for Neutronic Calculation with Temperature Feedback}},
journal = {accepted at Nuclear Engineering and Design},
volume = {},
pages = {},
year = {2024},
issn = {},
doi = {},
url = {},
author = {Lorenzo Loi and Stefano Riva and Carolina Introini and Francesca Giacobbo and Xiang Wang and Antonio Cammi},
}
```

