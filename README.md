# MultiPhysics in Python - OFELIA

**OFELIA: Openmc-FEnicsx for muLtiphysics tutorIAl**
**Authors**: Lorenzo Loi, Stefano Riva, Carolina Introini, Antonio Cammi

This repository collects the supporting code for OFELIA (Openmc-FEnicsx for muLtiphysics tutorIAl): an open-source tool to simulate multi-physics nuclear cases adopting [OpenMC (v 0.13.2)](https://openmc.org/) and [FEniCSx (v. 0.6.0)](https://fenicsproject.org/) for Python. Moreover, useful tutorials are available for the single-physics codes OpenMC and FEniCSx.

The aim of this repository consists in collecting tutorials for FEniCSx and OpenMC, and their coupling version (OFELIA), specific to nuclear reactors applications. In particular, the former FEniCSx collects the following:

- Fluid Dynamics: cavity, backward facing step and laminar flow over cylinder.
- Neutronics (MultiGroup Diffusion): ANL benchmarks and MultiPhysics application

whereas, the OpenMC tutorial implements the simulation of a TRIGA reactor.

## How to cite
If you use this set of codes in your research, please cite the following paper

- L. Loi, S. Riva, C. Introini, F. Giacobbo, X. Wang, and A. Cammi, “OFELIA: An OpenMC-FEniCSx coupling for neutronic calculation with temperature feedback,” Nuclear Engineering and Design, vol. 428, p. 113480, 2024.

Here reported in `.bib` format
```latex
@article{LOI2024113480,
title = {OFELIA: An OpenMC-FEniCSx coupling for neutronic calculation with temperature feedback},
journal = {Nuclear Engineering and Design},
volume = {428},
pages = {113480},
year = {2024},
issn = {0029-5493},
doi = {https://doi.org/10.1016/j.nucengdes.2024.113480},
url = {https://www.sciencedirect.com/science/article/pii/S0029549324005806},
author = {Lorenzo Loi and Stefano Riva and Carolina Introini and Francesca Giacobbo and Xiang Wang and Antonio Cammi},
keywords = {OpenMC, FEniCSx, Nuclear reactor, Monte Carlo, Thermal hydraulics},
}
```

## Contact Information

If interested, contact lorenzo.loi@polimi.it or stefano.riva@polimi.it
