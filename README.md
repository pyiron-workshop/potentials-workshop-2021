# Potentials Workshop

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pyiron/pyiron_potentialfit/HEAD)

## Installation guidelines

```
conda env update --file binder/environment.yml
```

Configuring the jupyter lab for nglview:

```
jupyter nbextension install nglview --py --sys-prefix
jupyter nbextension enable nglview --py --sys-prefix
jupyter labextension install @jupyter-widgets/jupyterlab-manager --no-build
jupyter labextension install nglview-js-widgets
jupyter labextension install @jupyterlab/toc
```

## Day 1
The scope of the first day is to become familiar with pyiron. We start with LAMMPS as molecular dynamics code to calculate moelcular dynamics trajectories, followed by the introduction of the pyiron tables object, S/PHI/nX as opensource DFT code and Master jobs like the Murnaghan job to calculate the energy volume curve, the calculation of elastic constants and finally the calculation of free energies with phonopy.  

## Day 2
The scope of the second day is to learn about the fitting of interatomic potentials, with a primary focus on atomicrex, the atomic cluster expansion and the neuronal network potentials. 

## Day 3
Finally on the third day we calculate the material properties defined on the first day with the interatomic potentials fitted on the second day and validate those with the results of other interatomic potentials available at atomistictools.org. 

