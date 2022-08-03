# SG-Eady-slice
MATLAB functions for solving the semi-geostrophic Eady slice equations using the geometric method.

## Getting started ##
To use these MATLAB functions you must first install

* [MATLAB-Voro](https://github.com/smr29git/MATLAB-Voro)

## Examples ##

See the MATLAB live script ``Examples.mlx`` and the corresponding PDF file ``Examples.pdf``. We have tested the code with several initial conditions including: random initial seeds; a stable perturbation of a steady shear flow with up to 1,000 seeds; an unstable perturbation of a steady shear flow with up to 7,000 seeds. 

## Licence ##

See LICENCE.md

## Software limitations ##

* The code is limited to 3D (2D code coming soon).
* The source measure is the Lebesgue measure on a cuboid.
* The support of the discrete target measure must be contained in the support of the source measure.
* The transport cost is either the quadratic cost or the periodic quadratic cost.

## Main contributors ##

* Charlie Egan, Heriot-Watt University and the Maxwell Institute for Mathematical Sciences
* [Steve Roper](https://www.gla.ac.uk/schools/mathematicsstatistics/staff/stevenroper/#), University of Glasgow
* [David Bourne](http://www.macs.hw.ac.uk/~db92/), Heriot-Watt University and the Maxwell Institute for Mathematical Sciences

This repository was created to accompany the following paper:

* Egan, C.P., Bourne, D.P., Cotter, C.J., Cullen, M.J.P., Pelloni, B., Roper, S.M. & Wilkinson, M. (to appear) A new implementation of the geometric method for solving the Eady slice equations, *Journal of Computational Physics*

Please consider citing this paper if you find our code useful.
