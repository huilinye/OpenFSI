# OpenFSI Manual
## Table of contents
- About OpenFSI
  - highlighted features
  - code structure
- Install
- Use OpenFSI
  - Prepare input files
  - Application
- Troubleshooting

## About OpenFSI

 OpenFSI is a highly efficient and portable fluid-structure interaction (FSI) simulation package based on immersed-boundary method. The fluid dynamics is accounted by the software [Palabos](http://www.palabos.org/). And the structure solver is implemented within the framework of [LAMMPS](https://lammps.sandia.gov/). In current version, there are 1D, 2D and 3D lattice model and 3D shell model in structure solvers. Using these models, we can model a broad FSI problems including swimming of micro-organisms with tails, flapping of 2D or 3D plate mimicking bird flying and fish swimming and biological flow with large numbers of blood cells.

### Highlighted features

- **Particle based lattice model** The solid is discretized into lattice structure, and the mechanical properties such as stretching and bending are described by series of potential functions that are applied on the lattice nodes.
- **Coupling high-performance software packages** i.e., LAMMPS and Palabos.
- **Highly portable** There are broad applications including 2D and 3D FSI problems.

## How to use it:
 Preparation: 

 1. First you should download Palabos source code from http://www.palabos.org/.
 2. Then you also need to download LAMMPS source code from https://lammps.sandia.gov/.
 3. Make sure you have installed the MPI library.
 4. Add the files in force_coupling and structure_potential into LAMMPS/src directory.

 Compiling:

 1. Compile the LAMMPS as a library. (please reference to https://lammps.sandia.gov/doc/Manual.html)
 2. Modify the Makefile under the eaxmple. 
 - You need to include the LAMMPS library path. 
 - You should make sure the lammps_palabos_coupling and IB_model are in the includePaths.
 - Remeber to include the palabos root directory.
 3. Compile the example.



