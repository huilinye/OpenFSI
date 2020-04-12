# OpenFSI Manual

## Table of contents
- [About OpenFSI](#About-OpenFSI)
  - [highlighted features](#Highlighted-features)
  - [code structure](#Code-structure)
- [Compile and Run](#Compile-and-Run)
  - [Preparation](#Preparation)
  - [Compile](#Compile)
  - [Run](#Run)
- [Use OpenFSI](#Use-OpenFSI)
  - Prepare input files
  - Application
- Troubleshooting

## About OpenFSI

 OpenFSI is a highly efficient and portable fluid-structure interaction (FSI) simulation package based on immersed-boundary method. The fluid dynamics is accounted by the software [Palabos](http://www.palabos.org/). And the structure solver is implemented within the framework of [LAMMPS](https://lammps.sandia.gov/). In current version, there are 1D, 2D and 3D lattice model and 3D shell model in structure solvers. Using these models, we can model a broad FSI problems including swimming of micro-organisms with tails, flapping of 2D or 3D plate mimicking bird flying and fish swimming and biological flow with large numbers of blood cells.

### Highlighted features

- **Particle based lattice model** The solid is discretized into lattice structure, and the mechanical properties such as stretching and bending are described by series of potential functions that are applied on the lattice nodes.
- **Coupling high-performance software packages** i.e., LAMMPS and Palabos.
- **Highly portable** There are broad applications including 2D and 3D FSI problems.

### Code structure
- `example`: Examples to show how OpenFSI runs 2D and 3D problems
- `src/fix_LB`: Revised `fix_lb` files embedded in LAMMPS for solving fluid dynamics using Latticel Boltzmann method (LBM)
- `src/IB_interface`: Implementation for immersed-boundary method including velocity and force couplings
- `src/lammps_palabos_coupling`: Technique to fulfill the consistant cpu mapping between LAMMPS and Palabos
- `src/main`: Main file to offer user interface to set up properties such as type of flow and boundary conditions
- `src/structure_potential`: Potentials to calculate lattice model and membrane model

## Compile and Run 
 
 ### Preparation: 
 - First you should download Palabos source code from http://www.palabos.org/.
 - Then you also need to download LAMMPS source code from https://lammps.sandia.gov/.
 - Make sure you have installed the MPI library.
 - Add the files in `src/fix_LB` and `src/structure_potential` into LAMMPS/src directory.

 ### Compiling:

- Compile the LAMMPS as a library

  Assume you install the lammps under the dirctory `$lammps_dir`
  Inside the directory `$lammps_dir/src/MAKE`, there are many options for the makefile. It is possibile to run lammps in serial
  , mpi or intel-optimised version, which depends on the computer architecture used in the simulation.
  For example, if you want to run the simulation in mpi mode, you should first observe that there is a makefile
  named `Makefile.mpi` in `$lammps_dir/src/MAKE`, then make it with library mode
  
  ```
  cd $lammps_dir/src
  make mode=lib mpi
  ```
  
  This will generate a lammps library named `liblammps_mpi.a`

- Install coupling interface 

  Create directory, e.g., `coupling_dir` to include the `src/IB_interface` and `src/lammps_palabos_coupling` files
  
- Modify Makefile

  In the `example` directory, each case is attached a Makefile for generating executable file. 
  First, you should provide the directory where Palabos locates, e.g., `palabos_dir`
 
  `palabosRoot  = $palabos_dir`
 
  You can set up the name of the project files for concrete problem. 2D and 3D prototypes are provided in `src/main`
 
  `projectFiles = $project_name.cpp`
 
  The libraries in lammps dirctory should be included
 
  `libraryPaths = $lammps_dir/src`
 
   Also, the dirctory containing the files in lammps and coupling is necessary
 
  `includePaths =  $lammps_dir/src $coupling_dir `
 
   Finally, the name of shared library generated in lammps is provided
 
   `libraries    = liblammps_mpi.a`
 
    It is strongly recommended not to change any setups including the compiler to use with MPI parallelism, general compiler flags and palabos compile setups.
 
- Compile

  Put the Makefile in your working dirctory `work_dir`, and then simply type
 
  `make`
 
   in the command window, a executable file named `$project_name` will generate
 
### Run
 
   Configure you lammps input file, e.g., `in.lammps`, which setups the system, including reading particle coordinates, angles information and 
   bond information, and parallelism pattern such as CPU cores that will be used to run the simulation. After checking 
   the input file, just type following command to run the simulation
 
   `./$project_name in.lammps > log.file`

   `log.file` will record the information displayed in the window.

## Use OpenFSI

###Prepare data and input files
One data file is necessary that describes the system including the number of the particles and bonds, the size of the domain,
coordinates of the particles, and coefficients of the bonds, etc. Following is an example from `example/2D`:

```
866 atoms
4572 bonds
0 angles
0 dihedrals
762 impropers
```