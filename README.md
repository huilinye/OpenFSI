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
 
  ```palabosRoot  = $palabos_dir```
 
  You can set up the name of the project files for concrete problem. 2D and 3D prototypes are provided in `src/main`
 
  ```projectFiles = $project_name.cpp```
 
  The libraries in lammps dirctory should be included
 
  ```libraryPaths = $lammps_dir/src```
 
   Also, the dirctory containing the files in lammps and coupling is necessary
 
  ```includePaths =  $lammps_dir/src $coupling_dir ```
 
   Finally, the name of shared library generated in lammps is provided
 
   ```libraries    = liblammps_mpi.a```
 
    It is strongly recommended not to change any setups including the compiler to use with MPI parallelism, general compiler flags and palabos compile setups.
 
- Compile

  Put the Makefile in your working dirctory `work_dir`, and then simply type
 
  ```make```
 
   in the command window, a executable file named `$project_name` will generate
 
### Run
 
   Configure you lammps input file, e.g., `in.lammps`, which setups the system, including reading particle coordinates, angles information and 
   bond information, and parallelism pattern such as CPU cores that will be used to run the simulation. After checking 
   the input file, just type following command to run the simulation
 
   ```./$project_name in.lammps > log.file```

   `log.file` will record the information displayed in the window.

## Use OpenFSI

### Prepare data and input files
- Data file

One data file is necessary that describes the system including the number of the particles and bonds, the size of the domain,
coordinates of the particles, and coefficients of the bonds, etc. Following is an example from `example/2D/2D_cylinder_beam.data`:

The first part claims the numbers and types of particles and relevant lattice properties, respectively.
```
866 atoms
4572 bonds
0 angles
0 dihedrals
762 impropers

866 atom types
4572 bond types
0 angle types
0 dihedral types
762 improper types
```
The second part define the boundary of the domain
```
0.000 240.000 xlo xhi
0.000 120.000 ylo yhi
0.000   1.000 zlo zhi
```
`lo` and `hi` represent the lower and upper bound of the domain. Because there is no physical 2D domain setup in LAMMPS, the third dimension should be claimed no matter what is 
the dimension of the problem. It is necessary to
exert 2D command in the lammps input file to exclude the effect of the third dimension.

Following is the mass part
```
Masses

1 0.04375914 
2 0.06035109 
3 0.02000890 
...
```
The first column is the type of the particle and the second one is the mass of the corresponding particle.

Then the Bond, angle, dihedral and improper coefficients are presented one by one. In bond part
```
Bond Coeffs #harmonic

1 0.00061920 0.00000000
2 0.00061805 0.00000000
3 0.00061692 0.00000000
...
```
Here, the harmonic bond is shown.
The first column is the type of the bond; second column tells the bond coefficient; and last column is the equlibrium length of the bond.
For this 2D problem, there is no angle and dihedral potential applied, therefore these parts can be leaved blank in the
data file. The improper part is like
```
Improper Coeffs #neohookean

1 0.010000 1.000000 0.795583
2 0.010000 1.000000 0.795640
3 0.010000 1.000000 0.795695
...
``` 
Following is the coordinates of the particles in the system
```
1 1 1 40.00000000 60.00000000 0.50000000
2 2 2 60.00000000 60.00000000 0.50000000
3 2 3 130.00000000 60.00000000 0.50000000
...
```
First column is the ID of the particle; second and third column are type ID and molecule ID, respectively. Here,
molecule ID is illustrated to represent the body of the solid structure. For example, if there are two beams in the
system, then one is denoted molecule 1 and the other is described as molecule 2. The third, fourth and fifth columns are
the coordinate of the corresponding particle in x, y and z directions, respectively.

The last part is the patterns of the lattice properties, e.g., bond, angle, and improper. We give an example of bond
```
Bonds

1 1 1 7
2 2 7 8
3 3 8 9
...
```
The first column is the ID of the bond, and the second one is the bond type. The third and fourth are the IDs of the particles which are
are bonded.

- Input files

   Two input files are required: one having to do with the setups of the LAMMPS, and the other is designed to 
   control the flow conditions.

   - LAMMPS input file
    
    Following is a typical LAMMPS input file
    ```
    units          lj
    dimension 3
    
    boundary    p p p
    processors 4 2 1
    bond_style     harmonic
    improper_style neohookean

    read_data 2D_cylinder_beam.data

    timestep 0.010
    dump 1 all custom 100  beam.lammpstrj id mol type x y z fx fy fz
	```
	The `units` gives the unit of this system. Here, `lj` unit is non-dimensionalized unit. Each parameter in
	the LAMMPS will be normalized by the reference. Other units can be found in LAMMPS manual. `dimension` defines
	the dimension of the problem. It is recommended to set it to 3. `boundary` points to the boundary of the particle
	system. Note that this is not the boundary of the flow field. The MPI pattern can be set using the `processors`
	command, which is followed by the number of processors used in x, y, z directions, respectively. Before reading
	the data file `read_data 2D_cylinder_beam.data`, the bond style and improper style should be defined. Also, the timestep
	is necessary to setup here. Finally, the command `dump` can output the particle system with LAMMPS trajectory file, which
	contains the particle ID, type and coordinates.
	
   - Flow input file
   
    The flow conditions are contained in the xml file named `param.xml` for the initialization of Palabos. There are three parts in the flow input file.
	The first one is the `<geometry> ` part. It gives the domain size used in the flow field. It is usually the same as
	the size in LAMMPS system like
	```
	<geometry>
    <Viscosity> 0.02 </Viscosity>
    <!-- Reynolds number -->
    <x_length> 240 </x_length>
    <y_length> 120 </y_length>
    <!-- size of the flow domain -->
    </geometry>
	```
	Also, you can define the viscosity of the fluid model that the Reynolds number can be fixed. 
	
	The second
	part is the `<fluid>` part. 
	```
    <fluid>
    <Shear_flag> 0 </Shear_flag>
	<!-- shear flow flag -->
    <Shear_Top_Velocity> 0.1 </Shear_Top_Velocity>
    <Shear_Bot_Velocity> -0.1 </Shear_Bot_Velocity>
	<Poiseuille_flag> 1 </Poiseuille_flag>
	<!-- Poiseuille flow flag -->
	<Poiseuille_bodyforce> 0.000015 </Poiseuille_bodyforce>
	<Uniform_flag> 0 </Uniform_flag>
	<!-- uniform flow flag -->
	<U_uniform> 0.01 </U_uniform>
    </fluid>
	```
	It provides three flow type: simple shear flow, Poiseuille flow and uniform flow. If simple shear 
	flow is applied, the velocities of the top and bottom wall should be explicitly provided. Also, if it is Poiseuille flow,
	the bodyforce should be presented. The uniform flow can be characterized by the uniform velocity.
	
	The Third part is `<simulation>` part. It shows the simulation time and output control.
	```
	<simulation>
    
    <Total_timestep>  400000 </Total_timestep>
    <!-- Maximum timestep for simulation. -->
    <Output_fluid_file> 200      </Output_fluid_file>
    <Output_check_file>  20000  </Output_check_file>

    <CouplingType>      2   </CouplingType>
    <!-- CouplingType: 1: velocity coupling(default); 2: force coupling -->
	
    </simulation>
	```
	Also, in this part, you can determine which type of the IB method will be used. 1 represents the velocity coupling,
	and 2 is the force coupling.