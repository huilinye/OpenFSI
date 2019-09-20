# OpenFSI

 Preparation: 

 1. First you should download Palabos source code from http://www.palabos.org/.
 2. Then you also need to download LAMMPS source code from https://lammps.sandia.gov/.
 3. Make sure you have installed the MPI library.
 4. Add the files in force_coupling and structure_potential into LAMMPS/src directory.

 Compiling:

 1. Compile the LAMMPS as a library. (please reference to https://lammps.sandia.gov/doc/Manual.html)
 2. Modify the Makefile under the eaxmple. 
 2.1 You need to include the LAMMPS library path. 
 2.2 You should make sure the lammps_palabos_coupling and IB_model are in the includePaths.
 2.3 Remeber to include the palabos root directory.
 3. Compile the example.
 


