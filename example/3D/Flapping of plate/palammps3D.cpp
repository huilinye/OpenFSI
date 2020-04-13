/* This file is part of the Palabos_Lammps coupling program.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * Copyright 2018 Huilin Ye University of Connecticut
 * Author: Huilin Ye (huilin.ye@uconn.edu)
*/
#include "palabos3D.h"
#include "palabos3D.hh"

#include "ibm3D.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "mpi.h"
#include "lammps.h"
#include "input.h"
#include "library.h"
#include "lammpsWrapper.h"

#include "latticeDecomposition.h"
//#include "nearestTwoNeighborLattices3D.h"

using namespace plb;
using namespace std;

typedef double T;
//#define DESCRIPTOR descriptors::ForcedN2D3Q19Descriptor
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//#define DYNAMICS BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())
#define DYNAMICS GuoExternalForceBGKdynamics<T, DESCRIPTOR>(parameters.getOmega())

#define NMAX 150

// initial parameters
	plint Resolution = 0;
	T ReynoldsNumber = 0;
	T Viscosity = 0.;
	T x_length = 0;
	T y_length = 0;
	T z_length = 0;
	
	plint Shear_flag = 0;
	T Shear_Top_Velocity = 0;
	T Shear_Bot_Velocity = 0;
	plint Poiseuille_flag = 0;
	T Poiseuille_bodyforce = 0;
	plint Uniform_flag = 0;
	T U_uniform = 0;
	
	
	plint Total_timestep = 0;
	plint Output_fluid_file = 0;
	plint Output_check_file = 0;
	plint CouplingType = 1;  //1: velocity coupling(default); 2: force coupling

	
	

const T pi = (T)4.*std::atan((T)1.);

template <typename T>
class ShearTopVelocity {
public:
    ShearTopVelocity( IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T(); 
        u[1] = Shear_Top_Velocity; 
        u[2] = T();
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
	
};

template <typename T>
class ShearBottomVelocity {
public:
    ShearBottomVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T();
        u[1] = Shear_Bot_Velocity;
        u[2] = T();
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
	
};



void bishearSetup(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                            IncomprFlowParam<T> const& parameters,
                            OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )  //velocity along y direction
{
	
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D top    = Box3D(0,    nx-1, 0, ny-1, nz-1, nz-1);  //z-direction
    Box3D bottom = Box3D(0,    nx-1, 0, ny-1,    0,    0);
    //Box3D left   = Box3D(1,    nx-2,    0,    0, 1, nz-2);
    //Box3D right  = Box3D(1,    nx-2, ny-1,    ny-1, 1, nz-2);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left, boundary::outflow );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right, boundary::outflow );
    
    setBoundaryVelocity(lattice, top, ShearTopVelocity<T>(parameters,NMAX));
    setBoundaryVelocity(lattice, bottom, ShearBottomVelocity<T>(parameters,NMAX));
	
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(),(T)1.0, Array<T,3>(0.0,0.0,0.0));

    lattice.initialize();
}

void squarePoiseuilleSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                            IncomprFlowParam<T> const& parameters,
                            OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D top    = Box3D(0,    nx-1, 0, ny-1, nz-1, nz-1);  //z direction
    Box3D bottom = Box3D(0,    nx-1, 0, ny-1,    0, 0);
    
    //Box3D inlet    = Box3D(0,    nx-1, 0,    0, 0,    nz-1);  //y direction
    //Box3D outlet = Box3D(0,    nx-1, ny-1,    ny-1, 0, nz-1);
    
    Box3D back   = Box3D(0,    0,    0,    ny-1, 0, nz-1);
    Box3D front  = Box3D(nx-1, nx-1, 0,    ny-1, 0, nz-1);
    
   
    // channel flow
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet);
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, front );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, back );
    
    //setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    //setBoundaryVelocity(lattice, outlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    
    setBoundaryVelocity(lattice, top, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, front, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, back, Array<T,3>((T)0.0,(T)0.0,(T)0.0));
    
        

    //initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX));
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(),(T)1.0, Array<T,3>(0.0,0.0,0.0));

    lattice.initialize();
}

void uniformSetup(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                            IncomprFlowParam<T> const& parameters,
                            OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )  //velocity along y direction
{
	
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    
    Box3D left   = Box3D(0,    nx-1,    0,    0, 0, nz-1);
    Box3D right  = Box3D(0,    nx-1, ny-1,    ny-1, 0, nz-1);
	Box3D front   = Box3D(0,    0,    1,    ny-2, 0, nz-1);
    Box3D back  = Box3D(nx-1,    nx-1, 1,    ny-2, 0, nz-1);
	Box3D top    = Box3D(1,    nx-2, 1, ny-2, nz-1, nz-1);  //z-direction
    Box3D bottom = Box3D(1,    nx-2, 1, ny-2,    0,    0);
	
    
	boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left,boundary::dirichlet);
	setBoundaryVelocity(lattice, left, Array<T,3>(0.0,U_uniform,0.0));
	
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top, boundary::freeslip );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom, boundary::freeslip);
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, front, boundary::freeslip );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, back, boundary::freeslip);
	setBoundaryVelocity(lattice, top, Array<T,3>((T)0.0,(T)U_uniform,(T)0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>((T)0.0,(T)U_uniform,(T)0.0));
	setBoundaryVelocity(lattice, front, Array<T,3>((T)0.0,(T)U_uniform,(T)0.0));
    setBoundaryVelocity(lattice, back, Array<T,3>((T)0.0,(T)U_uniform,(T)0.0));
    
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right, boundary::outflow);
	setBoundaryDensity(lattice,right, (T)1.0);
	  //set right outflow boundary
        /* Box3D globalDomain(lattice->getBoundingBox());
        std::vector<MultiBlock3D*> bcargs;
        bcargs.push_back(lattice);
        bcargs.push_back(rhoBar);
        bcargs.push_back(j);
        T outsideDensity = 1.0;
        int bcType = 1;
        integrateProcessingFunctional(new VirtualOutlet<T,DESCRIPTOR>(outsideDensity, globalDomain, bcType),
                param.outlet, bcargs, 2);
        setBoundaryVelocity(*lattice, param.outlet, uBoundary); */
   

	
	
            
	
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(),(T)1.0, Array<T,3>(0.0,U_uniform,0.0));

    lattice.initialize();
}


void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<float>(*computeDensity(lattice), "density", (T)1.0);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);	
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}

void readParameters(XMLreader const& document)
{
    
    document["geometry"]["Resolution"].read(Resolution);
    document["geometry"]["Viscosity"].read(Viscosity);
    document["geometry"]["x_length"].read(x_length);
    document["geometry"]["y_length"].read(y_length);
    document["geometry"]["z_length"].read(z_length);

    document["fluid"]["Shear_flag"].read(Shear_flag);
    document["fluid"]["Shear_Top_Velocity"].read(Shear_Top_Velocity);
    document["fluid"]["Shear_Bot_Velocity"].read(Shear_Bot_Velocity);
	document["fluid"]["Poiseuille_flag"].read(Poiseuille_flag);
    document["fluid"]["Poiseuille_bodyforce"].read(Poiseuille_bodyforce);
	document["fluid"]["Uniform_flag"].read(Uniform_flag);
    document["fluid"]["U_uniform"].read(U_uniform);
	

    document["simulation"]["Total_timestep"].read(Total_timestep);
    document["simulation"]["Output_fluid_file"].read(Output_fluid_file);
    document["simulation"]["Output_check_file"].read(Output_check_file);

    document["simulation"]["CouplingType"].read(CouplingType);
    
}

   

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

	//read parameters from external file
	string paramXmlFileName;
	paramXmlFileName = "param.xml";
    XMLreader document(paramXmlFileName);
    readParameters(paramXmlFileName);
		
	ReynoldsNumber = 1.0/Viscosity;	
    IncomprFlowParam<T> parameters(
	        1.,
            1.,
            ReynoldsNumber,
			Resolution,
            Resolution,
            x_length,        // lx
            y_length,        // ly
            z_length         // lz
    );
   
    writeLogFile(parameters, "Flow conditions");

    LammpsWrapper wrapper(argv,global::mpi().getGlobalCommunicator());
    char * inlmp = argv[1];
    wrapper.execFile(inlmp);
   
    //MultiTensorField3D<T,3> vel(parameters.getNx(),parameters.getNy(),parameters.getNz());
    pcout<<"Nx,Ny,Nz "<<parameters.getNx()<<" "<<parameters.getNy()<<" "<<parameters.getNz()<<endl;
    LatticeDecomposition lDec(parameters.getNx(),parameters.getNy(),parameters.getNz(),
                              wrapper.lmp);
    SparseBlockStructure3D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
    plint envelopeWidth = 2;

    MultiBlockLattice3D<T, DESCRIPTOR> 
      lattice (MultiBlockManagement3D (blockStructure, threadAttribution, envelopeWidth ),
               defaultMultiBlockPolicy3D().getBlockCommunicator(),
               defaultMultiBlockPolicy3D().getCombinedStatistics(),
               defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
               new DYNAMICS );
    
    //Cell<T,DESCRIPTOR> &cell = lattice.get(550,5500,550);
    pcout<<"dx "<<parameters.getDeltaX()<<" dt  "<<parameters.getDeltaT()<<" tau "<<parameters.getTau()<<endl;
    //pcout<<"51 works"<<endl;

    // set periodic boundary conditions.
	if (Shear_flag == 1){
    lattice.periodicity().toggle(0,true);
    lattice.periodicity().toggle(1,true);
	}
	
	if(Poiseuille_flag == 1){
	lattice.periodicity().toggle(1,true);
	}
	
	/* if (Uniform_flag == 1){
	lattice.periodicity().toggle(0,true);
    lattice.periodicity().toggle(2,true);
	} */
	
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    
	if (Poiseuille_flag == 1)
    squarePoiseuilleSetup(lattice, parameters, *boundaryCondition);   //poiseuille flow boundary condition
    if (Shear_flag == 1)
	bishearSetup(lattice, parameters, *boundaryCondition); //bi-shear flow boundary condition
    	
	if (Uniform_flag == 1){
	uniformSetup(lattice, parameters, *boundaryCondition); 
	}
	
    // Loop over main time iteration.
    util::ValueTracer<T> converge(parameters.getLatticeU(),parameters.getResolution(),1.0e-3);
      //coupling between lammps and palabos
    
    /* for (plint iT=0;iT<4e3;iT++){    //warm up
        lattice.collideAndStream();
    } */
	
    T timeduration = T();
    global::timer("mainloop").start();
	writeVTK(lattice, parameters, 0);
	
    for (plint iT=0; iT<Total_timestep+1; ++iT) {
   
        if (iT%Output_fluid_file ==0 && iT >0){
            pcout<<"Saving VTK file..."<<endl;
            writeVTK(lattice, parameters, iT);
        }
        if (iT%Output_check_file ==0 && iT >0){
            pcout<<"Timestep "<<iT<<" Saving checkPoint file..."<<endl;
            saveBinaryBlock(lattice,"checkpoint.dat");
        }
        // lammps to calculate force
        wrapper.execCommand("run 1 pre no post no");
		
        // Clear and spread fluid force
		if (Shear_flag == 1){
        Array<T,3> force(0,0.,0);
        setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
		}
		
		if (Poiseuille_flag == 1){
        Array<T,3> force(0,Poiseuille_bodyforce,0);
        setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
		}
		
		
		if (CouplingType == 1){
        //-----classical ibm coupling-------------//
        spreadForce3D_fix(lattice,wrapper);
        ////// Lattice Boltzmann iteration step.
        //interpolateVelocity3D_fix(lattice,wrapper);
		lattice.collideAndStream();
         //Interpolate and update solid position
        interpolateVelocity3D_fix(lattice,wrapper);
		}
		else {
        //-----force FSI ibm coupling-------------//
        forceCoupling3D(lattice,wrapper);
        lattice.collideAndStream();
		}
		
    }

    timeduration = global::timer("mainloop").stop();
    pcout<<"total execution time "<<timeduration<<endl;
    delete boundaryCondition;
}
