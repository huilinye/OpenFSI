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
#include "palabos2D.h"
#include "palabos2D.hh"

#include "ibm2D.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "mpi.h"
#include "lammps.h"
#include "input.h"
#include "library.h"
#include "lammpsWrapper.h"

#include "latticeDecomposition2D.h"
//#include "nearestTwoNeighborLattices3D.h"

using namespace plb;
using namespace std;

typedef double T;
//#define DESCRIPTOR descriptors::ForcedN2D3Q19Descriptor
#define DESCRIPTOR descriptors::ForcedD2Q9Descriptor
//#define DYNAMICS BGKdynamics<T, DESCRIPTOR>(parameters.getOmega())
#define DYNAMICS GuoExternalForceBGKdynamics<T, DESCRIPTOR>(parameters.getOmega())

#define NMAX 150

// initial parameters
	plint Resolution = 0;
	T ReynoldsNumber = 0;
	T Viscosity = 0.;
	T x_length = 0;
	T y_length = 0;
	//T z_length = 0;
	
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
    plint StaticAtomType = 0;  //fix atom with specific type
	
	

const T pi = (T)4.*std::atan((T)1.);

/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
	T Ly = parameters.getNy()-1;
    return 0.06*6*y*(Ly-y)/(Ly*Ly);
}

/// Linearly decreasing pressure profile
T poiseuillePressure(plint iX, IncomprFlowParam<T> const& parameters) {
    T Lx = parameters.getNx()-1;
    T Ly = parameters.getNy()-1;
    return 0.06*12.*parameters.getLatticeNu()/ (Ly*Ly) * (Lx-(T)iX);
}

/// Convert pressure to density according to ideal gas law
T poiseuilleDensity(plint iX, IncomprFlowParam<T> const& parameters) {
    return poiseuillePressure(iX,parameters)*DESCRIPTOR<T>::invCs2 + (T)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to initialize a pressure boundary to constant density
template<typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_)
        : density(density_)
    { }
    T operator()(plint iX, plint iY) const {
        return density;
    }
private:
    T density;
};

/// A functional, used to create an initial condition for the density and velocity
template<typename T>
class PoiseuilleVelocityAndDensity {
public:
    PoiseuilleVelocityAndDensity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        rho = poiseuilleDensity(iX,parameters);
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

template <typename T>
class ShearTopVelocity {
public:
    ShearTopVelocity( IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, Array<T,2>& u) const  {
        
        u[0] = Shear_Top_Velocity; 
		u[1] = T(); 
        //u[2] = T();
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
    void operator()(plint iX, plint iY, Array<T,2>& u) const  {
        
        u[0] = Shear_Bot_Velocity;
		u[1] = T();
        //u[2] = T();
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
	
};



void bishearSetup(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                            IncomprFlowParam<T> const& parameters,
                            OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )  //velocity along y direction
{
	
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    //const plint nz = parameters.getNz();
    Box2D top    = Box2D(0,    nx-1, ny-1, ny-1);  //y-direction
    Box2D bottom = Box2D(0,    nx-1,    0,    0);
    //Box3D left   = Box3D(1,    nx-2,    0,    0, 1, nz-2);
    //Box3D right  = Box3D(1,    nx-2, ny-1,    ny-1, 1, nz-2);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left, boundary::outflow );
    //boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right, boundary::outflow );
    
    setBoundaryVelocity(lattice, top, ShearTopVelocity<T>(parameters,NMAX));
    setBoundaryVelocity(lattice, bottom, ShearBottomVelocity<T>(parameters,NMAX));
	
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(),(T)1.0, Array<T,2>(0.0,0.0));

    lattice.initialize();
}

void poiseSetup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters,
                    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    Box2D outlet(nx-1,nx-1, 1, ny-2);

    // Create Velocity boundary conditions everywhere
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, 0, 1, ny-2) );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, 0, 0) );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, ny-1, ny-1) );
    // .. except on right boundary, where we prefer an outflow condition
    //    (zero velocity-gradient).
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(nx-1, nx-1, 1, ny-2), boundary::outflow );

    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );
    setBoundaryDensity (
            lattice, outlet,
            ConstantDensity<T>(1.) );
			
    initializeAtEquilibrium (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocityAndDensity<T>(parameters) );
    lattice.initialize();
}

void poise_force_Setup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters,
                    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    Box2D top    = Box2D(0,    nx-1, ny-1, ny-1);  //y-direction
    Box2D bottom = Box2D(0,    nx-1,    0,    0);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );

    setBoundaryVelocity(lattice, top, Array<T,2>((T)0.0,(T)0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,2>((T)0.0,(T)0.0));
	initializeAtEquilibrium(lattice, lattice.getBoundingBox(),(T)1.0, Array<T,2>(0.0,0.0));
	
    lattice.initialize();
}


void uniformSetup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters,
                    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    Box2D outlet(nx-1,nx-1, 0, ny-1);
    Box2D inlet(0,0, 0, ny-1);
	Box2D top    = Box2D(1,    nx-2, ny-1, ny-1);  //y-direction
    Box2D bottom = Box2D(1,    nx-2,    0,    0);
	
	boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top, boundary::freeslip );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom, boundary::freeslip );
	boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet );
	boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet, boundary::normalOutflow );
    
    
    
    setBoundaryVelocity(lattice, inlet,Array<T,2>(U_uniform,0.0));
    
    initializeAtEquilibrium (
            lattice, lattice.getBoundingBox(),
            (T)1.0, Array<T,2>(U_uniform,0.0) );
			
			
			
    lattice.initialize();
}

void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<float>(*computeDensity(lattice), "density", (T)1.0);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);	
    //vtkOut.writeData<2,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}

void readParameters(XMLreader const& document)
{
    
    document["geometry"]["Resolution"].read(Resolution);
    document["geometry"]["Viscosity"].read(Viscosity);
    document["geometry"]["x_length"].read(x_length);
    document["geometry"]["y_length"].read(y_length);
    //document["geometry"]["z_length"].read(z_length);

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
	document["simulation"]["StaticAtomType"].read(StaticAtomType);
    
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
	        1.0,
            ReynoldsNumber,
			Resolution,
            x_length,        // lx
            y_length        // ly
    );
   
    //writeLogFile(parameters, "Flow conditions");

    LammpsWrapper wrapper(argv,global::mpi().getGlobalCommunicator());
    char * inlmp = argv[1];
    wrapper.execFile(inlmp);
   
    //MultiTensorField3D<T,3> vel(parameters.getNx(),parameters.getNy(),parameters.getNz());
    pcout<<"Nx,Ny "<<parameters.getNx()<<" "<<parameters.getNy()<<endl;
    LatticeDecomposition2D lDec(parameters.getNx(),parameters.getNy(), 1.0,
                              wrapper.lmp);
    SparseBlockStructure2D blockStructure = lDec.getBlockDistribution();
    ExplicitThreadAttribution* threadAttribution = lDec.getThreadAttribution();
    plint envelopeWidth = 2;

    MultiBlockLattice2D<T, DESCRIPTOR> 
      lattice (MultiBlockManagement2D (blockStructure, threadAttribution, envelopeWidth ),
               defaultMultiBlockPolicy2D().getBlockCommunicator(),
               defaultMultiBlockPolicy2D().getCombinedStatistics(),
               defaultMultiBlockPolicy2D().getMultiCellAccess<T,DESCRIPTOR>(),
               new DYNAMICS );
    
    //Cell<T,DESCRIPTOR> &cell = lattice.get(550,5500,550);
    pcout<<"dx "<<parameters.getDeltaX()<<" dt  "<<parameters.getDeltaT()<<" tau "<<parameters.getTau()<<endl;
    //pcout<<"51 works"<<endl;

    // set periodic boundary conditions.
	if (Shear_flag == 1){
    lattice.periodicity().toggle(0,true);
    //lattice.periodicity().toggle(1,true);
	}
	
	//if(Poiseuille_flag == 1){
	//lattice.periodicity().toggle(0,true);
	//}
	
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition2D<T,DESCRIPTOR>();
    
	
    if (Shear_flag == 1)
	bishearSetup(lattice, parameters, *boundaryCondition); //bi-shear flow boundary condition
    
	if (Poiseuille_flag == 1)
	poiseSetup(lattice, parameters, *boundaryCondition);  //velocity distribution
    //if (Poiseuille_flag == 1)
	//poise_force_Setup(lattice, parameters, *boundaryCondition);

    if (Uniform_flag == 1)
	uniformSetup(lattice, parameters, *boundaryCondition); 

    

    // Loop over main time iteration.
    util::ValueTracer<T> converge(parameters.getLatticeU(),parameters.getResolution(),1.0e-3);
      //coupling between lammps and palabos
    
    /* for (plint iT=0;iT<4e3;iT++){    //warm up
        lattice.collideAndStream();
    } */
	
    T timeduration = T();
    global::timer("mainloop").start();
	//writeVTK(lattice, parameters, 0);
	
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
        //wrapper.execCommand("run 1 pre no post no");
		wrapper.execCommand("run 1");
		
        // Clear and spread fluid force
		if (Shear_flag == 1){
        Array<T,2> force(0,0.);
        setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
		}
		
		//if (Poiseuille_flag == 1){
        //Array<T,2> force(Poiseuille_bodyforce,0);
        //setExternalVector(lattice,lattice.getBoundingBox(),DESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
		//}
		
		if (CouplingType == 1){
        //-----classical ibm coupling-------------//
        spreadForce2D_fix(lattice,wrapper,StaticAtomType);
        ////// Lattice Boltzmann iteration step.
        //interpolateVelocity3D_fix(lattice,wrapper);
		lattice.collideAndStream();
         //Interpolate and update solid position
        interpolateVelocity2D_fix(lattice,wrapper,StaticAtomType);
		}else{
			forceCoupling2D(lattice,wrapper);
			lattice.collideAndStream();
		}
		
		
    }

    timeduration = global::timer("mainloop").stop();
    pcout<<"total execution time "<<timeduration<<endl;
    delete boundaryCondition;
}
