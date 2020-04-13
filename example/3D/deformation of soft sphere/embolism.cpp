/* This file is part of the Palabos library.
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
	T x_length = 0;
	T y_length = 0;
	T z_length = 0;
	
	plint Shear_flag = 0;
	T Shear_Top_Velocity = 0;
	T Shear_Bot_Velocity = 0;
	plint Poiseuille_flag = 0;
	T Poiseuille_bodyforce = 0;
	
	plint Total_timestep = 0;
	plint Output_fluid_file = 0;
	plint Output_check_file = 0;
	plint CouplingType = 1;  //1: velocity coupling(default); 2: force coupling

	
	

const T pi = (T)4.*std::atan((T)1.);

static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN){
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T nu = parameters.getLatticeNu();
    const T uMax = parameters.getLatticeU();

    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }
    for (plint iN = 1; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= ((T)1 / (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
    }

    T alpha = -(T)8 * uMax * pi * pi * pi / (a*a*(pi*pi*pi-(T)32*sum)); // alpha = -dp/dz / mu
    T deltaP = - (alpha * nu);
    return deltaP;
}

T poiseuilleVelocity(plint iX, plint iY, IncomprFlowParam<T> const& parameters, plint maxN){
    const T a = parameters.getNx()-1;
    const T b = parameters.getNy()-1;

    const T x = (T)iX - a / (T)2;
    const T y = (T)iY - b / (T)2;

    const T alpha = - poiseuillePressure(parameters,maxN) / parameters.getLatticeNu();

    T sum = T();

    for (plint iN = 0; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }
    for (plint iN = 1; iN < maxN; iN += 2){
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
             / ( std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }

    sum *= ((T)4 * alpha * a *a /std::pow(pi,(T)3));
    sum += (alpha / (T)2 * (x * x - a*a / (T)4));
    
    return sum;
}

template <typename T>
class SquarePoiseuilleDensityAndVelocity {
public:
    SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

template <typename T>
class SquarePoiseuilleVelocity {
public:
    SquarePoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
        : parameters(parameters_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
        u[0] = T();
        u[1] = T();
        u[2] = poiseuilleVelocity(iX, iY, parameters, maxN);
    }
private:
    IncomprFlowParam<T> parameters;
    plint maxN;
};

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
    

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );
    
    
    
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

T computeRMSerror ( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters )
{
    MultiTensorField3D<T,3> analyticalVelocity(lattice);
    setToFunction( analyticalVelocity, analyticalVelocity.getBoundingBox(),
                   SquarePoiseuilleVelocity<T>(parameters, NMAX) );
    MultiTensorField3D<T,3> numericalVelocity(lattice);
    computeVelocity(lattice, numericalVelocity, lattice.getBoundingBox());

           // Divide by lattice velocity to normalize the error
    return 1./parameters.getLatticeU() *
           // Compute RMS difference between analytical and numerical solution
           std::sqrt( computeAverage( *computeNormSqr(
                          *subtract(analyticalVelocity, numericalVelocity)
                     ) ) );
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
    document["geometry"]["ReynoldsNumber"].read(ReynoldsNumber);
    document["geometry"]["x_length"].read(x_length);
    document["geometry"]["y_length"].read(y_length);
    document["geometry"]["z_length"].read(z_length);

    document["fluid"]["Shear_flag"].read(Shear_flag);
    document["fluid"]["Shear_Top_Velocity"].read(Shear_Top_Velocity);
    document["fluid"]["Shear_Bot_Velocity"].read(Shear_Bot_Velocity);
	document["fluid"]["Poiseuille_flag"].read(Poiseuille_flag);
    document["fluid"]["Poiseuille_bodyforce"].read(Poiseuille_bodyforce);

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
		
		
    IncomprFlowParam<T> parameters(
	        Shear_Top_Velocity,
            Shear_Top_Velocity,
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
	
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    
	if (Poiseuille_flag == 1)
    squarePoiseuilleSetup(lattice, parameters, *boundaryCondition);   //poiseuille flow boundary condition
    if (Shear_flag == 1)
	bishearSetup(lattice, parameters, *boundaryCondition); //bi-shear flow boundary condition
    
	
	
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
        spreadForce3D(lattice,wrapper);
        ////// Lattice Boltzmann iteration step.
        lattice.collideAndStream();
         //Interpolate and update solid position
        interpolateVelocity3D(lattice,wrapper);
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
