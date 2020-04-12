/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   Written by Ying Li (yingli@engr.uconn.edu)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include "bond_wlc.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams
using namespace LAMMPS_NS;
using namespace std; // permanently use the standard namespace
/* ---------------------------------------------------------------------- */

BondWlc::BondWlc(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondWlc::~BondWlc()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(kT);
    memory->destroy(rnorm0);
    memory->destroy(rmax);
    memory->destroy(mu0);
    memory->destroy(m);
    memory->destroy(gammaC);
    memory->destroy(gammaT);
  }
}

/* ---------------------------------------------------------------------- */

void BondWlc::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,rnorm,rfactor1,rfactor11,rfactor2,rfactor3;
  double p,kp,r0;
  double a11,a12,a21,a22;

  double unwrap_1[3];
  double unwrap_2[3];
  
  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  imageint *image = atom->image;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
	

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);       // length of the bond
    rnorm = r/rmax[type];  // x
	
	//output bond length
	/* if(update->ntimestep == 0){
  stringstream output_filename;
  output_filename << "Bond_length.txt";
  ofstream output_file;
  
  output_file.open(output_filename.str().c_str(),ofstream::app);
   
  output_file << "Bond   " << n << " " << "Length  " << r << "" <<"bond_node"<<i1<<" "<<i2<<"\n";
   } */
  
  
    rfactor1 = 1.0/4.0/(1.0-rnorm)/(1.0-rnorm)-1.0/4.0+rnorm;
	rfactor11 = 1.0/4.0/(1.0-rnorm0[type])/(1.0-rnorm0[type])-1.0/4.0+rnorm0[type];
    rfactor2 = 2.0*rnorm*rnorm-rnorm-1.0+1.0/(1.0-rnorm);
    rfactor3 = rnorm0[type]/2.0/pow(1.0-rnorm0[type],3.0)-1.0/4.0/(1.0-rnorm0[type])/(1.0-rnorm0[type])+1.0/4.0;
    
    // use force equlibrium and macro property equation to compute kp and p
	r0 = rnorm0[type]*rmax[type];
	
	// set zero force in initial state
	/* a11=sqrt(3.0)*kT[type]/4.0/r0*rfactor3;
	a12=sqrt(3.0)*(m[type]+1)/4.0/pow(r0,m[type]+1.0);
	a21=-kT[type]*rfactor11;
	a22=1.0/pow(r0,m[type]);
	kp = -a21*mu0[type]/(a11*a22-a12*a21);
	p  = (a11*a22-a12*a21)/mu0[type]/a22; */
	
      
    //p = 0.001118/r0;
	//p=0.0204*r0;
	//p = 0.0006;     //adhesion for WBC
	p = 0.008;
    kp = (mu0[type]-sqrt(3.0)*kT[type]*rfactor3/4.0/p/r0)*4.0*pow(r0,m[type]+1.0)/sqrt(3.0)/(m[type]+1.0); 

    // force & energy

    if (r > 0.0) fbond = -(kT[type]*rfactor1/4.0/p-kp/pow(r,m[type]))/r;
    else fbond = 0.0;

    if (eflag) {
       if (m[type]==1.0) ebond = kT[type]*rmax[type]*rfactor2/4.0/p-kp*log(r);
       else ebond = kT[type]*rmax[type]*rfactor2/4.0/p+kp/(m[type]-1.0)/pow(r,m[type]-1.0); 
    }
    
    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondWlc::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(kT,n+1,"bond:kT");
  memory->create(rnorm0,n+1,"bond:rnorm0");
  memory->create(rmax,n+1,"bond:rmax");
  memory->create(mu0,n+1,"bond:mu0");
  memory->create(m,n+1,"bond:m");
  memory->create(gammaC,n+1,"bond:gammaC");
  memory->create(gammaT,n+1,"bond:gammaT");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondWlc::coeff(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double kT_one = force->numeric(FLERR,arg[1]);
  double rnorm0_one = force->numeric(FLERR,arg[2]);
  double rmax_one = force->numeric(FLERR,arg[3]);
  double mu0_one = force->numeric(FLERR,arg[4]);
  double m_one = force->numeric(FLERR,arg[5]);
  double gammaC_one = force->numeric(FLERR,arg[6]);
  double gammaT_one = force->numeric(FLERR,arg[7]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    kT[i] = kT_one;
    rnorm0[i] = rnorm0_one;
    rmax[i] = rmax_one;
    mu0[i] = mu0_one;
    m[i] = m_one;
    gammaC[i] = gammaC_one;
    gammaT[i] = gammaT_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondWlc::equilibrium_distance(int i)
{
  return rnorm0[i]*rmax[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondWlc::write_restart(FILE *fp)
{
  fwrite(&kT[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&rnorm0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&rmax[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&mu0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&m[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&gammaC[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&gammaT[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondWlc::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&kT[1],sizeof(double),atom->nbondtypes,fp);
    fread(&rnorm0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&rmax[1],sizeof(double),atom->nbondtypes,fp);
    fread(&mu0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&m[1],sizeof(double),atom->nbondtypes,fp);
    fread(&gammaC[1],sizeof(double),atom->nbondtypes,fp);
    fread(&gammaT[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&kT[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&rnorm0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&rmax[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&mu0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&m[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gammaC[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gammaT[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondWlc::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g\n",i,kT[i],rnorm0[i],rmax[i],mu0[i],m[i],gammaC[i],gammaT[i]);
}

/* ---------------------------------------------------------------------- */

double BondWlc::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  double r = sqrt(rsq);
  double rnorm = r/rmax[type];
  double rfactor1 = 1.0/4.0/(1.0-rnorm)/(1.0-rnorm)-1.0/4.0+rnorm;
  double rfactor2 = 2.0*rnorm*rnorm-rnorm-1.0+1.0/(1.0-rnorm);
  double rfactor3 = rnorm0[type]/2.0/pow(1.0-rnorm0[type],3.0)-1.0/4.0/(1.0-rnorm0[type])/(1.0-rnorm0[type])+1.0/4.0;
  
  // keep the p*r0 constant & use mu0 to calculate the kp
  double r0 = rnorm0[type]*rmax[type];
  double p = 0.001118/r0;
  double kp = (mu0[type]-sqrt(3.0)*kT[type]*rfactor3/4.0/p/r0)*4.0*pow(r0,m[type]+1.0)/sqrt(3.0)/(m[type]+1.0); 
    
  fforce = 0;
  if (r > 0.0) fforce = -(kT[type]*rfactor1/4.0/p-kp/pow(r,m[type]))/r;
  
  if (m[type]==1.0) return kT[type]*rmax[type]*rfactor2/4.0/p-kp*log(r);
  else return kT[type]*rmax[type]*rfactor2/4.0/p+kp/(m[type]-1.0)/pow(r,m[type]-1.0);
}
