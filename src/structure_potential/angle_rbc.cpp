/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   Designed for conserving the local area, total area and total volume
   Contributing author: Ying Li (yingli@engr.uconn.edu)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "angle_rbc.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "group.h"
#include "modify.h"

#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std; // permanently use the standard namespace
#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleRbc::AngleRbc(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleRbc::~AngleRbc()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(Cq);
    memory->destroy(q);
    memory->destroy(ka);
    memory->destroy(Atot0);
    memory->destroy(kv);
    memory->destroy(Vtot0);
    memory->destroy(kd);
    memory->destroy(A0);
  }
  /* memory->destroy(cm);
  memory->destroy(cminit);
  memory->destroy(cmall);
  memory->destroy(ntot_tmp);
  memory->destroy(ntot);
  memory->destroy(Atot_tmp);
  memory->destroy(Atot);
  memory->destroy(Vtot_tmp);
  memory->destroy(Vtot);
  memory->destroy(cmimagex);
  memory->destroy(cmimagey);
  memory->destroy(cmimagez); */
  
}

/* ---------------------------------------------------------------------- */

void AngleRbc::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f2[3],f3[3];
  double g1[3],g2[3],g3[3];
  double h1[3],h2[3],h3[3];
  // added for calculating the global area & volume and local area
  double a21x,a21y,a21z,a32x,a32y,a32z,a13x,a13y,a13z;
  double zetax,zetay,zetaz,Ak,Vk;
  double tcx,tcy,tcz;

  double alpha_h, alpha_a;
  double beta_a, beta_v;
  double unwrap[3];
  double cmdx,cmdy,cmdz;
  double angle_point, angle_point_temp;

  

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  // compute parameters for molecule
  tagint imol;
  
  nmolecules = molecules_in_group(idlo,idhi);
  
  memory->create(cm,nmolecules,3,"angle/molecule:cm");
  memory->create(cminit,nmolecules,3,"angle/molecule:cminit");
  memory->create(cmall,nmolecules,3,"angle/molecule:cmall");
  memory->create(ntot_tmp,nmolecules,"angle/molecule:ntot_tmp");
  memory->create(ntot,nmolecules,"angle/molecule:ntot");
  memory->create(Atot_tmp,nmolecules,"angle/molecule:Atot_tmp");
  memory->create(Atot,nmolecules,"angle/molecule:Atot");
  memory->create(Vtot_tmp,nmolecules,"angle/molecule:Vtot_tmp");
  memory->create(Vtot,nmolecules,"angle/molecule:Vtot");
  memory->create(cmimagex,nmolecules,"angle/molecule:cmimagex");
  memory->create(cmimagey,nmolecules,"angle/molecule:cmimagey");
  memory->create(cmimagez,nmolecules,"angle/molecule:cmimagez");
  
  tagint *molecule = atom->molecule;
  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int newton_bond = force->newton_bond;
  
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  // need to calculate the center-of-mass for each molecule
  for (int i = 0; i < nmolecules; i++){
	  cm[i][0] = cm[i][1] = cm[i][2] = 0.0;
      ntot_tmp[i] = 0;
	  if(update->ntimestep == 0){
		  cminit[i][0] = 0.0;
          cminit[i][1] = 0.0;
          cminit[i][2] = 0.0;
		  cmimagex[i] = 0;
	      cmimagey[i] = 0;
	      cmimagez[i] = 0;
	  }
  }
  
  for (int i = 0; i < nlocal; i++){
	  if (mask[i]) {
	  imol = molecule[i];
	  
	  if (molmap) imol = molmap[imol-idlo];
      else  imol--; 
	  domain->unmap(x[i],image[i],unwrap);
	  cm[imol][0] += unwrap[0];
	  cm[imol][1] += unwrap[1];
      cm[imol][2] += unwrap[2];
      ntot_tmp[imol] += 1;  
  }
  }
  
  
  
  // communicate the glocal values
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(cm[0],cmall[0],3*nmolecules,
                MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(ntot_tmp,ntot,nmolecules,
                MPI_INT,MPI_SUM,world);	
  
  // take the average value for center-of-mass
  for (int i = 0; i < nmolecules; i++) {
  cmall[i][0] /= ntot[i];
  cmall[i][1] /= ntot[i];
  cmall[i][2] /= ntot[i];
  }
  
  // output the center of mass 
  if(update->ntimestep %100 == 0 && comm->me == 0 ){
  stringstream output_filename;
  output_filename << "center_of_mass.dat";
  ofstream output_file;
  
  /// Open file
  
  output_file.open(output_filename.str().c_str(),ofstream::app);
  for (int i = 0; i < nmolecules; i++) {
  output_file <<i<<" "<<"molecule"<<" "<< "center  " << cmall[i][0] << " " << 
  cmall[i][1] << " "<<cmall[i][2] <<" "<<"totnumber"<<ntot[i]<< "\n";
  }
  }
  
  //compute the initial center of mass
  /* if(update->ntimestep == 0){
	  
	 for (int i = 0; i < nmolecules; i++) {
    cminit[i][0] = cmall[i][0];
    cminit[i][1] = cmall[i][1];
    cminit[i][2] = cmall[i][2];
  } 
	  
  } */
  //wrap the value of the center of mass
  // using the wrap coordinates of atoms as initial center of mass
  for (int i = 0; i < nlocal; i++){
	  if (mask[i]) {
	  imol = molecule[i];
	  
	  if (molmap) imol = molmap[imol-idlo];
      else imol--;
	  
	  cminit[imol][0] = x[i][0];
      cminit[imol][1] = x[i][1];
      cminit[imol][2] = x[i][2];
	  
  }
  }
  for (int i = 0; i < nmolecules; i++) {
	  
    cmdx = cmall[i][0] - cminit[i][0];
    cmdy = cmall[i][1] - cminit[i][1];
    cmdz = cmall[i][2] - cminit[i][2];
	
	cmimagex[i] = (int)floor(fabs(cmdx)/xprd);
	cmimagey[i] = (int)floor(fabs(cmdy)/yprd);
	cmimagez[i] = (int)floor(fabs(cmdz)/zprd);
	
	if(cmdx < 0.0) cmimagex[i] = 0-cmimagex[i];
	if(cmdy < 0.0) cmimagey[i] = 0-cmimagey[i];
	if(cmdz < 0.0) cmimagez[i] = 0-cmimagez[i];
	
	cmall[i][0] -= 1.0*cmimagex[i]*xprd;
	cmall[i][1] -= 1.0*cmimagey[i]*yprd;
	cmall[i][2] -= 1.0*cmimagez[i]*zprd;
  } 
  
  
  
  
  
  
  
  // need to calculate the global area and volume

  for (int i = 0; i < nmolecules; i++)
	  Atot_tmp[i] = Vtot_tmp[i] = 0.0;
       
	  
  for (n = 0; n < nanglelist; n++) 
	{
      
	  
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
	
	if (mask[i1]) {
	  imol = molecule[i1];
	  
	  if (molmap) imol = molmap[imol-idlo];
      else imol--;
	  }
		
	
    // 1st vector
    a21x = x[i2][0] - x[i1][0];
    a21y = x[i2][1] - x[i1][1];
    a21z = x[i2][2] - x[i1][2];
	
    // 2nd vector

    a32x = x[i3][0] - x[i2][0];
    a32y = x[i3][1] - x[i2][1];
    a32z = x[i3][2] - x[i2][2];
    
    // 3nd vector

    a13x = x[i1][0] - x[i3][0];
    a13y = x[i1][1] - x[i3][1];
    a13z = x[i1][2] - x[i3][2];
	
	
    // normal out-of-plane vector
    // zeta is cross product of a21 and a31 
    
    zetax = -a21y*a13z + a21z*a13y; 
    zetay = -a21z*a13x + a21x*a13z; 
    zetaz = -a21x*a13y + a21y*a13x;
    
    // calculate the local area
    Ak = sqrt(zetax*zetax + zetay*zetay + zetaz*zetaz)/2.0;
    
	
  
    // calculate the glocal area
    Atot_tmp[imol] += Ak;  

    // find the center-of-mass
	tcx = (x[i1][0] + x[i2][0] + x[i3][0])/3.0- cmall[imol][0];
    tcy = (x[i1][1] + x[i2][1] + x[i3][1])/3.0- cmall[imol][1];
    tcz = (x[i1][2] + x[i2][2] + x[i3][2])/3.0- cmall[imol][2];
	
	
	domain->minimum_image(tcx,tcy,tcz);
	
	angle_point_temp = zetax*tcx+zetay*tcy+zetaz*tcz;
	if(angle_point_temp <= 0.0) angle_point = -1.0;
	else angle_point = 1.0;
	
	
	
    // calculate the local volume
    Vk = angle_point*(zetax*tcx + zetay*tcy + zetaz*tcz)/6.0;
    
    // calculate the glocal volume
    Vtot_tmp[imol] += Vk; 
  }
  
  // communicate the glocal area & volume
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(Atot_tmp,Atot,nmolecules,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(Vtot_tmp,Vtot,nmolecules,MPI_DOUBLE,MPI_SUM,world);
  
  // output the global volume and area 
  /* if(update->ntimestep %100 == 0 && comm->me == 0){
  stringstream output_filename;
  output_filename << "Total_Volume_Area.dat";
  ofstream output_file;
  
  /// Open file
  
  output_file.open(output_filename.str().c_str(),ofstream::app);
  for (int i = 0; i < nmolecules; i++) {
  output_file <<i<<" "<<"molecule"<<" "<< "Volume  " << Vtot[i] << " " << "Area  " << Atot[i] << "\n";
  }
  } */
  // assign the local nodal forces
  
  for (n = 0; n < nanglelist; n++) 
  {
      
	  
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
	
	if (mask[i1] ) {
	imol = molecule[i1];
	
	  if (molmap) imol = molmap[imol-idlo];
      else imol--;
	}
		
    // 1st vector

    a21x = x[i2][0] - x[i1][0];
    a21y = x[i2][1] - x[i1][1];
    a21z = x[i2][2] - x[i1][2];
	

    // 2nd vector

    a32x = x[i3][0] - x[i2][0];
    a32y = x[i3][1] - x[i2][1];
    a32z = x[i3][2] - x[i2][2];
    
	
    // 3nd vector

    a13x = x[i1][0] - x[i3][0];
    a13y = x[i1][1] - x[i3][1];
    a13z = x[i1][2] - x[i3][2];
	
    // normal out-of-plane vector
    // zeta is cross product of a21 and a31 
    
    zetax = -a21y*a13z + a21z*a13y; 
    zetay = -a21z*a13x + a21x*a13z; 
    zetaz = -a21x*a13y + a21y*a13x;
    
    // calculate the local area
    Ak = sqrt(zetax*zetax + zetay*zetay + zetaz*zetaz)/2.0;
    
    // find the center-of-mass
	tcx = (x[i1][0] + x[i2][0] + x[i3][0])/3.0- cmall[imol][0];
    tcy = (x[i1][1] + x[i2][1] + x[i3][1])/3.0- cmall[imol][1];
    tcz = (x[i1][2] + x[i2][2] + x[i3][2])/3.0- cmall[imol][2];
	
   
	domain->minimum_image(tcx,tcy,tcz);
	
	angle_point_temp = zetax*tcx+zetay*tcy+zetaz*tcz;
	if(angle_point_temp <= 0.0) angle_point = -1.0;
	else angle_point = 1.0;
    
    // calculate the components along x,y,z directions
    // h1,2,3 are cross products of zeta and a32, a13, a21, respectively
    h1[0] = zetay*a32z - zetaz*a32y;
    h1[1] = zetaz*a32x - zetax*a32z;
    h1[2] = zetax*a32y - zetay*a32x;
    h2[0] = zetay*a13z - zetaz*a13y;
    h2[1] = zetaz*a13x - zetax*a13z;
    h2[2] = zetax*a13y - zetay*a13x;
    h3[0] = zetay*a21z - zetaz*a21y;
    h3[1] = zetaz*a21x - zetax*a21z;
    h3[2] = zetax*a21y - zetay*a21x;
    
    // g1,2,3 are cross products of tc and a32, a13, a21, respectively
    g1[0] = tcy*a32z - tcz*a32y + zetax/3.0;
    g1[1] = tcz*a32x - tcx*a32z + zetay/3.0;
    g1[2] = tcx*a32y - tcy*a32x + zetaz/3.0;
    g2[0] = tcy*a13z - tcz*a13y + zetax/3.0;
    g2[1] = tcz*a13x - tcx*a13z + zetay/3.0;
    g2[2] = tcx*a13y - tcy*a13x + zetaz/3.0;
    g3[0] = tcy*a21z - tcz*a21y + zetax/3.0;
    g3[1] = tcz*a21x - tcx*a21z + zetay/3.0;
    g3[2] = tcx*a21y - tcy*a21x + zetaz/3.0;
	
	// g1,2,3 are cross products of tc and a32, a13, a21, respectively
    /* g1[0] = angle_point*(tcy*a32z - tcz*a32y) + zetax/3.0;
    g1[1] = angle_point*(tcz*a32x - tcx*a32z) + zetay/3.0;
    g1[2] = angle_point*(tcx*a32y - tcy*a32x) + zetaz/3.0;
    g2[0] = angle_point*(tcy*a13z - tcz*a13y) + zetax/3.0;
    g2[1] = angle_point*(tcz*a13x - tcx*a13z) + zetay/3.0;
    g2[2] = angle_point*(tcx*a13y - tcy*a13x) + zetaz/3.0;
    g3[0] = angle_point*(tcy*a21z - tcz*a21y) + zetax/3.0;
    g3[1] = angle_point*(tcz*a21x - tcx*a21z) + zetay/3.0;
    g3[2] = angle_point*(tcx*a21y - tcy*a21x) + zetaz/3.0; */
    
    //calculate all the coefficients
    alpha_h = q[type]*Cq[type]/4.0/pow(Ak,q[type]+2.0); 
    alpha_a = -kd[type]*(Ak-A0[type])/4.0/Ak/A0[type];
    beta_a = -ka[type]*(Atot[imol]-Atot0[type])/4.0/Ak/Atot0[type];
    beta_v = -kv[type]*(Vtot[imol]-Vtot0[type])/6.0/Vtot0[type];
    
    // combine the force together
    f1[0] = h1[0]*(alpha_h+alpha_a+beta_a)+g1[0]*beta_v;
    f1[1] = h1[1]*(alpha_h+alpha_a+beta_a)+g1[1]*beta_v;
    f1[2] = h1[2]*(alpha_h+alpha_a+beta_a)+g1[2]*beta_v;
    f2[0] = h2[0]*(alpha_h+alpha_a+beta_a)+g2[0]*beta_v;
    f2[1] = h2[1]*(alpha_h+alpha_a+beta_a)+g2[1]*beta_v;
    f2[2] = h2[2]*(alpha_h+alpha_a+beta_a)+g2[2]*beta_v; 
    f3[0] = h3[0]*(alpha_h+alpha_a+beta_a)+g3[0]*beta_v;
    f3[1] = h3[1]*(alpha_h+alpha_a+beta_a)+g3[1]*beta_v;
    f3[2] = h3[2]*(alpha_h+alpha_a+beta_a)+g3[2]*beta_v;   

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }
    
    // calculate the energy
    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];
    
    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];
    
    if (eflag) eangle = Cq[type]/pow(Ak,q[type])+kd[type]*(Ak-A0[type])*(Ak-A0[type])/2.0/A0[type];

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
  /* for (n = 0; n < atom->nlocal; n++) {
//   check the angle force in the LAMMPS
      if(update->ntimestep % 1000 == 0){
      printf("     angle_F_x          angle_F_y          angle_F_z\n");
      printf("%16.12f %16.12f %16.12f \n",f[n][0],f[n][1],f[n][2]);	
	}  
	 } */
  memory->destroy(cm);
  memory->destroy(cminit);
  memory->destroy(cmall);
  memory->destroy(ntot_tmp);
  memory->destroy(ntot);
  memory->destroy(Atot_tmp);
  memory->destroy(Atot);
  memory->destroy(Vtot_tmp);
  memory->destroy(Vtot);
  memory->destroy(cmimagex);
  memory->destroy(cmimagey);
  memory->destroy(cmimagez);
  
}

/* ---------------------------------------------------------------------- */

void AngleRbc::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(Cq,n+1,"angle:Cq");
  memory->create(q,n+1,"angle:q");
  memory->create(ka,n+1,"angle:kd");
  memory->create(Atot0,n+1,"angle:Atot0");
  memory->create(kv,n+1,"angle:kv");
  memory->create(Vtot0,n+1,"angle:Vtot0");
  memory->create(kd,n+1,"angle:kd");
  memory->create(A0,n+1,"angle:A0");  

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleRbc::coeff(int narg, char **arg)
{
  if (narg != 9) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nangletypes,ilo,ihi);

  double Cq_one = force->numeric(FLERR,arg[1]);
  double q_one = force->numeric(FLERR,arg[2]);
  double ka_one = force->numeric(FLERR,arg[3]);
  double Atot0_one = force->numeric(FLERR,arg[4]);
  double kv_one = force->numeric(FLERR,arg[5]);
  double Vtot0_one = force->numeric(FLERR,arg[6]);
  double kd_one = force->numeric(FLERR,arg[7]);
  double A0_one = force->numeric(FLERR,arg[8]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    Cq[i] = Cq_one;
    q[i] = q_one;
    ka[i] = ka_one;
    Atot0[i] = Atot0_one;
    kv[i] = kv_one;
    Vtot0[i] = Vtot0_one;
    kd[i] = kd_one;
    A0[i] = A0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleRbc::equilibrium_angle(int i)
{
  return MY_PI/3;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleRbc::write_restart(FILE *fp)
{
  fwrite(&Cq[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&q[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ka[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&Atot0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kv[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&Vtot0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kd[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&A0[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleRbc::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&Cq[1],sizeof(double),atom->nangletypes,fp);
    fread(&q[1],sizeof(double),atom->nangletypes,fp);
    fread(&ka[1],sizeof(double),atom->nangletypes,fp);
    fread(&Atot0[1],sizeof(double),atom->nangletypes,fp);
    fread(&kv[1],sizeof(double),atom->nangletypes,fp);
    fread(&Vtot0[1],sizeof(double),atom->nangletypes,fp);
    fread(&kd[1],sizeof(double),atom->nangletypes,fp);
    fread(&A0[1],sizeof(double),atom->nangletypes,fp);
  }
  MPI_Bcast(&Cq[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&q[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ka[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Atot0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kv[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Vtot0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kd[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&A0[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleRbc::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g %g\n",i,Cq[i],q[i],ka[i],Atot0[i],kv[i],Vtot0[i],kd[i],A0[i]);
}

/* ---------------------------------------------------------------------- */

double AngleRbc::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double a21x = x[i2][0] - x[i1][0];
  double a21y = x[i2][1] - x[i1][1];
  double a21z = x[i2][2] - x[i1][2];
  domain->minimum_image(a21x,a21y,a21z);
  
  double a32x = x[i3][0] - x[i2][0];
  double a32y = x[i3][1] - x[i2][1];
  double a32z = x[i3][2] - x[i2][2];
  domain->minimum_image(a32x,a32y,a32z);

  double a13x = x[i1][0] - x[i3][0];
  double a13y = x[i1][1] - x[i3][1];
  double a13z = x[i1][2] - x[i3][2];
  domain->minimum_image(a13x,a13y,a13z);
  
  double zetax = -a21y*a13z + a21z*a13y; 
  double zetay = -a21z*a13x + a21x*a13z; 
  double zetaz = -a21x*a13y + a21y*a13x;
  
  double Ak = sqrt(zetax*zetax + zetay*zetay + zetaz*zetaz)/2.0;

  return Cq[type]/pow(Ak,q[type])+kd[type]*(Ak-A0[type])*(Ak-A0[type])/2.0/A0[type];
}

