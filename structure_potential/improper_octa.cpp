/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/*
Use four node improper type to compute eight node cube volume 
and corresponding force. -----------Huilin Ye(09/11/2018) huilin.ye@uconn.edu
*/

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include "improper_octa.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

ImproperOcta::ImproperOcta(LAMMPS *lmp) : Improper(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

ImproperOcta::~ImproperOcta()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(mu);
    memory->destroy(lamda);
	memory->destroy(v0);
	memory->destroy(ele5);
	memory->destroy(ele6);
	memory->destroy(ele7);
	memory->destroy(ele8);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperOcta::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  int j1,j2,j3,j4,m,po;
  int i,j,k;
  int new_cell[4]={0};
  double xa[8], ya[8], za[8];
  tagint *tag = atom->tag;
  double **x = atom->x;
  double **f = atom->f;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  // bond information for searching
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  int co_node,se_node,s1_node,s2_node,s3_node;
  int i1_list,co_list;

  for (n = 0; n < nimproperlist; n++) {
    i1 = improperlist[n][0];
    i2 = improperlist[n][1];
    i3 = improperlist[n][2];
    i4 = improperlist[n][3];
    type = improperlist[n][4];
		
	// obtain four new nodes
	j1 = atom->map(ele5[type]);
	j2 = atom->map(ele6[type]);
	j3 = atom->map(ele7[type]);
	j4 = atom->map(ele8[type]);
	//printf("local nodes! %d %d %d %d\n",tag[i1],tag[i2],tag[i3],tag[i4]);
	//printf("new nodes! %d %d %d %d\n",tag[j1],tag[j2],tag[j3],tag[j4]);
	//printf("element information! %f %d %d %d %d\n",v0[type],ele5[type],ele6[type],ele7[type],ele8[type]);
	
	//if (m!=4) error->all(FLERR,"Cannot find nodes!");
	//printf("coordinates! %f %f %f %f\n",x[j1][0],x[j2][0],x[j3][0],x[j4][0]);
	
	xa[0] = x[i1][0]; xa[1] = x[i2][0]; xa[2] = x[i3][0]; xa[3] = x[i4][0];
	xa[4] = x[j1][0]; xa[5] = x[j2][0]; xa[6] = x[j3][0]; xa[7] = x[j4][0];
	
	ya[0] = x[i1][1]; ya[1] = x[i2][1]; ya[2] = x[i3][1]; ya[3] = x[i4][1];
	ya[4] = x[j1][1]; ya[5] = x[j2][1]; ya[6] = x[j3][1]; ya[7] = x[j4][1];
	
	za[0] = x[i1][2]; za[1] = x[i2][2]; za[2] = x[i3][2]; za[3] = x[i4][2];
	za[4] = x[j1][2]; za[5] = x[j2][2]; za[6] = x[j3][2]; za[7] = x[j4][2];

 double ddve[24][24] = {0.0};
 int ci[8],cci[8];
 double epx[8][8] = {0.0}, epy[8][8] = {0.0}, epz[8][8] = {0.0};
 double eepx[8][8] = {0.0}, eepy[8][8] = {0.0}, eepz[8][8] = {0.0};
 
  
  /* stringstream output_filename;
  output_filename << "check.dat";
  ofstream output_file;
  
  /// Open file
  
  output_file.open(output_filename.str().c_str(),ofstream::app);
  int ppc; */
  
  for (j = 0; j < 24; j++){		 
		 for (k = 0; k < 24; k++){
  ddve[j][k] = 0.0;
  }
}
  
 for (i = 0; i < 4; i++){
	 
	 for (j = 0; j < 8; j++){
		 ci[j] = nc[j][i]-1;   //minus 1 is a must, because index in c++ is from 0 rather than 1;
	 }

  	 
	 for (j = 0; j < 8; j++){		 
		 for (k = 0; k < 8; k++){
			 ddve[8+j][16+k] += kab[ci[j]][ci[k]]*xa[i];
			 ddve[16+j][8+k] += -kab[ci[j]][ci[k]]*xa[i];
			 ddve[16+j][k] += kab[ci[j]][ci[k]]*ya[i];
			 ddve[j][16+k] += -kab[ci[j]][ci[k]]*ya[i];
			 ddve[j][8+k] += kab[ci[j]][ci[k]]*za[i];
			 ddve[8+j][k] += -kab[ci[j]][ci[k]]*za[i];
		 }
	 }
	 
	 
	 
	 for (j = 0; j < 8; j++){
		 cci[j] = nc[j][i+4]-1;
	 }
		 
	 
	/*  for (j = 0; j < 8; j++){		 
		 for (k = 0; k < 8; k++){
			 eepx[j][k] += kab[cci[j]][cci[k]]*xa[i+4];
			 eepy[j][k] += kab[cci[j]][cci[k]]*ya[i+4];
			 eepz[j][k] += kab[cci[j]][cci[k]]*za[i+4];
		 }
	 } */
	 
	 for (j = 0; j < 8; j++){		 
		 for (k = 0; k < 8; k++){
			 ddve[8+j][16+k] += -kab[cci[j]][cci[k]]*xa[i+4];
			 ddve[16+j][8+k] += kab[cci[j]][cci[k]]*xa[i+4];
			 ddve[16+j][k] += -kab[cci[j]][cci[k]]*ya[i+4];
			 ddve[j][16+k] += kab[cci[j]][cci[k]]*ya[i+4];
			 ddve[j][8+k] += -kab[cci[j]][cci[k]]*za[i+4];
			 ddve[8+j][k] += kab[cci[j]][cci[k]]*za[i+4];
		 }
	 }
 }
 
 
	
 
 //printf("Finish second derivation!\n");
 double xyz[24];
 double dve[24], fpe[24];
 double ve, Jordan,pp;
 for (i = 0; i < 8; i++ ){
	 xyz[i] = xa[i];
	 xyz[i+8] = ya[i];
	 xyz[i+16] = za[i];
 }
 
 /* //===============
	 for (j = 0; j < 24; j++){		 
		 
  output_file <<xyz[j]<<"\n";
  }
   
	 
	//====================  */
 for (j = 0; j < 24; j++){
	 dve[j] = 0.0;
 }
     
 for (j = 0; j < 24; j++){
	 for (k = 0; k < 24; k++){
		 dve[j] += ddve[j][k]*xyz[k];
	 }
	 
 }
 
 
	
 ve = 0.0;
 for (j = 0; j < 24; j++){ 
		ve += dve[j]*xyz[j];
 }
 ve = ve/72.0;
 //printf("Ve! %f\n", ve);
 
 for (j = 0; j < 24; j++){ 
 dve[j] = dve[j]/24.0; 
 }
 
 for (j = 0; j < 24; j++){
	 for (k = 0; k < 24; k++){
		 ddve[j][k] = ddve[j][k]/12.0;
	 } 
 }
 //printf("Finish derivation! volume: %f\n", ve);
 Jordan = ve/v0[type];
 //printf("Jordan! %f\n", Jordan);
 if (Jordan < 1e-7) error->one(FLERR,"Negative volume!");
 
 pp = v0[type]/Jordan*(-mu[type]+lamda[type]*log(Jordan));
 //printf("pp! %f\n", pp);
 for (i = 0; i < 24; i++){
	 fpe[i] = -pp*dve[i]/v0[type];
 }
 /* //===============
	 for (j = 0; j < 24; j++){		 
		 
  output_file <<fpe[j]<<"\n";
  }
   
	 output_file.close ();
	//====================  */
  
    f[i1][0] += fpe[0]; f[i2][0] += fpe[1]; f[i3][0] += fpe[2]; f[i4][0] += fpe[3];
	f[j1][0] += fpe[4]; f[j2][0] += fpe[5]; f[j3][0] += fpe[6]; f[j4][0] += fpe[7];
	
	f[i1][1] += fpe[8]; f[i2][1] += fpe[9]; f[i3][1] += fpe[10]; f[i4][1] += fpe[11];
	f[j1][1] += fpe[12]; f[j2][1] += fpe[13]; f[j3][1] += fpe[14]; f[j4][1] += fpe[15];
	
	f[i1][2] += fpe[16]; f[i2][2] += fpe[17]; f[i3][2] += fpe[18]; f[i4][2] += fpe[19];
	f[j1][2] += fpe[20]; f[j2][2] += fpe[21]; f[j3][2] += fpe[22]; f[j4][2] += fpe[23];
           
            /* f[EleNode[j]][0] += fpe[j];
            f[EleNode[j]][1] += fpe[j+8];
            f[EleNode[j]][2] += fpe[j+16]; */
           
//printf("finish force!\n");
   
}
}

/* ---------------------------------------------------------------------- */

void ImproperOcta::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(mu,n+1,"improper:mu");
  memory->create(lamda,n+1,"improper:lamda");
  memory->create(v0,n+1,"improper:v0");
  memory->create(ele5,n+1,"improper:ele5");
  memory->create(ele6,n+1,"improper:ele6");
  memory->create(ele7,n+1,"improper:ele7");
  memory->create(ele8,n+1,"improper:ele8");
  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void ImproperOcta::coeff(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR,"Incorrect args for improper coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nimpropertypes,ilo,ihi);

  double mu_one = force->numeric(FLERR,arg[1]);
  double lamda_one = force->numeric(FLERR,arg[2]);
  double v0_one = force->numeric(FLERR,arg[3]);
  int ele5_one = force->inumeric(FLERR,arg[4]);
  int ele6_one = force->inumeric(FLERR,arg[5]);
  int ele7_one = force->inumeric(FLERR,arg[6]);
  int ele8_one = force->inumeric(FLERR,arg[7]);
 

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    mu[i] = mu_one;
    lamda[i] = lamda_one;
	v0[i] = v0_one;
	ele5[i] = ele5_one;
	ele6[i] = ele6_one;
	ele7[i] = ele7_one;
	ele8[i] = ele8_one;
    setflag[i] = 1;
    count++;
  }
  
  //prepare for the force calculation
 kab[0][0] = 0; kab[0][1] = 0; kab[0][2] = 0; kab[0][3] = 0; kab[0][4] = 0; kab[0][5] = 0; kab[0][6] = 0; kab[0][7] = 0;
 kab[1][0] = 0; kab[1][1] = 0; kab[1][2] = -1; kab[1][3] = -1.0; kab[1][4] = 1.0; kab[1][5] = 1.0; kab[1][6] = 0; kab[1][7] = 0;
 kab[2][0] = 0; kab[2][1] = 1; kab[2][2] = 0; kab[2][3] = -1.0; kab[2][4] = 0; kab[2][5] = 0; kab[2][6] = 0; kab[2][7] = 0;		  
 kab[3][0] = 0; kab[3][1] = 1; kab[3][2] = 1; kab[3][3] = 0; kab[3][4] = -1; kab[3][5] = 0; kab[3][6] = 0; kab[3][7] = -1;
 kab[4][0] = 0; kab[4][1] = -1; kab[4][2] = 0; kab[4][3] = 1; kab[4][4] = 0; kab[4][5] = -1; kab[4][6] = 0; kab[4][7] = 1;
 kab[5][0] = 0; kab[5][1] = -1; kab[5][2] = 0; kab[5][3] = 0; kab[5][4] = 1; kab[5][5] = 0; kab[5][6] = 0; kab[5][7] = 0;
 kab[6][0] = 0; kab[6][1] = 0; kab[6][2] = 0; kab[6][3] = 0; kab[6][4] = 0; kab[6][5] = 0; kab[6][6] = 0; kab[6][7] = 0;
 kab[7][0] = 0; kab[7][1] = 0; kab[7][2] = 0; kab[7][3] = 1; kab[7][4] = -1; kab[7][5] = 0; kab[7][6] = 0; kab[7][7] = 0;

 nc[0][0] = 1; nc[0][1] = 4; nc[0][2] = 3; nc[0][3] = 2; nc[0][4] = 5; nc[0][5] = 8; nc[0][6] = 7; nc[0][7] = 6;
 nc[1][0] = 2; nc[1][1] = 1; nc[1][2] = 4; nc[1][3] = 3; nc[1][4] = 6; nc[1][5] = 5; nc[1][6] = 8; nc[1][7] = 7;
 nc[2][0] = 3; nc[2][1] = 2; nc[2][2] = 1; nc[2][3] = 4; nc[2][4] = 7; nc[2][5] = 6; nc[2][6] = 5; nc[2][7] = 8;		  
 nc[3][0] = 4; nc[3][1] = 3; nc[3][2] = 2; nc[3][3] = 1; nc[3][4] = 8; nc[3][5] = 7; nc[3][6] = 6; nc[3][7] = 5;
 nc[4][0] = 5; nc[4][1] = 8; nc[4][2] = 7; nc[4][3] = 6; nc[4][4] = 1; nc[4][5] = 4; nc[4][6] = 3; nc[4][7] = 2;
 nc[5][0] = 6; nc[5][1] = 5; nc[5][2] = 8; nc[5][3] = 7; nc[5][4] = 2; nc[5][5] = 1; nc[5][6] = 4; nc[5][7] = 3;
 nc[6][0] = 7; nc[6][1] = 6; nc[6][2] = 5; nc[6][3] = 8; nc[6][4] = 3; nc[6][5] = 2; nc[6][6] = 1; nc[6][7] = 4;
 nc[7][0] = 8; nc[7][1] = 7; nc[7][2] = 6; nc[7][3] = 5; nc[7][4] = 4; nc[7][5] = 3; nc[7][6] = 2; nc[7][7] = 1;

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
  //printf("Finish initial step!\n");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void ImproperOcta::write_restart(FILE *fp)
{
  fwrite(&mu[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&lamda[1],sizeof(double),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void ImproperOcta::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&mu[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&lamda[1],sizeof(double),atom->nimpropertypes,fp);
  }
  MPI_Bcast(&mu[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&lamda[1],atom->nimpropertypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void ImproperOcta::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp,"%d %g %g\n",i,mu[i],lamda[i]);
}
