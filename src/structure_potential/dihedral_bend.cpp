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

/* ----------------------------------------------------------------------
   Contributing author: Ying Li (yingli@engr.uconn.edu)
   [ based on dihedral_helix.cpp Paul Crozier (SNL) ]
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "dihedral_bend.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std; // permanently use the standard namespace
#define TOLERANCE 0.05
#define SMALL     0.001
#define SMALLER   0.00001

/* ---------------------------------------------------------------------- */

DihedralBend::DihedralBend(LAMMPS *lmp) : Dihedral(lmp) {}

/* ---------------------------------------------------------------------- */

DihedralBend::~DihedralBend()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(phi0);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBend::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double a21x,a21y,a21z,a31x,a31y,a31z;
  double a34x,a34y,a34z,a24x,a24y,a24z;
  double xix,xiy,xiz,zetax,zetay,zetaz;
  double A1,A2;
  double tc1x,tc1y,tc1z,tc2x,tc2y,tc2z;
  double tmp1x,tmp1y,tmp1z,tmp2x,tmp2y,tmp2z;
  double sign,c,s,c0,s0,phi,p,siinv;
  double beta,b11,b12,b22;
  double a32x,a32y,a32z;
  double edihedral,f1[3],f2[3],f3[3],f4[3];
  double g1[3],g2[3],g3[3],g4[3],g5[3];
  double h1[3],h2[3],h3[3],h4[3],h5[3];
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z; 
  
  double unwrap_1[3];
  double unwrap_2[3];
  double unwrap_3[3];
  double unwrap_4[3];
  int pro_pbc_pointer[3];
  int xbox[4],ybox[4],zbox[4];
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  
  edihedral = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  imageint *image = atom->image;
  int newton_bond = force->newton_bond;

  
  
  
  
  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];
    
	
	   
	   
    // 1st bond

    a21x = x[i2][0] - x[i1][0];
    a21y = x[i2][1] - x[i1][1];
    a21z = x[i2][2] - x[i1][2];
	

    // 2nd bond

    a31x = x[i3][0] - x[i1][0];
    a31y = x[i3][1] - x[i1][1];
    a31z = x[i3][2] - x[i1][2];
	
	

    // 3rd bond

    a34x = x[i3][0] - x[i4][0];
    a34y = x[i3][1] - x[i4][1];
    a34z = x[i3][2] - x[i4][2];
	
	
    
    // 4rd bond

    a24x = x[i2][0] - x[i4][0];
    a24y = x[i2][1] - x[i4][1];
    a24z = x[i2][2] - x[i4][2];
	
    
    //normal vector calculation
    // xi is the cross product of a21 and a31
    xix = a21y*a31z - a21z*a31y;
    xiy = a21z*a31x - a21x*a31z;
    xiz = a21x*a31y - a21y*a31x;
    
    //zeta is the cross product of a34 and a24
    zetax = a34y*a24z - a34z*a24y;
    zetay = a34z*a24x - a34x*a24z;
    zetaz = a34x*a24y - a34y*a24x;
    
    A1 = sqrt(xix*xix + xiy*xiy + xiz*xiz)/2.0;
    A2 = sqrt(zetax*zetax + zetay*zetay + zetaz*zetaz)/2.0;
    
    // center-of-mass calculation
       

	tc1x = (x[i1][0] + x[i2][0] +x[i3][0])/3.0;
    tc1y = (x[i1][1] + x[i2][1] +x[i3][1])/3.0; 
    tc1z = (x[i1][2] + x[i2][2] +x[i3][2])/3.0;
    
	
    tc2x = (x[i4][0] + x[i2][0] +x[i3][0])/3.0;
    tc2y = (x[i4][1] + x[i2][1] +x[i3][1])/3.0;   
    tc2z = (x[i4][2] + x[i2][2] +x[i3][2])/3.0;
	
	
    // calculate the sign
    tmp1x = xix - zetax;
    tmp1y = xiy - zetay;
    tmp1z = xiz - zetaz;
	
    tmp2x = tc1x - tc2x;
    tmp2y = tc1y - tc2y;
    tmp2z = tc1z - tc2z;
	domain->minimum_image(tmp2x,tmp2y,tmp2z);
	
    sign = tmp1x*tmp2x + tmp1y*tmp2y + tmp1z*tmp2z;
    
    // calculate the c and s
    c = (xix*zetax + xiy*zetay +xiz*zetaz)/A1/A2/4.0;
    
    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me;
      MPI_Comm_rank(world,&me);
      if (screen) {
	      char str[128];
        sprintf(str,"Dihedral problem: %d " BIGINT_FORMAT " " 
                TAGINT_FORMAT " " TAGINT_FORMAT " " 
                TAGINT_FORMAT " " TAGINT_FORMAT,
                me,update->ntimestep,
		    atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
	      error->warning(FLERR,str,0);
	      fprintf(screen,"  1st atom: %d %g %g %g\n",
		            me,x[i1][0],x[i1][1],x[i1][2]);
	      fprintf(screen,"  2nd atom: %d %g %g %g\n",
		            me,x[i2][0],x[i2][1],x[i2][2]);
	      fprintf(screen,"  3rd atom: %d %g %g %g\n",
		            me,x[i3][0],x[i3][1],x[i3][2]);
	      fprintf(screen,"  4th atom: %d %g %g %g\n",
		            me,x[i4][0],x[i4][1],x[i4][2]);
      }
    }

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    // energy
    // p = k [ 1 - cos(phi-phi0) ]
    // pd = dp/dc

    // calculate the phi and s
    phi = acos(c);
    if (sign < 0.0) phi *= -1.0;
	
	
  
    s = sin(phi);
	//s = sqrt(1.0-c*c);
    if (fabs(s) < SMALLER) {
		s = SMALLER;
	}
	//if(s < 0.0 && s>SMALLER) s=-1.0*SMALLER;
	//if(s>0.0 && s<SMALLER) s=SMALLER;
	siinv=1.0/s;

    double dphi = phi-phi0[type];
    p = k[type]*(1-cos(dphi));
    
    if (eflag) edihedral = p; 

    // force
    c0 = cos(phi0[type]);
    s0 = sin(phi0[type]);
   
    //beta = k[type]*sin(dphi)*siinv;
	beta = k[type]*(sin(phi)*c0-c*s0)*siinv;
    b11 = -beta*c/A1/A1/4.0;
    b12 =  beta/A1/A2/4.0;
    b22 = -beta*c/A2/A2/4.0;
    
    // additional bond
    a32x = x[i3][0] - x[i2][0];
    a32y = x[i3][1] - x[i2][1];
    a32z = x[i3][2] - x[i2][2];   
	
	

    // cross products!!!
    g1[0] = xiy*a32z - xiz*a32y;
    g1[1] = xiz*a32x - xix*a32z;
    g1[2] = xix*a32y - xiy*a32x;
    
    g2[0] = xiy*a31z - xiz*a31y;
    g2[1] = xiz*a31x - xix*a31z;
    g2[2] = xix*a31y - xiy*a31x;
    
    g3[0] = xiy*a21z - xiz*a21y;
    g3[1] = xiz*a21x - xix*a21z;
    g3[2] = xix*a21y - xiy*a21x;         

    g4[0] = xiy*a34z - xiz*a34y;
    g4[1] = xiz*a34x - xix*a34z;
    g4[2] = xix*a34y - xiy*a34x; 

    g5[0] = xiy*a24z - xiz*a24y;
    g5[1] = xiz*a24x - xix*a24z;
    g5[2] = xix*a24y - xiy*a24x;

    h1[0] = zetay*a32z - zetaz*a32y;
    h1[1] = zetaz*a32x - zetax*a32z;
    h1[2] = zetax*a32y - zetay*a32x;
    
    h2[0] = zetay*a31z - zetaz*a31y;
    h2[1] = zetaz*a31x - zetax*a31z;
    h2[2] = zetax*a31y - zetay*a31x;
    
    h3[0] = zetay*a21z - zetaz*a21y;
    h3[1] = zetaz*a21x - zetax*a21z;
    h3[2] = zetax*a21y - zetay*a21x;         

    h4[0] = zetay*a34z - zetaz*a34y;
    h4[1] = zetaz*a34x - zetax*a34z;
    h4[2] = zetax*a34y - zetay*a34x; 

    h5[0] = zetay*a24z - zetaz*a24y;
    h5[1] = zetaz*a24x - zetax*a24z;
    h5[2] = zetax*a24y - zetay*a24x;
    
    // calculate the force vector
    f1[0] = b11*g1[0] + b12*h1[0];
    f1[1] = b11*g1[1] + b12*h1[1];
    f1[2] = b11*g1[2] + b12*h1[2];

    f2[0] = -b11*g2[0] + b12*(g4[0]-h2[0]) + b22*h4[0];
    f2[1] = -b11*g2[1] + b12*(g4[1]-h2[1]) + b22*h4[1];
    f2[2] = -b11*g2[2] + b12*(g4[2]-h2[2]) + b22*h4[2];

    f3[0] = b11*g3[0] + b12*(-g5[0]+h3[0]) - b22*h5[0];
    f3[1] = b11*g3[1] + b12*(-g5[1]+h3[1]) - b22*h5[1];
    f3[2] = b11*g3[2] + b12*(-g5[2]+h3[2]) - b22*h5[2];

    f4[0] = -b12*g1[0] - b22*h1[0];
    f4[1] = -b12*g1[1] - b22*h1[1];
    f4[2] = -b12*g1[2] - b22*h1[2];

    // apply force to each of 4 atoms

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

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }
    
    // 1st bond

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];

    // 3rd bond

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];

    if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,edihedral,f1,f3,f4,
	       vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
  }
  /* for (n = 0; n < atom->nlocal; n++) {
//   check the dihedral_bend force in the LAMMPS
      if(update->ntimestep % 1000 == 0){
      printf("     bend_F_x          bend_F_y          bend_F_z\n");
      printf("%16.12f %16.12f %16.12f \n",f[n][0],f[n][1],f[n][2]);	
	}  
	 } */
}

/* ---------------------------------------------------------------------- */

void DihedralBend::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(k,n+1,"dihedral:k");
  memory->create(phi0,n+1,"dihedral:phi0");

  memory->create(setflag,n+1,"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void DihedralBend::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for dihedral coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->ndihedraltypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double phi0_one= force->numeric(FLERR,arg[2]);

  // require k >= 0
  if (k_one < 0.0)
    error->all(FLERR,"Incorrect coefficient arg for dihedral coefficients");
                       
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    phi0[i] = phi0_one*MY_PI/180.0;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void DihedralBend::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&phi0[1],sizeof(double),atom->ndihedraltypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void DihedralBend::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&phi0[1],sizeof(double),atom->ndihedraltypes,fp);
  }
  MPI_Bcast(&k[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi0[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
 
  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;
}

