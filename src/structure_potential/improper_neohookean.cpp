/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   Modified by Teng Zhang at 2/7/2017 to simulate area change
   u_a = -2k1*lnJ + k2*(lnJ)^2    k1: mu ; k2: lammda ; chi: initial area of element
   J = A/A0
------------------------------------------------------------------------- */

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include "improper_neohookean.h"
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

ImproperNeohookean::ImproperNeohookean(LAMMPS *lmp) : Improper(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

ImproperNeohookean::~ImproperNeohookean()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k1);
	memory->destroy(k2);
    memory->destroy(chi);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperNeohookean::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double eimproper,f1[3],f2[3],f3[3],f4[3];
  double domega,a, area, tk;


  eimproper = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nimproperlist; n++) {
    i1 = improperlist[n][0];
    i2 = improperlist[n][1];
    i3 = improperlist[n][2];
    i4 = improperlist[n][3];
    type = improperlist[n][4];

    // geometry of 4-body

	/*
	1,2,3,4 is the real number for calculating area
	i,j,k,l is the number in the improper angle
	
	(k)4------- 3(j)
	  !      !
	  !      !
	1(i)!------ !2(l)
	
	vb1  ---> x1-x3
	vb3  ---> x2-x4
	
	*/
    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];
	
	area = 1.0/2.0* (vb3y*vb1x - vb1y*vb3x);
//	printf("i1 %d i2 %d i3 %d i4 %d\n", i1,i2,i3,i4);
//	printf("area %12.9f\n", area);

    // force & energy
	
    domega = area/chi[type];
    // du/dA=du/dJ*(1/A0)
	
    tk = 1.0/domega*(k2[type]*log(domega) - k1[type]);

    if (eflag) eimproper = (k2[type]*log(domega) - 2.0*k1[type])*log(domega)*chi[type];

    a = -tk * 2.0;

    f1[0] = 1.0/2.0*a*vb3y;
    f1[1] = -1.0/2.0*a*vb3x;
    f1[2] = 0.0;

    f2[0] = -1.0/2.0*a*vb3y;
    f2[1] = +1.0/2.0*a*vb3x;
    f2[2] = 0.0;

    f4[0] = -1.0/2.0*a*vb1y;
    f4[1] = 1.0/2.0*a*vb1x;
    f4[2] = 0.0;

    f3[0] = 1.0/2.0*a*vb1y;
    f3[1] = -1./2.0*a*vb1x;
    f3[2] = 0.0;

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

    if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,eimproper,f1,f3,f4,
               vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperNeohookean::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(k1,n+1,"improper:k1");
  memory->create(k2,n+1,"improper:k2");
  memory->create(chi,n+1,"improper:chi");

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void ImproperNeohookean::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for improper coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nimpropertypes,ilo,ihi);

  double k1_one = force->numeric(FLERR,arg[1]);
  double k2_one = force->numeric(FLERR,arg[2]);
  double chi_one = force->numeric(FLERR,arg[3]);

  // convert chi from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k1[i] = k1_one;
	k2[i] = k2_one;
    chi[i] = chi_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void ImproperNeohookean::write_restart(FILE *fp)
{
  fwrite(&k1[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&k2[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&chi[1],sizeof(double),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void ImproperNeohookean::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k1[1],sizeof(double),atom->nimpropertypes,fp);
	fread(&k2[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&chi[1],sizeof(double),atom->nimpropertypes,fp);
  }
  MPI_Bcast(&k1[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k2[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&chi[1],atom->nimpropertypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void ImproperNeohookean::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp,"%d %g %g\n",i,k1[i],k2[i],chi[i]);
}
