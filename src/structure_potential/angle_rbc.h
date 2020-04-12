/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS

AngleStyle(rbc,AngleRbc)

#else

#ifndef LMP_ANGLE_RBC_H
#define LMP_ANGLE_RBC_H

#include <stdio.h>
#include "angle.h"

namespace LAMMPS_NS {

class AngleRbc : public Angle {
 public:
  AngleRbc(class LAMMPS *);
  virtual ~AngleRbc();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, int, int, int);

 protected:
  double *Cq;
  int    *q;
  double *ka, *kv, *kd;
  double *Atot0, *Vtot0;
  double *A0;
  
  int nmolecules;
  tagint idlo, idhi;
  
  double **cm,**cmall,**cminit;
  int *ntot_tmp,*ntot;
  double *Atot,*Vtot,*Atot_tmp,*Vtot_tmp;
  int *cmimagex,*cmimagey,*cmimagez;
  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
