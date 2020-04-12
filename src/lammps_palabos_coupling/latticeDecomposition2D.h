/*
 * This file is part of the OpenFSI package.
 *
 * OpenFSI is free package: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2019 Huilin Ye University of Connecticut
 *
 * Author: Huilin Ye (huilin.ye@uconn.edu)
 * Note: This file is written based on lattice decomposition from LIGGGHTS parallelization
 */

#ifndef LATTICE_DECOMPOSITION2D_H
#define LATTICE_DECOMPOSITION2D_H

#include "mpi.h"
#include "lammps.h"

namespace plb{

class LatticeDecomposition2D {
public:
  LatticeDecomposition2D(plb::plint nx_, plb::plint ny_,plb::plint nz_, LAMMPS_NS::LAMMPS *lmp_);
  
  ~LatticeDecomposition2D();

  plb::SparseBlockStructure2D getBlockDistribution();
  plb::ExplicitThreadAttribution* getThreadAttribution();
private:
  plb::plint nx,ny,nz;
  LAMMPS_NS::LAMMPS &lmp;
  plb::plint npx,npy,npz;
  std::vector<plb::plint> xVal, yVal,zVal;
  plb::SparseBlockStructure2D *blockStructure;
  plb::ExplicitThreadAttribution *threadAttribution;
};

};
#include "latticeDecomposition2D.hh"

#endif /* LATTICE_DECOMPOSITION2D_H */
