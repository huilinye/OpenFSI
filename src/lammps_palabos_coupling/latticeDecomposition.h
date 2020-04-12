/*
 * This file is part of the OpenFSI software.
 *
 * OpenFSI is free software: you can redistribute it and/or modify
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
 * This file is written based on the DEMCoupling in LIGGGHTS by the 
 * author Philippe Seil (philippe.seil@jku.at) 2014 Johannes Kepler University Linz
 *
 * Copyright 2020 University of Connecticut
 *
 * Author: Huilin Ye (huilin.ye@uconn.edu)
 */

#ifndef LATTICE_DECOMPOSITION_H
#define LATTICE_DECOMPOSITION_H

#include "mpi.h"
#include "lammps.h"

namespace plb{

class LatticeDecomposition {
public:
  LatticeDecomposition(plb::plint nx_, plb::plint ny_, plb::plint nz_, LAMMPS_NS::LAMMPS *lmp_);
  
  ~LatticeDecomposition();

  plb::SparseBlockStructure3D getBlockDistribution();
  plb::ExplicitThreadAttribution* getThreadAttribution();
private:
  plb::plint nx,ny,nz;
  LAMMPS_NS::LAMMPS &lmp;
  plb::plint npx,npy,npz;
  std::vector<plb::plint> xVal, yVal, zVal;
  plb::SparseBlockStructure3D *blockStructure;
  plb::ExplicitThreadAttribution *threadAttribution;
};

};
#include "latticeDecomposition.hh"

#endif /* LATTICE_DECOMPOSITION_H */
