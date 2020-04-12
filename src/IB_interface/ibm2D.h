#ifndef IBM_LBM_2D_H
#define IBM_LBM_2D_H

#include "lammpsWrapper.h"

namespace plb {
  						   
  template<typename T, template<typename U> class Descriptor>
  void interpolateVelocity2D_fix(MultiBlockLattice2D<T,Descriptor> &lattice,
                           LammpsWrapper &wrapper, plint StaticAtomType);
  

  template<typename T, template<typename U> class Descriptor>
  void spreadForce2D_fix(MultiBlockLattice2D<T,Descriptor> &lattice,
                           LammpsWrapper &wrapper, plint StaticAtomType);
						   
  template<typename T, template<typename U> class Descriptor>
  void forceCoupling2D(MultiBlockLattice2D<T,Descriptor> &lattice,
                           LammpsWrapper &wrapper);
						   
  
}; /* namespace plb */

#include "ibm2D.hh"

#endif 
