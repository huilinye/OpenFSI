#ifndef LAMMPS_WRAPPER_H
#define LAMMPS_WRAPPER_H

//#include "palabos3D.h"
//#include "palabos3D.hh"

// necessary LAMMPS includes

#include "mpi.h"
#include "lammps.h"
#include "input.h"
#include "library.h"

class LammpsWrapper {
public:
  LammpsWrapper(char **argv, MPI_Comm communicator);
  //LammpsWrapper(int narg, char **argv, MPI_Comm communicator);
  void execFile(char* const fname);
  void execCommand(std::stringstream const &cmd);
  void execCommand(char* const cmd);
  //void run(plb::plint nSteps);
  //void runUpto(plb::plint nSteps);
  void run(long int nSteps);
  void runUpto(long int nSteps);
  int getNumParticles();
  void setVariable(char const *name, double value);
  void setVariable(char const *name, std::string &value);

  //private:
  LAMMPS_NS::LAMMPS *lmp;
  //LAMMPS_NS::LAMMPS **lmp;
  //int n;
};
#include "lammpsWrapper.hh"
#endif /* LAMMPS_WRAPPER_H */
