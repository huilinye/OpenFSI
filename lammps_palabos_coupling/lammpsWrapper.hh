/*
 * This file is part of the Palabos_Lammps coupling program.
 *
 * Palabos_Lammps is free software: you can redistribute it and/or modify
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
 * Copyright 2018 Huilin Ye University of Connecticut
 *
 * Author: Huilin Ye (huilin.ye@uconn.edu)
 * Note: This file is written based on the IBM3D file implemented by Jifu Tan.
 */

#include "lammpsWrapper.h"
#include "mpi.h"
#include <sstream>
//#include "palabos3D.h"
//#include "palabos3D.hh"

LammpsWrapper::LammpsWrapper(char **argv, MPI_Comm communicator)
  : lmp(0)
{
  // todo: get LAMMPS to recognize command line options
  /*int argc_lmp = 1;
  char **argv_lmp = 0;
  argv_lmp = new char*[1];
  argv_lmp[0] = argv[0];*/
  //--------works for none output-----------//
  int argc_lmp = 5;
  char **argv_lmp = 0;
  argv_lmp = new char*[5];
  argv_lmp[0] = argv[0];
  argv_lmp[1]="-sc";
  argv_lmp[2]="none";
  argv_lmp[3]="-log";
  argv_lmp[4]="none";

  lmp = new LAMMPS_NS::LAMMPS(argc_lmp,argv_lmp,communicator);

  //    delete[] argv_lmp[0];
  delete[] argv_lmp;
}
void LammpsWrapper::execFile(char* const fname)
{
  lmp->input->file(fname);
}
void LammpsWrapper::execCommand(std::stringstream const &cmd)
{
  lmp->input->one(cmd.str().c_str());
}
void LammpsWrapper::execCommand(char* const cmd)
{
  lmp->input->one(cmd);
}
int LammpsWrapper::getNumParticles()
{
  return lammps_get_natoms(lmp);  
}
void LammpsWrapper::setVariable(char const *name, double value)
{
  std::stringstream cmd;
  cmd << "variable " << name << " equal " << value;
  //plb::pcout << cmd.str() << std::endl;
  execCommand(cmd);
}
void LammpsWrapper::setVariable(char const *name, std::string &value)
{
  std::stringstream cmd;
  cmd << "variable " << name << " string " << value;
  //plb::pcout << cmd.str() << std::endl;
  execCommand(cmd);
}
void LammpsWrapper::run(long int nSteps)
{
  std::stringstream cmd;
  cmd << "run " << nSteps;
  execCommand(cmd);
}
void LammpsWrapper::runUpto(long int nSteps)
{
  std::stringstream cmd;
  cmd << "run " << nSteps << " upto";
  execCommand(cmd);
}

                              
                           
                           
