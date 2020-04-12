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
 * Note: This file is written based on Jifu Tan(https://github.com/TJFord/palabos-lammps).
 */

#ifndef IBM_LBM_2D_HH
#define IBM_LBM_2D_HH

#include "atom.h"
#include "modify.h"
#include "group.h"
#include "molecule.h"
#include "fix.h"
#include "fix_fcm.h"
#include "update.h"
#include <algorithm>
#include <map>
#include "pointers.h"

namespace plb {
	
	template<typename T>
  void weight(T r, std::vector<T> & w){
      T q = sqrt(1 + 4*r*(1-r));
      w[0] = (3 - 2*r - q)/8.0;
      w[1] = (3 - 2*r + q)/8.0;
      w[2] = (1 + 2*r + q)/8.0;
      w[3] = (1 + 2*r - q)/8.0;
  }

  //******************************
  //  interpolation velocity
  //******************************
  template<typename T, template<typename U> class Descriptor>
  class Interpolation2D_fix: public BoxProcessingFunctional2D_L<T,Descriptor>{
    public:
      Interpolation2D_fix(LammpsWrapper &wrapper_, plint StaticAtomType_):wrapper(wrapper_), StaticAtomType(StaticAtomType_){
        
        dt = wrapper.lmp->update->dt;
      }
      virtual void process(Box2D domain, BlockLattice2D<T,Descriptor> &lattice){
        Dot2D offset = lattice.getLocation();
        TensorField2D<T,Descriptor<T>::d> velocity(lattice.getNx(),lattice.getNy());
        
        //plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        //T rx,ry,rz,wgt,rho;
        Array<T,2> us(0.,0.);
        Array<T,2> uf;
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        int *mask = wrapper.lmp->atom->mask;
		int *type = wrapper.lmp->atom->type;
        plint nlocal = wrapper.lmp->atom->nlocal;
        //std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
		
		plint ix,iy;
        plint ixp,iyp;
        T dx1,dy1;
        plint isten,ii,jj;
        T r,rsq,weightx,weighty;
        T Ffp[16];
        
		
        plint env = 2;
        for(ix=domain.x0-env;ix<=domain.x1+env;ix++)
        for(iy=domain.y0-env;iy<=domain.y1+env;iy++){
        //for(iz=domain.z0-env;iz<=domain.z1+env;iz++)
          lattice.get(ix,iy).computeVelocity(velocity.get(ix,iy));
        }
		
        for (plint iS=0; iS<nlocal; iS++){
			
          //if(mask[iS] && type[iS] < MAX_TYPE){
			if(mask[iS] && type[iS] != StaticAtomType){  
			  us[0] = us[1] = 0.0;
			  ix = (int)ceil(x[iS][0]-offset.x);
			  iy = (int)ceil(x[iS][1]-offset.y);
			  //iz = (int)ceil(x[iS][2]-offset.z);
			  
			  dx1 = x[iS][0] - offset.x -ix + 1.0;
			  dy1 = x[iS][1] - offset.y -iy + 1.0;

             isten = 0;			 
    for (ii=-1;ii<3;ii++ ){
				rsq = (-dx1+ii)*(-dx1+ii);
    if(rsq>=4)
      weightx=0.0;
    else{
      r=sqrt(rsq);
      if(rsq>1){
	weightx=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
      } else{
	weightx=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
      }
    }
    for(jj=-1; jj<3; jj++){
      rsq=(-dy1+jj)*(-dy1+jj);
      if(rsq>=4)
	weighty=0.0;
      else{
	r=sqrt(rsq);
	if(rsq>1){
	  weighty=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
	} else{
	  weighty=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
	}
      }
      
	
	ixp = ix+ii;
	iyp = iy+jj;
	
				
				if(ixp<domain.x0-env || ixp>domain.x1+env) continue;//ixp = domain.x1+2;
	            if(iyp<domain.y0-env || iyp>domain.y1+env) continue;
	            
		  
				Ffp[isten] = weightx*weighty;
                 
                 int shift = 0;
                  uf = velocity.get(ixp+shift,iyp+shift);
				  
				  Cell<T,Descriptor>& cell  = lattice.get(ixp+shift,iyp+shift);
                  T *bodyforce=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                  
                  us[0] += Ffp[isten]*(uf[0] + 0.5*bodyforce[0]/lattice.get(ixp,iyp).computeDensity());  //update velocity
                  us[1] += Ffp[isten]*(uf[1] + 0.5*bodyforce[1]/lattice.get(ixp,iyp).computeDensity());
                  
				  
				  isten++;
                }
	
			}
			
            v[iS][0]=us[0];
            v[iS][1]=us[1];
            //v[iS][2]=us[2];
			
            //Euler method to update position
            x[iS][0] += v[iS][0]*dt;
            x[iS][1] += v[iS][1]*dt;
            
          
        
      }
		}
	  }
      virtual Interpolation2D_fix<T,Descriptor> * clone() const{
        return new Interpolation2D_fix(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::nothing; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      plint StaticAtomType;
      T dt;
  };


  template<typename T, template<typename U> class Descriptor>
  void interpolateVelocity2D_fix(MultiBlockLattice2D<T,Descriptor> &lattice, LammpsWrapper &wrapper, plint StaticAtomType)
  {
    //plint envelopeWidth = 3;
    //applyProcessingFunctional(new Interpolation3D_fix<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice, envelopeWidth); 
    applyProcessingFunctional(new Interpolation2D_fix<T,Descriptor>(wrapper,StaticAtomType), lattice.getBoundingBox(),lattice); 
  }
//*****************************
// interpolation velocity ends
//*****************************



//*********************************
//spreding fsi force to fluid nodes  using fix_lb style IB
//*********************************
  template<typename T, template<typename U> class Descriptor>
  class Spreading2D_fix: public BoxProcessingFunctional2D_L<T,Descriptor>{
      public:
      Spreading2D_fix(LammpsWrapper &wrapper_, plint StaticAtomType_):wrapper(wrapper_), StaticAtomType(StaticAtomType_){
        
      }
	  
	  
      virtual void process(Box2D domain, BlockLattice2D<T,Descriptor> &lattice){
        Dot2D offset = lattice.getLocation();
        
        T **x = wrapper.lmp->atom->x;
        T **f = wrapper.lmp->atom->f;
        int *mask = wrapper.lmp->atom->mask;
		int *type = wrapper.lmp->atom->type;
        plint nlocal = wrapper.lmp->atom->nlocal;
        //std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        int ix,iy;
        int ixp,iyp;
        double dx1,dy1;
        int isten,ii,jj;
        double r,rsq,weightx,weighty;
        double Ffp[16];
        //int k;
        //double unode[3];
		
		
        plint env =2;
        for (plint iS=0; iS<nlocal; iS++){
			
          //if(mask[iS] && type[iS]<MAX_TYPE){
			  if(mask[iS] && type[iS] != StaticAtomType){
			  
			  ix = (int)ceil(x[iS][0]-offset.x);
			  iy = (int)ceil(x[iS][1]-offset.y);
			  
			  
			  dx1 = x[iS][0] - offset.x -ix + 1.0;
			  dy1 = x[iS][1] - offset.y -iy + 1.0;
			  //dz1 = x[iS][2] - offset.z -iz + 1.0;
			  
             
             isten = 0;			 
            for (ii=-1;ii<3;ii++ ){
				rsq = (-dx1+ii)*(-dx1+ii);
    if(rsq>=4)
      weightx=0.0;
    else{
      r=sqrt(rsq);
      if(rsq>1){
	weightx=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
      } else{
	weightx=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
      }
    }
    for(jj=-1; jj<3; jj++){
      rsq=(-dy1+jj)*(-dy1+jj);
      if(rsq>=4)
	weighty=0.0;
      else{
	r=sqrt(rsq);
	if(rsq>1){
	  weighty=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
	} else{
	  weighty=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
	}
      }
      
	ixp = ix+ii;
	iyp = iy+jj;
	
				
				if(ixp<domain.x0-env || ixp>domain.x1+env ) continue;//ixp = domain.x1+2;
	            if(iyp<domain.y0-env || iyp>domain.y1+env) continue;//iyp = domain.y1+2;
	            
				
				Ffp[isten] = weightx*weighty;
				int shift = 0;
                Cell<T,Descriptor>& cell  = lattice.get(ixp+shift,iyp+shift);
                T *ff=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                
                ff[0] += Ffp[isten]*f[iS][0]; 
                ff[1] += Ffp[isten]*f[iS][1]; 
                
                cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,ff );
				
				isten++;
              }
          
			}
			
			
		  }
		}
	  }
        
      
      virtual Spreading2D_fix<T,Descriptor> * clone() const{
        return new Spreading2D_fix(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::staticVariables; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      plint StaticAtomType;

  };

  template<typename T, template<typename U> class Descriptor>
  void spreadForce2D_fix(MultiBlockLattice2D<T,Descriptor> &lattice,
                   LammpsWrapper &wrapper, plint StaticAtomType ){
    //plint envelopeWidth = 2;
    applyProcessingFunctional(new Spreading2D_fix<T,Descriptor>(wrapper,StaticAtomType), lattice.getBoundingBox(),lattice); 
  } 
  
//*********************************
//spreding force ends
//*********************************

//********************************
//force coupling
//********************************
  template<typename T, template<typename U> class Descriptor>
  class ForceFSI2D: public BoxProcessingFunctional2D_L<T,Descriptor>{
    public:
      ForceFSI2D(LammpsWrapper &wrapper_):wrapper(wrapper_){
        plint i,ifix(0),nfix;
        nfix = wrapper.lmp->modify->nfix;
        for (i=0;i<nfix;i++)
          if (strcmp(wrapper.lmp->modify->fix[i]->style,"fcm")==0) ifix=i;
          
        f_fcm = static_cast<LAMMPS_NS::FixFCM *>(wrapper.lmp->modify->fix[ifix]);
        f_fcm->grow_arrays(wrapper.lmp->atom->nmax);
        f_fcm->init();
        groupbit = f_fcm->groupbit;//new code
		dt = wrapper.lmp->update->dt;
		
		
	   
	   
      }
      virtual void process(Box2D domain, BlockLattice2D<T,Descriptor> &lattice){
        Dot2D offset = lattice.getLocation();
        TensorField2D<T,Descriptor<T>::d> velocity(lattice.getNx(),lattice.getNy());
        
        plint xl,yl,ix,iy,ii,jj;
        T rx,ry,wgt;
        T rho,dtfm;
        Array<T,2> us(0.,0.);
        Array<T,2> fsi(0.,0.);
        Array<T,2> uf(0.,0.);
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        T **f = wrapper.lmp->atom->f;
        T **fe = f_fcm->fexternal;
        int *mask = wrapper.lmp->atom->mask;
		double *mass = wrapper.lmp->atom->mass;
        int *type = wrapper.lmp->atom->type;
		molecule = wrapper.lmp->atom->molecule;
        
        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0);
        for(ix=domain.x0-2;ix<=domain.x1+2;ix++)
        for(iy=domain.y0-2;iy<=domain.y1+2;iy++){
        
          lattice.get(ix,iy).computeVelocity(velocity.get(ix,iy));
          
        }
        for (plint iS=0; iS<nlocal; iS++){
          if (mask[iS] & groupbit ){
            xl = floor(x[iS][0]); 
            yl = floor(x[iS][1]); 
            
            rx = x[iS][0] - xl;
            ry = x[iS][1] - yl;
            
            weight<T>(rx,wx);
            weight<T>(ry,wy);
           
            us[0] = us[1] = 0.0;
            rho=0.0;
            for (ii=0;ii<4;ii++ )
              for (jj=0;jj<4;jj++ ){
                
                  ix = xl-1 + ii - offset.x ;
                  iy = yl-1 + jj - offset.y ;
                  
                  if ( ix < domain.x0-2) continue; 
                  if ( iy < domain.y0-2) continue;
                  
                 
                  uf = velocity.get(ix,iy);
                  wgt = wx[ii]*wy[jj];
                  us[0] += wgt*uf[0];
                  us[1] += wgt*uf[1];
                  
                  rho += wgt*lattice.get(ix,iy).computeDensity();
                }
            
            fsi[0]=rho*(us[0]-v[iS][0]);       
            fsi[1]=rho*(us[1]-v[iS][1]);          
            
            
            fe[iS][0] = fsi[0]; // no accumulation for external force
            fe[iS][1] = fsi[1];
            
			
			//update position and velocity of atom
			if (molecule[iS] != 1) {
		 dtfm = dt / mass[type[iS]];
        v[iS][0] += dtfm * (f[iS][0] + fe[iS][0]);
        v[iS][1] += dtfm * (f[iS][1] + fe[iS][1]);
        //v[iS][2] += dtfm * f[iS][2];
        x[iS][0] += dt * v[iS][0];
        x[iS][1] += dt * v[iS][1];
        
			}

            for (ii=0;ii<4;ii++ )
              for (jj=0;jj<4;jj++ ){
                
                  ix = xl-1 + ii - offset.x ;
                  iy = yl-1 + jj - offset.y ;
                  
                 				  
				if (ix > domain.x1 || ix < domain.x0-2) continue; 
                if (iy > domain.y1 || iy < domain.y0-2) continue;
                
                  Cell<T,Descriptor>& cell  = lattice.get(ix,iy);
                  T *ff=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                  wgt = wx[ii]*wy[jj];
                  ff[0] -= wgt*fsi[0]; 
                  ff[1] -= wgt*fsi[1];                  
                  cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,ff );
                }
          }//mask[is]
        }
      }
      virtual ForceFSI2D<T,Descriptor> * clone() const{
        return new ForceFSI2D(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::staticVariables;
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      class LAMMPS_NS::FixFCM *f_fcm;
      plint groupbit;
	  LAMMPS_NS::tagint *molecule;
	  T dt;
  };


  template<typename T, template<typename U> class Descriptor>
  void forceCoupling2D(MultiBlockLattice2D<T,Descriptor> &lattice, LammpsWrapper &wrapper)
  {
    
    applyProcessingFunctional(new ForceFSI2D<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  }
//*********************************
// force coupling ends
//*********************************

}; /* namespace plb */

#endif 
