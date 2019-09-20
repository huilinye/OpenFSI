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

#ifndef IBM_LBM_3D_HH
#define IBM_LBM_3D_HH

#include "atom.h"
#include "modify.h"
#include "fix.h"
#include "fix_fcm.h"
#include "update.h"
#include <algorithm>

namespace plb {
  
  template<typename T>
  void weight(T r, std::vector<T> & w){
      T q = sqrt(1 + 4*r*(1-r));
      w[0] = (3 - 2*r - q)/8.0;
      w[1] = (3 - 2*r + q)/8.0;
      w[2] = (1 + 2*r + q)/8.0;
      w[3] = (1 + 2*r - q)/8.0;
  }
  
  template<typename T>
  T phi2(T r){
	  r = fabs(r);
	  r = 1.0 - r; 
	  if (r>0.0) return r;
	  else return 0.0;
  }

  //******************************
  //  interpolation velocity
  //******************************
  template<typename T, template<typename U> class Descriptor>
  class Interpolation3D: public BoxProcessingFunctional3D_L<T,Descriptor>{
    public:
      Interpolation3D(LammpsWrapper &wrapper_):wrapper(wrapper_){
        
        dt = wrapper.lmp->update->dt;
      }
      virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice){
        Dot3D offset = lattice.getLocation();
        TensorField3D<T,Descriptor<T>::d> velocity(lattice.getNx(),lattice.getNy(),lattice.getNz());
        
        plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        T rx,ry,rz,wgt,rho;
        Array<T,3> us(0.,0.,0.);
        Array<T,3> uf;
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        int *mask = wrapper.lmp->atom->mask;
        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
		plint env = 2;
        for(ix=domain.x0-env;ix<=domain.x1+env;ix++)
        for(iy=domain.y0-env;iy<=domain.y1+env;iy++)
        for(iz=domain.z0-env;iz<=domain.z1+env;iz++){
          lattice.get(ix,iy,iz).computeVelocity(velocity.get(ix,iy,iz));
        }
		
        for (plint iS=0; iS<nlocal; iS++){
          if (mask[iS] ){
            xl = floor(x[iS][0]); 
            yl = floor(x[iS][1]); 
            zl = floor(x[iS][2]);
            rx = x[iS][0] - xl;
            ry = x[iS][1] - yl;
            rz = x[iS][2] - zl;
            weight<T>(rx,wx);
            weight<T>(ry,wy);
            weight<T>(rz,wz);
            us[0] = us[1] = us[2]=0.0;
			rho = 0.0;
            for (ii=0;ii<4;ii++ )
              for (jj=0;jj<4;jj++ )
                for (kk=0;kk<4;kk++ ){
                  ix = xl-1 + ii - offset.x ;
                  iy = yl-1 + jj - offset.y ;
                  iz = zl-1 + kk - offset.z ;
                  if (ix > domain.x1+2 || ix < domain.x0-2) continue; 
                  if (iy > domain.y1+2 || iy < domain.y0-2) continue;
                  if (iz > domain.z1+2 || iz < domain.z0-2) continue;
				  
                 
                  uf = velocity.get(ix,iy,iz);
				  Cell<T,Descriptor>& cell  = lattice.get(ix,iy,iz);
                  T *bodyforce=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                  wgt = wx[ii]*wy[jj]*wz[kk];
                  us[0] += wgt*(uf[0] + 0.5*bodyforce[0]/lattice.get(ix,iy,iz).computeDensity());  //update velocity
                  us[1] += wgt*(uf[1] + 0.5*bodyforce[1]/lattice.get(ix,iy,iz).computeDensity());
                  us[2] += wgt*(uf[2] + 0.5*bodyforce[2]/lattice.get(ix,iy,iz).computeDensity());
				  
				  
				  rho += wgt*lattice.get(ix,iy,iz).computeDensity();
                }
								
            v[iS][0]=us[0];
            v[iS][1]=us[1];
            v[iS][2]=us[2];
            //Euler method to update position
            x[iS][0] += v[iS][0]*dt;
            x[iS][1] += v[iS][1]*dt;
            x[iS][2] += v[iS][2]*dt;
          }
        }
      }
      virtual Interpolation3D<T,Descriptor> * clone() const{
        return new Interpolation3D(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::nothing; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      
      T dt;
  };


  template<typename T, template<typename U> class Descriptor>
  void interpolateVelocity3D(MultiBlockLattice3D<T,Descriptor> &lattice, LammpsWrapper &wrapper)
  {
    
    applyProcessingFunctional(new Interpolation3D<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  }
//*****************************
// interpolation velocity ends
//*****************************

  //******************************
  //  interpolation velocity
  //******************************
  template<typename T, template<typename U> class Descriptor>
  class Interpolation3D_fix: public BoxProcessingFunctional3D_L<T,Descriptor>{
    public:
      Interpolation3D_fix(LammpsWrapper &wrapper_):wrapper(wrapper_){
        
        dt = wrapper.lmp->update->dt;
      }
      virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice){
        Dot3D offset = lattice.getLocation();
        TensorField3D<T,Descriptor<T>::d> velocity(lattice.getNx(),lattice.getNy(),lattice.getNz());
        
        //plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        //T rx,ry,rz,wgt,rho;
        Array<T,3> us(0.,0.,0.);
        Array<T,3> uf;
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        int *mask = wrapper.lmp->atom->mask;
        plint nlocal = wrapper.lmp->atom->nlocal;
        //std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
		
		plint ix,iy,iz;
        plint ixp,iyp,izp;
        T dx1,dy1,dz1;
        plint isten,ii,jj,kk;
        T r,rsq,weightx,weighty,weightz;
        T Ffp[64];
        //int k;
        //double unode[3];
		
		
        plint env = 2;
        for(ix=domain.x0-env;ix<=domain.x1+env;ix++)
        for(iy=domain.y0-env;iy<=domain.y1+env;iy++)
        for(iz=domain.z0-env;iz<=domain.z1+env;iz++){
          lattice.get(ix,iy,iz).computeVelocity(velocity.get(ix,iy,iz));
        }
		
        for (plint iS=0; iS<nlocal; iS++){
			
          if(mask[iS] ){
			  
			  us[0] = us[1] = us[2]=0.0;
			  ix = (int)ceil(x[iS][0]-offset.x);
			  iy = (int)ceil(x[iS][1]-offset.y);
			  iz = (int)ceil(x[iS][2]-offset.z);
			  
			  dx1 = x[iS][0] - offset.x -ix + 1.0;
			  dy1 = x[iS][1] - offset.y -iy + 1.0;
			  dz1 = x[iS][2] - offset.z -iz + 1.0;
	
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
      for(kk=-1; kk<3; kk++){
	rsq=(-dz1+kk)*(-dz1+kk);
	if(rsq>=4)
	  weightz=0.0;
	else{
	  r=sqrt(rsq);
	  if(rsq>1){
	    weightz=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
	  } else{
	    weightz=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
	  }
	}
	
	ixp = ix+ii;
	iyp = iy+jj;
	izp = iz+kk;
				
				if(ixp<domain.x0-env || ixp>domain.x1+env) continue;//ixp = domain.x1+2;
	            if(iyp<domain.y0-env || iyp>domain.y1+env) continue;//iyp = domain.y1+2;
	            if(izp<domain.z0-env || izp>domain.z1+env) continue;//izp = domain.z1+2;
		  
				Ffp[isten] = weightx*weighty*weightz;
                //std::cout<<" Ffp "<<isten <<" value "<<Ffp[isten]<<std::endl; 
                 int shift = 0;
                  uf = velocity.get(ixp+shift,iyp+shift,izp+shift);
				  
				  Cell<T,Descriptor>& cell  = lattice.get(ixp+shift,iyp+shift,izp+shift);
                  T *bodyforce=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                  
                  us[0] += Ffp[isten]*(uf[0] + 0.5*bodyforce[0]/lattice.get(ixp,iyp,izp).computeDensity());  //update velocity
                  us[1] += Ffp[isten]*(uf[1] + 0.5*bodyforce[1]/lattice.get(ixp,iyp,izp).computeDensity());
                  us[2] += Ffp[isten]*(uf[2] + 0.5*bodyforce[2]/lattice.get(ixp,iyp,izp).computeDensity()); 
				  
				  isten++;
				  //rho += wgt*lattice.get(ix,iy,iz).computeDensity();
                }
	}
			}
			
            v[iS][0]=us[0];
            v[iS][1]=us[1];
            v[iS][2]=us[2];
			//std::cout<<" usx "<<us[0] <<" usy "<<us[1]<<" usz "<<us[2]<<std::endl;
            //Euler method to update position
            x[iS][0] += v[iS][0]*dt;
            x[iS][1] += v[iS][1]*dt;
            x[iS][2] += v[iS][2]*dt;
          
        
      }
		}
	  }
      virtual Interpolation3D_fix<T,Descriptor> * clone() const{
        return new Interpolation3D_fix(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::nothing; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      
      T dt;
  };


  template<typename T, template<typename U> class Descriptor>
  void interpolateVelocity3D_fix(MultiBlockLattice3D<T,Descriptor> &lattice, LammpsWrapper &wrapper)
  {
    
    applyProcessingFunctional(new Interpolation3D_fix<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  }
//*****************************
// interpolation velocity ends
//*****************************

//******************************
  //  interpolation velocity
  //******************************
  template<typename T, template<typename U> class Descriptor>
  class Interpolation3D_fix2: public BoxProcessingFunctional3D_L<T,Descriptor>{
    public:
      Interpolation3D_fix2(LammpsWrapper &wrapper_):wrapper(wrapper_){
        
        dt = wrapper.lmp->update->dt;
      }
      virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice){
        Dot3D offset = lattice.getLocation();
        TensorField3D<T,Descriptor<T>::d> velocity(lattice.getNx(),lattice.getNy(),lattice.getNz());
        
        //plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        //T rx,ry,rz,wgt,rho;
        //Array<T,3> uf;
		Array<T,3> ufp[8];
		Array<T,3> us(0.,0.,0.);
		//T ufp[8][3];
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        int *mask = wrapper.lmp->atom->mask;
        plint nlocal = wrapper.lmp->atom->nlocal;
        //std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
		
		plint ix,iy,iz;
        plint ixp,iyp,izp;
        T dx1,dy1,dz1;
        plint isten,ii,jj,k;
        T r,rsq,weightx,weighty,weightz;
        T FfP[8];
        //int k;
        //double unode[3];
		
		
        plint env = 2;
        for(ix=domain.x0-env;ix<=domain.x1+env;ix++)
        for(iy=domain.y0-env;iy<=domain.y1+env;iy++)
        for(iz=domain.z0-env;iz<=domain.z1+env;iz++){
          lattice.get(ix,iy,iz).computeVelocity(velocity.get(ix,iy,iz));
        }
		
        for (plint iS=0; iS<nlocal; iS++){
			
          if(mask[iS] ){
			  
			  
			  ix = (int)ceil(x[iS][0]-offset.x);
			  iy = (int)ceil(x[iS][1]-offset.y);
			  iz = (int)ceil(x[iS][2]-offset.z);
			  
			  dx1 = x[iS][0] - offset.x -ix + 1.0;
			  dy1 = x[iS][1] - offset.y -iy + 1.0;
			  dz1 = x[iS][2] - offset.z -iz + 1.0;
			  
  //--------------------------------------------------------------------------
  // Calculate the interpolation weights
  //--------------------------------------------------------------------------
  FfP[0] = (1.-dx1)*(1.-dy1)*(1.-dz1);
  FfP[1] = (1.-dx1)*(1.-dy1)*dz1;
  FfP[2] = (1.-dx1)*dy1*(1.-dz1);
  FfP[3] = (1.-dx1)*dy1*dz1;
  FfP[4] = dx1*(1.-dy1)*(1.-dz1);
  FfP[5] = dx1*(1.-dy1)*dz1;
  FfP[6] = dx1*dy1*(1.-dz1);
  FfP[7] = dx1*dy1*dz1;
			  
	
	ixp = ix+1;
	iyp = iy+1;
	izp = iz+1;
				
				if(ixp<domain.x0-2 || ixp>domain.x1+2 ) continue;//ixp = domain.x1+2;
	            if(iyp<domain.y0-2 || iyp>domain.y1+2) continue;//iyp = domain.y1+2;
	            if(izp<domain.z0-2 || izp>domain.z1+2) continue;//izp = domain.z1+2;
		  
				  
                  
					  //extract velocity
				  ufp[0] = velocity.get(ix,iy,iz);
				  ufp[1] = velocity.get(ix,iy,izp);
				  ufp[2] = velocity.get(ix,iyp,iz);
				  ufp[3] = velocity.get(ix,iyp,izp);
				  ufp[4] = velocity.get(ixp,iy,iz);
				  ufp[5] = velocity.get(ixp,iy,izp);
				  ufp[6] = velocity.get(ixp,iyp,iz);
				  ufp[7] = velocity.get(ixp,iyp,izp);
				     // extract bodyforce
				  
				  //T *
				  //Cell<T,Descriptor>& cell  = lattice.get(ixp+shift,iyp+shift,izp+shift);
                 // T *bodyforce=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  T *bodyforce0,*bodyforce1,*bodyforce2,*bodyforce3,*bodyforce4,*bodyforce5,*bodyforce6,*bodyforce7;
				  Cell<T,Descriptor>& cell0  = lattice.get(ix,iy,iz);
                  bodyforce0=cell0.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell1  = lattice.get(ix,iy,izp);
                  bodyforce1=cell1.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell2  = lattice.get(ix,iyp,iz);
                  bodyforce2=cell2.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell3  = lattice.get(ix,iyp,izp);
                  bodyforce3=cell3.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell4  = lattice.get(ixp,iy,iz);
                  bodyforce4=cell4.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell5  = lattice.get(ixp,iy,izp);
                  bodyforce5=cell5.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell6  = lattice.get(ixp,iyp,iz);
                  bodyforce6=cell6.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell7  = lattice.get(ixp,iyp,izp);
                  bodyforce7=cell7.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  
				  
                  //wgt = wx[ii]*wy[jj]*wz[kk];
				  //std::cout<<" bodyforcex "<<bodyforce[0] <<" bodyforcey "<<bodyforce[1]<<" bodyforcez "<<bodyforce[2]<<std::endl; 
				  //std::cout<<" ufx "<<uf[0] <<" ufy "<<uf[1]<<" ufz "<<uf[2]<<std::endl; 
				  for (k=0; k<3; k++){
                  //us[0] += Ffp[isten]*(uf[0] + 0.5*bodyforce[0]/lattice.get(ixp,iyp,izp).computeDensity());  //update velocity
                  
				  us[k] = (ufp[0][k] + 0.5*bodyforce0[k]/lattice.get(ix,iy,iz).computeDensity())*FfP[0]
				        + (ufp[1][k] + 0.5*bodyforce1[k]/lattice.get(ix,iy,izp).computeDensity())*FfP[1]
						+ (ufp[2][k] + 0.5*bodyforce2[k]/lattice.get(ix,iyp,iz).computeDensity())*FfP[2]
						+ (ufp[3][k] + 0.5*bodyforce3[k]/lattice.get(ix,iyp,izp).computeDensity())*FfP[3]
						+ (ufp[4][k] + 0.5*bodyforce4[k]/lattice.get(ixp,iy,iz).computeDensity())*FfP[4]
						+ (ufp[5][k] + 0.5*bodyforce5[k]/lattice.get(ixp,iy,izp).computeDensity())*FfP[5]
						+ (ufp[6][k] + 0.5*bodyforce6[k]/lattice.get(ixp,iyp,iz).computeDensity())*FfP[6]
						+ (ufp[7][k] + 0.5*bodyforce7[k]/lattice.get(ixp,iyp,izp).computeDensity())*FfP[7];
				                  
	              }
			
			
            v[iS][0]=us[0];
            v[iS][1]=us[1];
            v[iS][2]=us[2];
			//std::cout<<" usx "<<us[0] <<" usy "<<us[1]<<" usz "<<us[2]<<std::endl;
            //Euler method to update position
            x[iS][0] += v[iS][0]*dt;
            x[iS][1] += v[iS][1]*dt;
            x[iS][2] += v[iS][2]*dt;
          
        
      }
		}
	  }
      virtual Interpolation3D_fix2<T,Descriptor> * clone() const{
        return new Interpolation3D_fix2(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::nothing; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      
      T dt;
  };


  template<typename T, template<typename U> class Descriptor>
  void interpolateVelocity3D_fix2(MultiBlockLattice3D<T,Descriptor> &lattice, LammpsWrapper &wrapper)
  {
    //plint envelopeWidth = 3;
    //applyProcessingFunctional(new Interpolation3D_fix<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice, envelopeWidth); 
    applyProcessingFunctional(new Interpolation3D_fix2<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  }
//*****************************
// interpolation velocity ends
//*****************************

//*********************************
//spreding fsi force to fluid nodes
//*********************************
  template<typename T, template<typename U> class Descriptor>
  class Spreading3D: public BoxProcessingFunctional3D_L<T,Descriptor>{
      public:
      Spreading3D(LammpsWrapper &wrapper_):wrapper(wrapper_){
        
      }
      virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice){
        Dot3D offset = lattice.getLocation();
        plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        T rx,ry,rz,wgt;
        //Array<T,3> ff(0.,0.,0.);
        T **x = wrapper.lmp->atom->x;
        T **f = wrapper.lmp->atom->f;
        int *mask = wrapper.lmp->atom->mask;
        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        for (plint iS=0; iS<nlocal; iS++){
          
          if(mask[iS] ){
          xl = floor(x[iS][0]); 
          yl = floor(x[iS][1]); 
          zl = floor(x[iS][2]);
          rx = x[iS][0] - xl;
          ry = x[iS][1] - yl;
          rz = x[iS][2] - zl;
          weight<T>(rx,wx);
          weight<T>(ry,wy);
          weight<T>(rz,wz);
          for (ii=0;ii<4;ii++ )
            for (jj=0;jj<4;jj++ )
              for (kk=0;kk<4;kk++ ){
                ix = xl-1 + ii - offset.x ;
                iy = yl-1 + jj - offset.y ;
                iz = zl-1 + kk - offset.z ;
                
         
                
                if (ix > domain.x1+2 || ix < domain.x0-2) continue; 
                if (iy > domain.y1+2 || iy < domain.y0-2) continue;
                if (iz > domain.z1+2 || iz < domain.z0-2) continue;
				  //if ( ix == domain.x0-3) ix = ix+1; //ensure the interpolation nodes locate within subprocessor
                  //if ( iy == domain.y0-3) iy = iy+1;
                  //if ( iz == domain.z0-3) iz = iz+1;
				
                Cell<T,Descriptor>& cell  = lattice.get(ix,iy,iz);
                T *ff=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                wgt = wx[ii]*wy[jj]*wz[kk];
                ff[0] += wgt*f[iS][0]; 
                ff[1] += wgt*f[iS][1]; 
                ff[2] += wgt*f[iS][2]; 
                cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,ff );
              }
          }//mask[iS]
        }
        
      }
      virtual Spreading3D<T,Descriptor> * clone() const{
        return new Spreading3D(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::staticVariables; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      

  };

  template<typename T, template<typename U> class Descriptor>
  void spreadForce3D(MultiBlockLattice3D<T,Descriptor> &lattice,
                   LammpsWrapper &wrapper ){
    //plint envelopeWidth = 2;
    applyProcessingFunctional(new Spreading3D<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  } 
  
//*********************************
//spreding force ends
//*********************************

//*********************************
//spreding fsi force to fluid nodes  using fix_lb style IB
//*********************************
  template<typename T, template<typename U> class Descriptor>
  class Spreading3D_fix: public BoxProcessingFunctional3D_L<T,Descriptor>{
      public:
      Spreading3D_fix(LammpsWrapper &wrapper_):wrapper(wrapper_){
        
      }
      virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice){
        Dot3D offset = lattice.getLocation();
        //plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        //T rx,ry,rz,wgt;
        //Array<T,3> ff(0.,0.,0.);
        T **x = wrapper.lmp->atom->x;
        T **f = wrapper.lmp->atom->f;
        int *mask = wrapper.lmp->atom->mask;
        plint nlocal = wrapper.lmp->atom->nlocal;
        //std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        int ix,iy,iz;
        int ixp,iyp,izp;
        double dx1,dy1,dz1;
        int isten,ii,jj,kk;
        double r,rsq,weightx,weighty,weightz;
        double Ffp[64];
        //int k;
        //double unode[3];
		
		
        plint env =2;
        for (plint iS=0; iS<nlocal; iS++){
			
          if(mask[iS] ){
			  
			  
			  ix = (int)ceil(x[iS][0]-offset.x);
			  iy = (int)ceil(x[iS][1]-offset.y);
			  iz = (int)ceil(x[iS][2]-offset.z);
			  
			  dx1 = x[iS][0] - offset.x -ix + 1.0;
			  dy1 = x[iS][1] - offset.y -iy + 1.0;
			  dz1 = x[iS][2] - offset.z -iz + 1.0;
			  
             //unode[0] = 0.0; unode[1] = 0.0; unode[2] = 0.0;
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
      for(kk=-1; kk<3; kk++){
	rsq=(-dz1+kk)*(-dz1+kk);
	if(rsq>=4)
	  weightz=0.0;
	else{
	  r=sqrt(rsq);
	  if(rsq>1){
	    weightz=(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*rsq))/8.;
	  } else{
	    weightz=(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*rsq))/8.;
	  }
	}
	ixp = ix+ii;
	iyp = iy+jj;
	izp = iz+kk;
				
				if(ixp<domain.x0-env || ixp>domain.x1+env ) continue;//ixp = domain.x1+2;
	            if(iyp<domain.y0-env || iyp>domain.y1+env) continue;//iyp = domain.y1+2;
	            if(izp<domain.z0-env || izp>domain.z1+env) continue;//izp = domain.z1+2;
				
				Ffp[isten] = weightx*weighty*weightz;
				int shift = 0;
                Cell<T,Descriptor>& cell  = lattice.get(ixp+shift,iyp+shift,izp+shift);
                T *ff=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                //wgt = wx[ii]*wy[jj]*wz[kk];
                ff[0] += Ffp[isten]*f[iS][0]; 
                ff[1] += Ffp[isten]*f[iS][1]; 
                ff[2] += Ffp[isten]*f[iS][2]; 
                cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,ff );
				
				isten++;
              }
          }
			}
			
			
		  }
		}
	  }
        
      
      virtual Spreading3D_fix<T,Descriptor> * clone() const{
        return new Spreading3D_fix(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::staticVariables; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      

  };

  template<typename T, template<typename U> class Descriptor>
  void spreadForce3D_fix(MultiBlockLattice3D<T,Descriptor> &lattice,
                   LammpsWrapper &wrapper ){
    //plint envelopeWidth = 2;
    applyProcessingFunctional(new Spreading3D_fix<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  } 
  
//*********************************
//spreding force ends
//*********************************

//*********************************
//spreding fsi force to fluid nodes  using fix_lb style IB
//*********************************
  template<typename T, template<typename U> class Descriptor>
  class Spreading3D_fix2: public BoxProcessingFunctional3D_L<T,Descriptor>{
      public:
      Spreading3D_fix2(LammpsWrapper &wrapper_):wrapper(wrapper_){
        
      }
      virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice){
        Dot3D offset = lattice.getLocation();
        //plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        //T rx,ry,rz,wgt;
        //Array<T,3> ff(0.,0.,0.);
        T **x = wrapper.lmp->atom->x;
        T **f = wrapper.lmp->atom->f;
        int *mask = wrapper.lmp->atom->mask;
        plint nlocal = wrapper.lmp->atom->nlocal;
        //std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        int ix,iy,iz;
        int ixp,iyp,izp;
        double dx1,dy1,dz1;
        int isten,ii,jj,kk;
        double r,rsq,weightx,weighty,weightz;
        double FfP[8];
        int k;
        double unode[3];
		
		
        
        for (plint iS=0; iS<nlocal; iS++){
			
          if(mask[iS] ){
			  
			  
			  ix = (int)ceil(x[iS][0]-offset.x);
			  iy = (int)ceil(x[iS][1]-offset.y);
			  iz = (int)ceil(x[iS][2]-offset.z);
			  
			  dx1 = x[iS][0] - offset.x -ix + 1.0;
			  dy1 = x[iS][1] - offset.y -iy + 1.0;
			  dz1 = x[iS][2] - offset.z -iz + 1.0;
			  
             //--------------------------------------------------------------------------
  // Calculate the interpolation weights
  //--------------------------------------------------------------------------
  FfP[0] = (1.-dx1)*(1.-dy1)*(1.-dz1);
  FfP[1] = (1.-dx1)*(1.-dy1)*dz1;
  FfP[2] = (1.-dx1)*dy1*(1.-dz1);
  FfP[3] = (1.-dx1)*dy1*dz1;
  FfP[4] = dx1*(1.-dy1)*(1.-dz1);
  FfP[5] = dx1*(1.-dy1)*dz1;
  FfP[6] = dx1*dy1*(1.-dz1);
  FfP[7] = dx1*dy1*dz1;
			  
	
	ixp = ix+1;
	iyp = iy+1;
	izp = iz+1;
				
				if(ixp<domain.x0-2 || ixp>domain.x1+2 ) continue;//ixp = domain.x1+2;
	            if(iyp<domain.y0-2 || iyp>domain.y1+2) continue;//iyp = domain.y1+2;
	            if(izp<domain.z0-2 || izp>domain.z1+2) continue;//izp = domain.z1+2;
		  
				  
                  
					
				     // extract bodyforce
				  
				  //T *
				  //Cell<T,Descriptor>& cell  = lattice.get(ixp+shift,iyp+shift,izp+shift);
                 // T *bodyforce=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  T *bodyforce0,*bodyforce1,*bodyforce2,*bodyforce3,*bodyforce4,*bodyforce5,*bodyforce6,*bodyforce7;
				  Cell<T,Descriptor>& cell0  = lattice.get(ix,iy,iz);
                  bodyforce0=cell0.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell1  = lattice.get(ix,iy,izp);
                  bodyforce1=cell1.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell2  = lattice.get(ix,iyp,iz);
                  bodyforce2=cell2.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell3  = lattice.get(ix,iyp,izp);
                  bodyforce3=cell3.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell4  = lattice.get(ixp,iy,iz);
                  bodyforce4=cell4.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell5  = lattice.get(ixp,iy,izp);
                  bodyforce5=cell5.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell6  = lattice.get(ixp,iyp,iz);
                  bodyforce6=cell6.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				  Cell<T,Descriptor>& cell7  = lattice.get(ixp,iyp,izp);
                  bodyforce7=cell7.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
				
				for(k=0;k<3;k++){
					bodyforce0[k] += FfP[0]*f[iS][k];
					bodyforce1[k] += FfP[1]*f[iS][k];
					bodyforce2[k] += FfP[2]*f[iS][k];
					bodyforce3[k] += FfP[3]*f[iS][k];
					bodyforce4[k] += FfP[4]*f[iS][k];
					bodyforce5[k] += FfP[5]*f[iS][k];
					bodyforce6[k] += FfP[6]*f[iS][k];
					bodyforce7[k] += FfP[7]*f[iS][k];
				}
				
				cell0.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,bodyforce0 );
				cell1.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,bodyforce1 );
				cell2.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,bodyforce2 );
				cell3.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,bodyforce3 );
				cell4.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,bodyforce4 );
				cell5.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,bodyforce5 );
				cell6.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,bodyforce6 );
				cell7.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,bodyforce7 );
				
				
              
          
			
			
			
		  }
		}
	  }
        
      
      virtual Spreading3D_fix2<T,Descriptor> * clone() const{
        return new Spreading3D_fix2(*this);
      }
      void getTypeOfModification(std::vector<modif::ModifT> & modified) const {
        modified[0]=modif::staticVariables; 
      }
      virtual BlockDomain::DomainT appliesTo() const{
        return BlockDomain::bulk;
      }
    private:
      LammpsWrapper &wrapper;
      

  };

  template<typename T, template<typename U> class Descriptor>
  void spreadForce3D_fix2(MultiBlockLattice3D<T,Descriptor> &lattice,
                   LammpsWrapper &wrapper ){
    //plint envelopeWidth = 2;
    applyProcessingFunctional(new Spreading3D_fix2<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  } 
  
//*********************************
//spreding force ends
//*********************************

//********************************
//force coupling
//********************************
  template<typename T, template<typename U> class Descriptor>
  class ForceFSI3D: public BoxProcessingFunctional3D_L<T,Descriptor>{
    public:
      ForceFSI3D(LammpsWrapper &wrapper_):wrapper(wrapper_){
        plint i,ifix(0),nfix;
        nfix = wrapper.lmp->modify->nfix;
        for (i=0;i<nfix;i++)
          if (strcmp(wrapper.lmp->modify->fix[i]->style,"fcm")==0) ifix=i;
          
        f_fcm = static_cast<LAMMPS_NS::FixFCM *>(wrapper.lmp->modify->fix[ifix]);
        f_fcm->grow_arrays(wrapper.lmp->atom->nmax);
        f_fcm->init();
        groupbit = f_fcm->groupbit;//new code
      }
      virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice){
        Dot3D offset = lattice.getLocation();
        TensorField3D<T,Descriptor<T>::d> velocity(lattice.getNx(),lattice.getNy(),lattice.getNz());
        
        plint xl,yl,zl,ix,iy,iz,ii,jj,kk;
        T rx,ry,rz,wgt;
        T rho;
        Array<T,3> us(0.,0.,0.);
        Array<T,3> fsi(0.,0.,0.);
        Array<T,3> uf(0.,0.,0.);
        T **x = wrapper.lmp->atom->x;
        T **v = wrapper.lmp->atom->v;
        //T **f = wrapper.lmp->atom->f;
        T **fe = f_fcm->fexternal;
		T dampcoe = f_fcm->dampcoe;
		int ntimestep = wrapper.lmp->update->ntimestep;
        int *mask = wrapper.lmp->atom->mask;
		
		double invdampcoe = 1.0;
	    invdampcoe = 1.0-exp(-ntimestep/dampcoe);
	    //if (ntimestep < dampcoe) invdampcoe = 1.0-exp(-ntimestep/dampcoe);

        plint nlocal = wrapper.lmp->atom->nlocal;
        std::vector<T> wx(4,0.0),wy(4,0.0),wz(4,0.0);
        for(ix=domain.x0-2;ix<=domain.x1+2;ix++)
        for(iy=domain.y0-2;iy<=domain.y1+2;iy++)
        for(iz=domain.z0-2;iz<=domain.z1+2;iz++){
          lattice.get(ix,iy,iz).computeVelocity(velocity.get(ix,iy,iz));
          //density(ix,iy,iz)=lattice.get(ix,iy,iz).computeDensity();
        }
        for (plint iS=0; iS<nlocal; iS++){
          if (mask[iS] & groupbit ){
            xl = floor(x[iS][0]); 
            yl = floor(x[iS][1]); 
            zl = floor(x[iS][2]);
            rx = x[iS][0] - xl;
            ry = x[iS][1] - yl;
            rz = x[iS][2] - zl;
            weight<T>(rx,wx);
            weight<T>(ry,wy);
            weight<T>(rz,wz);
            us[0] = us[1] = us[2]=0.0;
            rho=0.0;
            for (ii=0;ii<4;ii++ )
              for (jj=0;jj<4;jj++ )
                for (kk=0;kk<4;kk++ ){
                  ix = xl-1 + ii - offset.x ;
                  iy = yl-1 + jj - offset.y ;
                  iz = zl-1 + kk - offset.z ;
                  if ( ix < domain.x0-2) continue; 
                  if ( iy < domain.y0-2) continue;
                  if ( iz < domain.z0-2) continue;
                 
                  uf = velocity.get(ix,iy,iz);
                  wgt = wx[ii]*wy[jj]*wz[kk];
                  us[0] += wgt*uf[0];
                  us[1] += wgt*uf[1];
                  us[2] += wgt*uf[2];
                  rho += wgt*lattice.get(ix,iy,iz).computeDensity();
                }
            // assume rho = 1.
            //fsi[0]=us[0]-v[iS][0];       //Eqn.(5) in K. Aidun. Int J. Numer. Meth. Fluids 2010:62:765-783   
            //fsi[1]=us[1]-v[iS][1];          
            //fsi[2]=us[2]-v[iS][2]; 
            // using rho value
            fsi[0]=rho*(us[0]-v[iS][0])*invdampcoe;       //Eqn.(5) in K. Aidun. Int J. Numer. Meth. Fluids 2010:62:765-783   
            fsi[1]=rho*(us[1]-v[iS][1])*invdampcoe;          
            fsi[2]=rho*(us[2]-v[iS][2])*invdampcoe; 
            
            fe[iS][0] = fsi[0]; // no accumulation for external force
            fe[iS][1] = fsi[1];
            fe[iS][2] = fsi[2];
            

            for (ii=0;ii<4;ii++ )
              for (jj=0;jj<4;jj++ )
                for (kk=0;kk<4;kk++ ){
                  ix = xl-1 + ii - offset.x ;
                  iy = yl-1 + jj - offset.y ;
                  iz = zl-1 + kk - offset.z ;
                 				  
				if (ix > domain.x1 || ix < domain.x0-2) continue; 
                if (iy > domain.y1 || iy < domain.y0-2) continue;
                if (iz > domain.z1 || iz < domain.z0-2) continue;
				  
                  Cell<T,Descriptor>& cell  = lattice.get(ix,iy,iz);
                  T *ff=cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                  wgt = wx[ii]*wy[jj]*wz[kk];
                  ff[0] -= wgt*fsi[0]; 
                  ff[1] -= wgt*fsi[1]; 
                  ff[2] -= wgt*fsi[2]; 
                  cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt,Descriptor<T>::ExternalField::sizeOfForce,ff );
                }
          }//mask[is]
        }
      }
      virtual ForceFSI3D<T,Descriptor> * clone() const{
        return new ForceFSI3D(*this);
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
  };


  template<typename T, template<typename U> class Descriptor>
  void forceCoupling3D(MultiBlockLattice3D<T,Descriptor> &lattice, LammpsWrapper &wrapper)
  {
    //plint envelopeWidth = 2;
    //applyProcessingFunctional(new Interpolation3D<T>(wrapper), velocity.getBoundingBox(),velocity, envelopeWidth); 
    applyProcessingFunctional(new ForceFSI3D<T,Descriptor>(wrapper), lattice.getBoundingBox(),lattice); 
  }
//*********************************
// force coupling ends
//*********************************

}; /* namespace plb */

#endif 
