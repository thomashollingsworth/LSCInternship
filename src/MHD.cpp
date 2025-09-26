//Methods for the MHD Equations: overloaded for 1D or 2D
#include "MHD.h"

namespace MHD{

Vector8 primitiveToConservative(const Vector8& prim,double gamma){
    //Convert a 1D vector of density,v,p to the conserved quantities density,momentum,energy
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    Vector8 conserved;
    conserved[0]= prim[0];
    conserved[1]=prim[0]*prim[1];
    conserved[2]=prim[0]*prim[2];
    conserved[3]=prim[0]*prim[3];
    conserved[4]=prim[4]/((gamma-1)) +0.5*prim[0]*(prim[1]*prim[1]+prim[2]*prim[2]+prim[3]*prim[3])+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]);
    conserved[5]=prim[5];
    conserved[6]=prim[6];
    conserved[7]=prim[7];

    return conserved;
}

Vector8 conservativeToPrimitive(const Vector8& conserved,double gamma){
    //Convert a 1D vector of density,mom., energy to density,velocity,pressure
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    Vector8 prim;
    prim[0]= conserved[0];
    prim[1]=conserved[1]/conserved[0];
    prim[2]=conserved[2]/conserved[0];
    prim[3]=conserved[3]/conserved[0];
    prim[4] = (gamma-1)*(conserved[4]-0.5*(conserved[1]*conserved[1]+conserved[2]*conserved[2]+conserved[3]*conserved[3])/conserved[0]-0.5*(conserved[5]*conserved[5]+conserved[6]*conserved[6]+conserved[7]*conserved[7]));
    prim[5]=conserved[5];
    prim[6]=conserved[6];
    prim[7]=conserved[7];

    return prim;
}

Vector8 flux_x(const Vector8& u0,double gamma){
    Vector8 flux;
    Vector8 prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[1];
    flux[1]=prim[0]*prim[1]*prim[1]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[5]*prim[5];
    flux[2]= prim[0]*prim[1]*prim[2]-prim[5]*prim[6];
    flux[3]= prim[0]*prim[1]*prim[3]-prim[5]*prim[7];
    flux[4]=(u0[4]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[1]-(prim[1]*prim[5]+prim[2]*prim[6]+prim[3]*prim[7])*prim[5];

    flux[5]=u0[8];

    flux[6]=prim[6]*prim[1]-prim[5]*prim[2];
    flux[7]=prim[7]*prim[1]-prim[5]*prim[3];

    
    return flux;
}



Vector9 primitiveToConservative(const Vector9& prim,double gamma){
    //Convert a 1D vector of density,vx,vy,p to the conserved quantities density,momentum_x,momentum_y,energy
    Vector9 conserved;
    conserved[0]= prim[0];
    conserved[1]=prim[0]*prim[1];
    conserved[2]=prim[0]*prim[2];
    conserved[3]=prim[0]*prim[3];
    conserved[4]=prim[4]/((gamma-1)) +0.5*prim[0]*(prim[1]*prim[1]+prim[2]*prim[2]+prim[3]*prim[3])+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]);
    conserved[5]=prim[5];
    conserved[6]=prim[6];
    conserved[7]=prim[7];
    conserved[8]=prim[8];

    return conserved;

}

Vector9 conservativeToPrimitive(const Vector9& conserved,double gamma){
    //Convert a 1D vector of density,mom_x,mom_y, energy to density,vx,vy,pressure
    Vector9 prim;
    prim[0]= conserved[0];
    prim[1]=conserved[1]/conserved[0];
    prim[2]=conserved[2]/conserved[0];
    prim[3]=conserved[3]/conserved[0];
    prim[4] = (gamma-1)*(conserved[4]-0.5*(conserved[1]*conserved[1]+conserved[2]*conserved[2]+conserved[3]*conserved[3])/conserved[0]-0.5*(conserved[5]*conserved[5]+conserved[6]*conserved[6]+conserved[7]*conserved[7]));
    prim[5]=conserved[5];
    prim[6]=conserved[6];
    prim[7]=conserved[7];
    prim[8]=conserved[8];

    return prim;

}

Vector9 flux_x(const Vector9& u0,double gamma,double c_h){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    Vector9 flux;
    Vector9 prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[1];
    flux[1]=prim[0]*prim[1]*prim[1]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[5]*prim[5];
    flux[2]= prim[0]*prim[1]*prim[2]-prim[5]*prim[6];
    flux[3]= prim[0]*prim[1]*prim[3]-prim[5]*prim[7];
    flux[4]=(u0[4]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[1]-(prim[1]*prim[5]+prim[2]*prim[6]+prim[3]*prim[7])*prim[5];

    flux[5]=u0[8];
    flux[8]=u0[5]*c_h*c_h;
 
    
    flux[6]=prim[6]*prim[1]-prim[5]*prim[2];
    flux[7]=prim[7]*prim[1]-prim[5]*prim[3];

    
    return flux;

}

Vector9 flux_y(const Vector9& u0,double gamma,double c_h){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    Vector9 flux;
    Vector9 prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[2];
    flux[1]= prim[0]*prim[2]*prim[1]-prim[6]*prim[5];
    flux[2]=prim[0]*prim[2]*prim[2]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[6]*prim[6];
    flux[3]= prim[0]*prim[2]*prim[3]-prim[6]*prim[7];
    flux[4]=(u0[4]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[2]-(prim[1]*prim[5]+prim[2]*prim[6]+prim[3]*prim[7])*prim[6];
    flux[5]=-(prim[6]*prim[1]-prim[2]*prim[5]);
    flux[6]=u0[8];
    // flux[6]=0;
    flux[7]=-(prim[6]*prim[3]-prim[2]*prim[7]);
    
    flux[8]=u0[6]*c_h*c_h;
    // flux[8]=0;
    
    return flux;

}

}