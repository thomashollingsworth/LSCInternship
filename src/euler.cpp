//Methods for the Euler Equations: overloaded for 1D or 2D
#include "euler.h"

namespace euler{

typedef std::array< double ,3> Vector3;

Vector3 primitiveToConservative(const Vector3& prim,double gamma){
    //Convert a 1D vector of density,v,p to the conserved quantities density,momentum,energy
    Vector3 conserved;
    conserved[0]= prim[0];
    conserved[1]=prim[0]*prim[1];
    conserved[2] = prim[2]/((gamma-1)) +0.5*prim[0]*prim[1]*prim[1];

    return conserved;
}

Vector3 conservativeToPrimitive(const Vector3& conserved,double gamma){
    //Convert a 1D vector of density,mom., energy to density,velocity,pressure
    Vector3 prim;
    prim[0]= conserved[0];
    prim[1]=conserved[1]/conserved[0];
    prim[2] = (gamma-1)*(conserved[2]-0.5*(conserved[1]*conserved[1])/conserved[0]);

    return prim;
}

Vector3 flux_x(const Vector3& u0,double gamma){
    Vector3 flux;
    Vector3 prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[1];
    flux[1]=prim[0]*prim[1]*prim[1] +prim[2];
    flux[2]= (u0[2]+prim[2])*prim[1];
    return flux;
}

typedef std::array< double ,4> Vector4;

Vector4 primitiveToConservative(const Vector4& prim,double gamma){
    //Convert a 1D vector of density,vx,vy,p to the conserved quantities density,momentum_x,momentum_y,energy
    Vector4 conserved;
    conserved[0]= prim[0];
    conserved[1]=prim[0]*prim[1];
    conserved[2]=prim[0]*prim[2];
    conserved[3] = prim[3]/((gamma-1)) +0.5*prim[0]*(prim[1]*prim[1]+prim[2]*prim[2]);

    return conserved;

}

Vector4 conservativeToPrimitive(const Vector4& conserved,double gamma){
    //Convert a 1D vector of density,mom_x,mom_y, energy to density,vx,vy,pressure
    Vector4 prim;
    prim[0]= conserved[0];
    prim[1]=conserved[1]/conserved[0];
    prim[2]=conserved[2]/conserved[0];
    prim[3] = (gamma-1)*(conserved[3]-0.5*(conserved[1]*conserved[1]+conserved[2]*conserved[2])/conserved[0]);

    return prim;

}

Vector4 flux_x(const Vector4& u0,double gamma){
    Vector4 flux;
    Vector4 prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[1];
    flux[1]=prim[0]*prim[1]*prim[1] +prim[3];
    flux[2]=prim[0]*prim[1]*prim[2];
    flux[3]= (u0[3]+prim[3])*prim[1];
    return flux;

}

Vector4 flux_y(const Vector4& u0,double gamma){
    Vector4 flux;
    Vector4 prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[2];
    flux[1]=prim[0]*prim[2]*prim[1];
    flux[2]=prim[0]*prim[2]*prim[2] +prim[3];
    flux[3]= (u0[3]+prim[3])*prim[2];
    return flux;

}

}