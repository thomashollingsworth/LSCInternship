#include "fluxmethods.h"


//Basic LF, Richtymer, FORCE fluxes that act on a general 1D Vector-like object with a general flux function


template<typename Vector>
Vector getFluxRI(const Vector& flux0, const Vector& u_0,const Vector& flux1,const Vector& u_1,const double dx,const double dt,double gamma){
    //returns the flux f_(i+1/2)
    const size_t size= u_0.size();
    Vector u_plus_half;
    for(size_t i=0;i<size;i++){
        u_plus_half[i]=0.5*(u_0[i]+u_1[i]) -0.5*dt/dx*(flux1[i]-flux0[i]);
    }
    return Flux(u_plus_half,gamma);
}
template<typename Vector>
Vector getFluxLF(const Vector& flux0,const Vector& u_0,const Vector& flux1,const Vector& u_1,const double dx,const double dt,double gamma){
    //returns the flux f_(i+1/2)
    const size_t size= u_0.size();
    Vector flux;
    for(size_t i=0;i<size;i++){
        flux[i]=0.5*dx/dt*(u_0[i]-u_1[i])+0.5*(flux0[i]+flux1[i]);
    }
    return flux;
}
template<typename Vector>
Vector getFluxFORCE(const Vector& flux0,const Vector& u_0,const Vector& flux1, const Vector& u_1,const double dx,const double dt,double gamma){
    const size_t size= u_0.size();
    Vector flux;
    Vector fluxLF=getFluxLF(flux0,u_0,flux1,u_1,dx,dt,gamma);
    Vector fluxRI=getFluxRI(flux0,u_0,flux1,u_1,dx,dt,gamma);
    for(size_t i=0;i<size;i++){
        flux[i]=0.5*(fluxLF[i]+fluxRI[i]);
    }
    return flux;
}

