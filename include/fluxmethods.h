#ifndef FLUXMETHODS_H
#define FLUXMETHODS_H

#include <array>


//Basic LF, Richtymer, FORCE numerical flux methods that act on a general statevector-like object and flux function (i.e. can be used for 1D/2D euler/MHD equations)
//Approximates the numerical flux at the interface between two states u_0 and u_1 given the flux function Flux, grid spacing dx, time step dt, and adiabatic index gamma
//Vector: type of the state vector (e.g. std::array<double,4> for 2D Euler, std::array<double,9> for 2D MHD)
//FluxFunc: type of the flux function (e.g. euler::flux, MHD_xflux, MHD_yflux)
//FluxFunc should be a callable that takes (const Vector&, double) and returns a Vector


template<typename Vector>
Vector getFluxRI(const Vector& flux0, const Vector& u_0,const Vector& flux1,const Vector& u_1,const double dx,const double dt,double gamma);


template<typename Vector>
Vector getFluxLF(const Vector& flux0,const Vector& u_0,const Vector& flux1,const Vector& u_1,const double dx,const double dt,double gamma);

template<typename Vector>
Vector getFluxFORCE(const Vector& flux0,const Vector& u_0,const Vector& flux1, const Vector& u_1,const double dx,const double dt,double gamma);

#endif