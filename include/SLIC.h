#ifndef SLIC_H
#define SLIC_H

//Performs SLIC method on a general std::array like object
//Can be used for euler or MHD or any N-param system
#include <tuple>

template<typename StateVector>
StateVector getDelta_plus(const StateVector& u_0, const StateVector& u_plus);
template<typename StateVector>
StateVector getDelta_minus(const StateVector& u_minus, const StateVector& u_0);
template<typename StateVector>
StateVector getDelta(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w);
template<typename StateVector>
StateVector get_r(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus);
template<typename StateVector>
StateVector getXi_Minbee(const StateVector& r);
template<typename StateVector>
std::tuple< StateVector,StateVector  > calc_ubar(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w);

template<typename StateVector>
std::tuple< StateVector,StateVector > calc_ubar_plus(const StateVector& ubarL, const StateVector& ubarR,const StateVector& fluxL, const StateVector& fluxR, const double dt, const double dx);
#endif