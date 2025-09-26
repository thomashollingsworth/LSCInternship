#ifndef TIMESTEP_H
#define TIMESTEP_H
//Overloaded function to handle 1D or 2D
#include <vector>
#include "euler.h"
#include "MHD.h"
#include <initializer_list>

namespace euler{

using Vector4=std::array<double, 4>;
using Vector3=std::array<double, 3>;

//2D Euler Version
double computeTimeStep(const std::vector<std::vector<Vector4 > >& Prim,double dx,double dy, double C,double gamma,double num_xcells,double num_ycells);
//1D Euler Version
double computeTimeStep(const std::vector<Vector3 >& Prim,double dx,double C,double gamma,double num_xcells);
}


namespace MHD{
using Vector8=std::array<double, 8>;
using Vector9=std::array<double, 9>;
//The following functions can be useed for 1D or 2D MHD

template<typename StateVector>
double get_cfx(const double gamma, const StateVector& prim);
template<typename StateVector>
double get_cfy(const double gamma, const StateVector& prim);
template<typename StateVector>
double get_cfz(const double gamma, const StateVector& prim);
double get_c_h(const std::vector<std::vector<Vector9> >& prim, const int num_xcells,const int num_ycells,const double gamma);
double get_c_h(const std::vector<Vector8>& prim, const int num_xcells,const double gamma);
double computeTimeStep2D(double dx,double dy, double C,double c_h);
double computeTimeStep1D(double dx, double C,double c_h);

}
#endif