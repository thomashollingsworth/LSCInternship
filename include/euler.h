//Methods for the Euler Equations: overloaded for 1D or 2D
#ifndef EULER_H
#define EULER_H
#include <array>

namespace euler{
std::array<double, 3> primitiveToConservative(std::array<double, 3> prim,double gamma);
std::array<double, 3> conservativeToPrimitive(std::array<double, 3> conserved,double gamma);
std::array<double, 3> flux(std::array<double, 3> u0,double gamma);

typedef std::array< double ,4> StateVector;
StateVector primitiveToConservative(const StateVector& prim,double gamma);
StateVector conservativeToPrimitive(const StateVector& conserved,double gamma);
StateVector flux(const StateVector& u0,double gamma);
}
#endif