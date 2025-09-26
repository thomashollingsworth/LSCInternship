//Methods for the MHD Equations: overloaded for 1D or 2D
#ifndef MHD_H
#define MHD_H
#include <array>

namespace MHD{
typedef std::array< double ,8> Vector8;
typedef std::array< double ,9> Vector9;

Vector8 primitiveToConservative(const Vector8& prim,double gamma);
Vector8 conservativeToPrimitive(const Vector8& conserved,double gamma);
Vector8 flux_x(const Vector8& u0,double gamma);

Vector9 primitiveToConservative(const Vector9& prim,double gamma);
Vector9 conservativeToPrimitive(const Vector9& conserved,double gamma);
Vector9 flux_x(const Vector9& u0,double gamma,double c_h);
Vector9 flux_y(const Vector9& u0,double gamma,double c_h);
}
#endif