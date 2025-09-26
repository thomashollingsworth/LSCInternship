#include "config.h"
#include <vector>
#include <tuple>
//Builds a config structure and initial data arrays x,y,u from specified initial conditions

//Initial conditions
enum class IC{
    Euler_ShockTube1D,
    Euler_ShockTube2DX,
    Euler_ShockTube2DY,
    MHD_BrioWu1D,
    MHD_BrioWu2DX,
    MHD_BrioWu2DY,
    OrszagTang,
    KelvinHelmholtz
};

//This generalises for x or y array
double dx(const int num_xcells, const double x0, const double xf);
std::vector<double> set_x(const int num_xcells, const double dx,const double x0);

template<typename StateVector>//Returns x,u
std::tuple<std::vector<double>, std::vector<StateVector>> initialise1D(Config& cfg,IC initial_conditions, const double num_xcells);

template<typename StateVector>
std::vector<std::vector<StateVector> > initialise2D(Config& cfg,IC initial_conditions, const double num_xcells, const double num_ycells);