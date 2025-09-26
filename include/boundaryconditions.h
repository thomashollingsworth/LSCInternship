#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H
#include <vector>
#include "config.h"

//update boundaries individually or all together for 1D/2D MHD/Euler
//update boundaries based on specific Config (which contains datatype and boundary type for each boundary)

template<typename StateVector> //2D array version
void updateBCs(const Config& cfg,std::vector<std::vector<StateVector> >& u);
template<typename StateVector> //1D array version (overloaded)
void updateBCs(const Config& cfg, std::vector<StateVector>& u);

template<typename Array> //Could be a 1D or 2D array
void update_bcs_x0(const Config& cfg, Array& u);
template<typename Array> //Could be a 1D or 2D array
void update_bcs_xf(const Config& cfg, Array& u);
template<typename StateVector>//y updates are only defined for 2D arrays of StateVectors
void update_bcs_y0(const Config& cfg,std::vector<std::vector<StateVector> >& u);
template<typename StateVector>//y updates are only defined for 2D arrays of StateVectors
void update_bcs_yf(const Config& cfg,std::vector<std::vector<StateVector> >& u);




#endif