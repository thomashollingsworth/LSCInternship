#ifndef SAVE_H
#define SAVE_H

#include <iostream>
#include <fstream>
//Overloaded to handle 1 or 2D data

//1D version
template<typename StateVector>//Could be a std::array<double, num_params> for any num_params
void save_to_file(const std::vector<StateVector >& Prim,std::vector<double> x, double num_xcells, double gamma, std::string filename,std::string dir);

//2D version
template<typename StateVector>
void save_to_file(const std::vector<std::vector<StateVector > >& Prim,const std::vector<double >& x, double num_xcells,const std::vector<double>& y, double num_ycells, double gamma, std::string filename, std::string dir);

#endif