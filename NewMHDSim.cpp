#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <cstdlib>
#include <filesystem>  // C++17
#include <initializer_list>



struct ConservedStateVector{
    //A structure to hold and access a single StateVector member
    //The data is stored contiguously (std::array<double,N>)
    //Each component can be accessed by its name or by the [] operator
    //Member function converts to PrimitiveStateVector
    
    std::array<double,9> data;

    //Enable access by index (makes iteration easier)
    double& operator[](size_t i) { return data[i]; }
    const double& operator[](size_t i) const { return data[i]; }
    
    //These functions are called on the structure instance to reference a single component 
    //(aids readibility + debugging)
    
    double& density()   {return data[0];}
    double& mom_x()     {return data[1];}
    double& mom_y()     {return data[2];}
    double& mom_z()     {return data[3];}
    double& B_x()       {return data[5];}
    double& B_y()       {return data[6];}
    double& B_z()       {return data[7];}
    double& psi()       {return data[8];}

    







}
