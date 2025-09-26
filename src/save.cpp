#include "save.h"
//Includes an overload for either 2D or 1D data can work for MHD or Euler
//No plotting just pure csv data

//2D Version
template<typename StateVector>//Could be a std::array<double, num_params> for any num_params i.e. can handle euler or MHD
void save_to_file(const std::vector<std::vector<StateVector > >& Prim,const std::vector<double >& x, double num_xcells,const std::vector<double>& y, double num_ycells, double gamma, std::string filename, std::string dir){
    std::filesystem::create_directories(dir);
    const size_t num_params= Prim[0][0].size();
    
    filename=dir+filename;

    std::ofstream outfile(filename);

    for(size_t j=2; j<num_ycells+2;j++){
        for(size_t i=2; i<num_xcells+2;i++){
            outfile << x[i] << " " << y[j] << " ";
            for(size_t k=0; k<num_params;k++){outfile << Prim[i][j][k] << " ";}
            outfile<< "\n";}
    outfile<<"\n"; //Blank line between y rows for correct gnuplot formatting
}
}

//1D Version
template<typename StateVector>
void save_to_file(const std::vector<StateVector >& Prim,std::vector<double> x, double num_xcells, double gamma, std::string filename,std::string dir){
    std::filesystem::create_directories(dir); 
    const size_t num_params= Prim[0][0].size();
    
    filename=dir+filename;

    std::ofstream outfile(filename);

    for(size_t i=2; i<num_xcells+2;i++){
            outfile << x[i] << " ";
            for(size_t k=0; k<num_params;k++){outfile << Prim[i][k] << " ";}
            outfile<< "\n";}
}
