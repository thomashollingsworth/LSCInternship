#include "timestep.h"


//Overloaded function to handle 1D or 2D 
namespace euler{

using Vector4=std::array<double, 4>;
using Vector3=std::array<double, 3>;

//2D Euler Version
double computeTimeStep(const std::vector<std::vector<Vector4 > >& Prim,double dx,double dy, double C,double gamma,double num_xcells,double num_ycells){
    Vector4 prim0=Prim[0][0];
    double max=std::abs(std::sqrt(prim0[1]*prim0[1]+prim0[2]*prim0[2]))+std::sqrt(gamma*prim0[3]/prim0[0]);
    Vector4 prim;
    //Calculate max wave speed for all real cells
    for(size_t i=2; i<num_xcells+2;i++){
        for(size_t j=2; j<num_ycells+2;j++) 
        prim=Prim[i][j];
        double new_max=std::abs(std::sqrt(prim[1]*prim[1]+prim[2]*prim[2]))+std::sqrt(gamma*prim[3]/prim[0]);
        if(new_max> max){max=new_max;}
    }
 return C*std::min(dx,dy)/max;
}

//1D Euler Version
double computeTimeStep(const std::vector<Vector3 >& Prim,double dx,double C,double gamma,double num_xcells){
    Vector3 prim0=Prim[0];
   
    double max=std::abs(prim0[1])+std::sqrt(gamma*prim0[2]/prim0[0]);
    
    Vector3 prim;
    
    for(size_t i=2; i<num_xcells+2;i++){
        prim=Prim[i];
        double new_max=std::abs(prim[1])+std::sqrt(gamma*prim[2]/prim[0]);
        if(new_max> max){max=new_max;}  
    }
 return C*dx/max;
}

}



namespace MHD{

using Vector8=std::array<double, 8>;
using Vector9=std::array<double, 9>;
//The following functions can be useed for 1D or 2D MHD

template<typename StateVector>
double get_cfx(const double gamma, const StateVector& prim){
    
    double factor= (gamma*prim[4]+prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])/prim[0];
    double c_fx=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*prim[5]*prim[5]/(prim[0]*prim[0]))));
    return c_fx;

}
template<typename StateVector>
double get_cfy(const double gamma, const StateVector& prim){
    double factor= (gamma*prim[4]+prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])/prim[0];
    double c_fy=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*prim[6]*prim[6]/(prim[0]*prim[0]))));
    return c_fy;
}
template<typename StateVector>
double get_cfz(const double gamma, const StateVector& prim){
    double factor= (gamma*prim[4]+prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])/prim[0];
    double c_fz=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*prim[7]*prim[7]/(prim[0]*prim[0]))));
    return c_fz;
}

//This function is overloaded to handle either 2D or 1D
double get_c_h(const std::vector<std::vector<Vector9> >& prim, const int num_xcells,const int num_ycells,const double gamma){
    //2D specific operation (c_h is just max wavespeed in system)
    double max_ch=0.;
    for(size_t i=2; i<num_xcells+2;i++){
        for(size_t j=2; j<num_ycells+2;j++){
        double c_fx=get_cfx(gamma,prim[i][j]);
        double c_fy = get_cfy(gamma,prim[i][j]);
        double c_fz = get_cfz(gamma,prim[i][j]);

        double new_max=std::max({std::abs(prim[i][j][1])+c_fx,std::abs(prim[i][j][2])+c_fy,std::abs(prim[i][j][3])+c_fz});
        if(new_max> max_ch){
            max_ch=new_max;
        }   
    }
}
return max_ch;
}

double get_c_h(const std::vector<Vector8>& prim, const int num_xcells,const double gamma){
    //1D specific operation (c_h is just max wavespeed in system)
    double max_ch=0.;
    for(size_t i=2; i<num_xcells+2;i++){
        double c_fx=get_cfx(gamma,prim[i]);
        double new_max=std::sqrt(prim[i][1]*prim[i][1]+prim[i][2]*prim[i][2]+prim[i][3]*prim[i][3])+c_fx;
        if(new_max> max_ch){
            max_ch=new_max;
        }   
    }
    return max_ch;
}

double computeTimeStep2D(double dx,double dy, double C,double c_h){
 //Making an adjustment here
return C* std::min(dx,dy)/(c_h);
}

double computeTimeStep1D(double dx, double C,double c_h){
 //Making an adjustment here
return C* dx/(c_h);
}

}