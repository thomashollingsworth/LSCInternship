#include "initialise.h"


double dx(const int num_xcells, const double x0, const double xf){
    return (xf-x0)/num_xcells;
}
std::vector<double> set_x(const int num_xcells, const double dx,const double x0){
    std::vector<double> x_array;
    for(size_t i=0;i<num_xcells;i++){
        x_array[i]=x0+i*dx;
    }
    return x_array;
}



template<typename StateVector>
std::tuple<std::vector<double>, std::vector<StateVector>> initialise1D(Config& cfg,IC initial_conditions, const double num_xcells){
    cfg.num_xcells=num_xcells;
    std::vector<double> x_array(cfg.num_xcells);
    std::vector<StateVector> u0(cfg.num_xcells+4);
    switch(initial_conditions){
        case IC::Euler_ShockTube1D:
        //Build config
        cfg.data_type=Config::DType::Euler1D;
        cfg.boundary_conditions.x0=Config::BC::Type::transmissive;
        cfg.boundary_conditions.xf=Config::BC::Type::transmissive;
        cfg.x0=0.;
        cfg.xf=800.;
        cfg.dx=dx(cfg.num_xcells,cfg.x0,cfg.xf);
    
        //Use config to construct initial conditions
        x_array=set_x(cfg.num_xcells,cfg.dx,cfg.x0);

        //Set initial condition in primitive coords
        for(size_t i=2;i<cfg.num_xcells+2;i++){
        if(x[i]<=0.5){
            prim[0]=1.;
            prim[1]=0.75;
            prim[2]=1.0;
        }
        else{
            prim[0]=0.125;
            prim[1]=0.;
            prim[2]=0.1;
        }}



  
    }
    


}