#ifndef CONFIG_H
#define CONFIG_H

struct Config{

    //Structure that stores the type of Boundaries for an array and is used when update_BC functions are called
    struct BC {
    enum class Type{
        periodic,
        transmissive,//2nd order Neumann in all variables
        reflective,
        conductive
    };
    Type x0 = Type::periodic; //Default assigned values
    Type xf = Type::periodic;
    //(Only have meaning for 2D systems)
    Type y0 = Type::periodic;
    Type yf = Type::periodic;};


    //Stores the type of data in each cell
    enum class DType{
        Euler1D,
        Euler2D,
        MHD1D,
        MHD2D};

    BC boundary_conditions;
    DType data_type= DType::Euler1D; //Useful for generalising BC updates and flux calcs
    
    
    //Parameters used in updates
    double dx=0;
    double dy=0;
    double dt=0;

    //Parameters used all the time
    int num_xcells=0;
    int num_ycells=0;

    //Parameters used in initialisation
    double x0=0;
    double xf=0;
    double y0=0;
    double yf=0;
};

#endif