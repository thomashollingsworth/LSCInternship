#include "boundaryconditions.h"


//update boundaries individually or all together for 1D/2D MHD/Euler
//update boundaries based on specific Config (which contains datatype and boundary type for each boundary)

template<typename StateVector> //2D array version
void updateBCs(const Config& cfg,std::vector<std::vector<StateVector> >& u){
    
    update_bcs_x0(cfg.boundary_conditions,cfg.data_type,u);
    update_bcs_xf(cfg.boundary_conditions,cfg.data_type,u);

    update_bcs_y0(cfg.boundary_conditions,cfg.data_type,u);
    update_bcs_yf(cfg.boundary_conditions,cfg.data_type,u);
    };//Updates all BCs for an array

template<typename StateVector> //1D array version (overloaded)
void updateBCs(const Config& cfg, std::vector<StateVector>& u){
    update_bcs_x0(cfg.conditions,cfg.data_type,u);
    update_bcs_xf(cfg.conditions,cfg.data_type,u);
    };//Updates all BCs for an array


template<typename Array> //Could be a 1D or 2D array
void update_bcs_x0(const Config& cfg, Array& u){
size_t num_xcells= u.size()-4; //For 2 ghost cells each side (makes checking BCs clearer)
size_t num_ycells;
if(cfg.data_type==Config::DType::Euler2D || cfg.data_type==Config::DType::MHD2D){// num_ycells is defined for 2D arrays only
    size_t num_ycells=u[0].size()-4;
}
switch(conditions.x0){
        case Config::BC::Type::periodic:
        u[0]=u[num_xcells];
        u[1]=u[num_xcells+1];
        break;

        case Config::BC::Type::transmissive:
        u[0]=u[3];
        u[1]=u[2];
        break;
        
        case Config::BC::Type::reflective:
        u[0]=u[3];
        u[1]=u[2];
        //Reflects normal velocity for both Euler and MHD, 1D or 2D
        switch(cfg.data_type){
            case(Config::DType::Euler1D || Config::DType::MHD2D):
                u[0][1]*=-1;
                u[1][1]*=-1;
                break;
            case(Config::DType::Euler2D || Config::DType::MHD2D):
                for(size_t j=0;j<num_ycells+4;j++){
                    u[0][j][1]*=-1;
                    u[1][j][1]*=-1;
                } break;};
        break;

        case Config::BC::Type::conductive: //Only valid for B_n(0)=0
        u[0]=u[3];
        u[1]=u[2];

        switch(cfg.data_type){
            case Config::DType::Euler1D:
                u[0][1]*=-1;
                u[1][1]*=-1;
                break;
            case Config::DType::Euler2D:
                for(size_t j=0;j<num_ycells+4;j++){
                    u[0][j][1]*=-1;
                    u[1][j][1]*=-1;
                } break;
            case Config::DType::MHD1D:
                u[0][1]*=-1;
                u[1][1]*=-1;
                //Reflecting B_n
                u[0][5]*=-1;
                u[1][5]*=-1
                break;
            case Config::DType::MHD2D:
                for(size_t j=0;j<num_ycells+4;j++){
                    u[0][j][1]*=-1;
                    u[1][j][1]*=-1;
                    u[0][j][5]*=-1;
                    u[1][j][5]*=-1;
                } break;};break;}
}

template<typename Array> //Could be a 1D or 2D array
void update_bcs_xf(const Config& cfg, Array& u){
size_t num_xcells= u.size()-4; //For 2 ghost cells each side (makes checking BCs clearer)
size_t num_ycells;
if(cfg.data_type==Config::DType::Euler2D || cfg.data_type==Config::DType::MHD2D){// num_ycells is defined for 2D arrays only
    size_t num_ycells=u[0].size()-4;
}
switch(conditions.xf){
        case Config::BC::Type::periodic:
        u[num_xcells+2]=u[2];
        u[num_xcells+3]=u[3];
        break;

        case Config::BC::Type::transmissive:
        u[num_xcells+2]=u[num_xcells+1];
        u[num_xcells+3]=u[num_xcells];
        break;
        
        case Config::BC::Type::reflective:
        u[num_xcells+2]=u[num_xcells+1];
        u[num_xcells+3]=u[num_xcells];
        //Reflects normal velocity for both Euler and MHD, 1D or 2D
        switch(cfg.data_type){
            case(Config::DType::Euler1D || Config::DType::MHD2D):
                u[num_xcells+2][1]*=-1;
                u[num_xcells+3][1]*=-1;
                break;
            case(Config::DType::Euler2D || Config::DType::MHD2D):
                for(size_t j=0;j<num_ycells+4;j++){
                    u[num_xcells+2][j][1]*=-1;
                    u[num_xcells+3][j][1]*=-1;
                } break;};
        break;

        case Config::BC::Type::conductive: //Only valid for B_n(0)=0
        u[num_xcells+2]=u[num_xcells+1];
        u[num_xcells+3]=u[num_xcells];
        switch(cfg.data_type){
            case Config::DType::Euler1D:
                u[num_xcells+2][1]*=-1;
                u[num_xcells+3][1]*=-1;
                break;
            case Config::DType::Euler2D:
                for(size_t j=0;j<num_ycells+4;j++){
                    u[num_xcells+2][j][1]*=-1;
                    u[num_xcells+3][j][1]*=-1;
                } break;
            case Config::DType::MHD1D:
                u[num_xcells+2][1]*=-1;
                u[num_xcells+3][1]*=-1;
                //Reflecting B_n
                u[num_xcells+2][5]*=-1;
                u[num_xcells+3][5]*=-1
                break;
            case Config::DType::MHD2D:
                for(size_t j=0;j<num_ycells+4;j++){
                    u[num_xcells+2][j][1]*=-1;
                    u[num_xcells+3][j][1]*=-1;
                    u[num_xcells+2][j][5]*=-1;
                    u[num_xcells+3][j][5]*=-1;
                } break;};break;}
}

template<typename StateVector>//y updates are only defined for 2D arrays of StateVectors
void update_bcs_y0(const Config& cfg,std::vector<std::vector<StateVector> >& u){
    size_t num_xcells= u.size()-4; //For 2 ghost cells each side (makes checking BCs clearer)
    size_t num_ycells=u[0].size()-4;

    switch(conditions.y0){
        case Config::BC::Type::periodic:
            for(size_t i=0; i<num_xcells+4;i++){
                u[i][0]=u[i][num_ycells];
                u[i][1]=u[i][num_ycells+1];
            }break;
        case Config::BC::Type::transmissive:
            for(size_t i=0; i<num_xcells+4;i++){
                u[i][0]=u[i][3];
                u[i][1]=u[i][2];
            }break;
        case Config::BC::Type::reflective:
            for(size_t i=0; i<num_xcells+4;i++){
                u[i][0]=u[i][3];
                u[i][1]=u[i][2];

                //reflect y-velocity
                u[i][0][2]*=-1;
                u[i][1][2]*=-1;
            }break;
        case Config::BC::Type::conductive:
            switch(cfg.data_type){
                case Config::DType::MHD2D:
                for(size_t i=0; i<num_xcells+4;i++){
                    u[i][0]=u[i][3];
                    u[i][1]=u[i][2];

                    //reflect y-velocity
                    u[i][0][2]*=-1;
                    u[i][1][2]*=-1;

                    //reflect B_y
                    u[i][0][6]*=-1;
                    u[i][1][6]*=-1;}break;
                case Config::DType::Euler2D:
                for(size_t i=0; i<num_xcells+4;i++){
                    u[i][0]=u[i][3];
                    u[i][1]=u[i][2];

                    //reflect y-velocity
                    u[i][0][2]*=-1;
                    u[i][1][2]*=-1;}break;
        break;}
    }
}

template<typename StateVector>//y updates are only defined for 2D arrays of StateVectors
void update_bcs_yf(const Config& cfg,std::vector<std::vector<StateVector> >& u){
    size_t num_xcells= u.size()-4; //For 2 ghost cells each side (makes checking BCs clearer)
    size_t num_ycells=u[0].size()-4;

    switch(conditions.yf){
        case Config::BC::Type::periodic:
            for(size_t i=0; i<num_xcells+4;i++){
                u[i][num_ycells+2]=u[i][2];
                u[i][num_ycells+3]=u[i][3];
            }break;
        case Config::BC::Type::transmissive:
            for(size_t i=0; i<num_xcells+4;i++){
                u[i][num_ycells+2]=u[i][num_ycells+1];
                u[i][num_ycells+3]=u[i][num_ycells];
            }break;
        case Config::BC::Type::reflective:
            for(size_t i=0; i<num_xcells+4;i++){
                u[i][num_ycells+2]=u[i][num_ycells+1];
                u[i][num_ycells+3]=u[i][num_ycells];

                //reflect y-velocity
                u[i][num_ycells+2][2]=u[i][num_ycells+1];
                u[i][num_ycells+3][2]=u[i][num_ycells];
            }break;
        case Config::BC::Type::conductive:
            switch(cfg.data_type){
                case Config::DType::MHD2D:
                for(size_t i=0; i<num_xcells+4;i++){
                    u[i][num_ycells+2]=u[i][num_ycells+1];
                    u[i][num_ycells+3]=u[i][num_ycells];

                    //reflect y-velocity
                    u[i][num_ycells+2][2]*=-1;
                    u[i][num_ycells+3][2]*=-1;

                    //reflect B_y
                    u[i][num_ycells+2][6]*=-1;
                    u[i][num_ycells+3][6]*=-1;}break;
                case Config::DType::Euler2D:
                for(size_t i=0; i<num_xcells+4;i++){
                    u[i][num_ycells+2]=u[i][num_ycells+1];
                    u[i][num_ycells+3]=u[i][num_ycells];

                    //reflect y-velocity
                    u[i][num_ycells+2][2]*=-1;
                    u[i][num_ycells+3][2]*=-1;}break;
        break;}
    }
}