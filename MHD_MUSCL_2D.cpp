#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <filesystem>  // C++17
#include <initializer_list>

typedef std::array< double ,9> StateVector;

StateVector primitiveToConservative(const StateVector& prim,double gamma){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    StateVector conserved;
    conserved[0]= prim[0];
    conserved[1]=prim[0]*prim[1];
    conserved[2]=prim[0]*prim[2];
    conserved[3]=prim[0]*prim[3];
    conserved[4]=prim[4]/((gamma-1)) +0.5*prim[0]*(prim[1]*prim[1]+prim[2]*prim[2]+prim[3]*prim[3])+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]);
    conserved[5]=prim[5];
    conserved[6]=prim[6];
    conserved[7]=prim[7];
    conserved[8]=prim[8];

    return conserved;

}

StateVector conservativeToPrimitive(const StateVector& conserved,double gamma){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    StateVector prim;
    prim[0]= conserved[0];
    prim[1]=conserved[1]/conserved[0];
    prim[2]=conserved[2]/conserved[0];
    prim[3]=conserved[3]/conserved[0];
    prim[4] = (gamma-1)*(conserved[4]-0.5*(conserved[1]*conserved[1]+conserved[2]*conserved[2]+conserved[3]*conserved[3])/conserved[0]-0.5*(conserved[5]*conserved[5]+conserved[6]*conserved[6]+conserved[7]*conserved[7]));
    prim[5]=conserved[5];
    prim[6]=conserved[6];
    prim[7]=conserved[7];
    prim[8]=conserved[8];

    return prim;

}
   
void save_to_file(std::vector<std::vector<StateVector > >& u,std::vector<double> x, double num_xcells,std::vector<double> y, double num_ycells, double gamma, std::string dir,std::string filename){
    //Convert to primitive variables
    std::vector<std::vector<StateVector > > output(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    for(size_t i=0;i<num_xcells+4;i++){
        for(size_t j=0;j<num_ycells+4;j++){
        output[i][j]=conservativeToPrimitive(u[i][j],gamma);
    }}

    //Output the data
    std::filesystem::create_directories(dir);
    
    filename=dir+filename;

    std::ofstream outfile(filename);

    for(size_t j=2; j<num_ycells+2;j++){
        for(size_t i=2; i<num_xcells+2;i++){
            outfile << x[i] << " " << y[j] << " " << output[i][j][0] << " "<< output[i][j][1]<< " "<< output[i][j][2]<< " "<< output[i][j][3]<< " "<< output[i][j][4]<< " "<< output[i][j][5]<< " "<< output[i][j][6]<< " "<< output[i][j][7]<< " "<< output[i][j][8]<<"\n";
}
outfile<<"\n"; //Blank line between y rows for correct gnuplot formatting


}
    std::string gnuplot_cmd = "gnuplot -e \"filename='" + filename + "'\" HLLC_plot_script.plt";
    std::system(gnuplot_cmd.c_str());
    std::cout << "Plot saved to " <<filename<< ".png" << std::endl;
    //std::filesystem::remove(filename);
}


void update_bcs(std::vector<std::vector<StateVector > >& u,double num_xcells,double num_ycells){
    // Periodic 2D
    u[0]=u[num_xcells];
    u[1]=u[num_xcells+1];
    u[num_xcells+2]= u[2];
    u[num_xcells+3]= u[3];
    for(size_t i=0;i<num_xcells+4;i++){
        //For every row set the four ghost cells
        u[i][0]=u[i][num_ycells];
        u[i][1]=u[i][num_ycells+1];
        u[i][num_ycells+2]=u[i][2];
        u[i][num_ycells+3] = u[i][3];

    // //Kelvin-Helmholtz Test (periodic left-right reflective top-bottom)
    
    // u[0]=u[num_xcells];
    // u[1]=u[num_xcells+1];
    // u[num_xcells+2]= u[2];
    // u[num_xcells+3]= u[3];

    // for (int i = 0; i < num_xcells + 2; ++i) {//reflect in u_y and B_y
    // u[i][0] = u[i][2];
    // u[i][1] = u[i][3]; // Bottom boundary
    // u[i][num_ycells + 2] = u[i][num_ycells];
    // u[i][num_ycells + 3] = u[i][num_ycells + 1]; // Top boundary
    // u[i][0][2] = -u[i][2][2];
    // u[i][0][6] = -u[i][2][6];
    // u[i][1][2] = -u[i][3][2];
    // u[i][1][6] = -u[i][3][6];
    // u[i][num_ycells + 2][2] = -u[i][num_ycells][2];
    // u[i][num_ycells + 2][6] = -u[i][num_ycells][6];
    // u[i][num_ycells + 3][2] = -u[i][num_ycells + 1][2];
    // u[i][num_ycells + 3][6] = -u[i][num_ycells + 1][6];
}


    } 




//Set the initial t=0 condition from a given xarray

std::vector<std::vector<StateVector > > set_u0(std::vector<double> x,std::vector<double> y,double gamma,double num_xcells,double num_ycells){
    StateVector prim;
    std::vector<std::vector<StateVector > > u0(x.size(),std::vector<StateVector >(y.size()));
    double pi = 4*std::atan(1);

    //Define intiial conditions here!
    for(size_t i=0; i<x.size();i++){
        for(size_t j=0;j<y.size();j++){
        
        //Shock Tube (1D tests)

        // //Move from large x to small x with field
        // if(x[i]<=400 && y[j]<=400){
        //     //bottom left corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[5]=0.75;
        //     prim[6]=-1.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }
        // else if(x[i]<=400 && y[j]>400){
        //     //top left corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[5]=0.75;
        //     prim[6]=-1.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }
        // else if(x[i]>400 && y[j]<400){
        //     //bottom right corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[5]=0.75;
        //     prim[6]=1.;
        //     prim[7]=0.;
        //     prim[8]=0.;
            
        // }else if(x[i]>400 && y[j]>400){
        //     //top right corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[5]=0.75;
        //     prim[6]=1.;
        //     prim[7]=0.;
        //     prim[8]=0.;
        // }
        
        
        ////Move from small y to large y no B field
        // if(x[i]<=400 && y[j]<=400){
        //     //bottom left corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[6]=0.;
        //     prim[5]=0.;
        //     prim[7]=0.;
        //     prim[8]=0.;
        // }
        // else if(x[i]<=400 && y[j]>400){
        //     //top left corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[6]=0.;
        //     prim[5]=0.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }
        // else if(x[i]>400 && y[j]<400){
        //     //bottom right corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[6]=0.;
        //     prim[5]=0.;
        //     prim[7]=0.;
        //     prim[8]=0.;
        // }else if(x[i]>400 && y[j]>400){
        //     //top right corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[6]=0.;
        //     prim[5]=0.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }


        // //Move from small y to large y with B field
        // if(x[i]<=400 && y[j]<=400){
        //     //bottom left corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[6]=0.75;
        //     prim[5]=1.;
        //     prim[7]=0.;
        //     prim[8]=0.;
        // }
        // else if(x[i]<=400 && y[j]>400){
        //     //top left corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[6]=0.75;
        //     prim[5]=-1.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }
        // else if(x[i]>400 && y[j]<400){
        //     //bottom right corner
        //     prim[0]=1.;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=1.;
        //     prim[6]=0.75;
        //     prim[5]=1.;
        //     prim[7]=0.;
        //     prim[8]=0.;
        // }else if(x[i]>400 && y[j]>400){
        //     //top right corner
        //     prim[0]=0.125;
        //     prim[1]=0.;
        //     prim[2]=0.;
        //     prim[3]=0.;
        //     prim[4]=0.1;
        //     prim[6]=0.75;
        //     prim[5]=-1.;
        //     prim[7]=0.0;
        //     prim[8]=0.;
        // }

        
        //Orszang-Tang vortex test
        prim[0]=gamma*gamma;
        prim[1]=-std::sin(2*pi*y[j]);
        prim[2]=std::sin(2*pi*x[i]);
        prim[3]=0.;
        prim[4]=gamma;
        prim[5]=-std::sin(2*pi*y[j]);
        prim[6]=std::sin(4*pi*x[i]);
        prim[7]=0;
        prim[8]=0;

        // //Orszang-Tang vortex test
        // prim[0]=gamma*gamma;
        // prim[1]=-std::sin(2*pi*y[j]);
        // prim[2]=std::sin(2*pi*x[i]);
        // prim[3]=0.;
        // prim[4]=gamma;
        // prim[5]=-std::sin(2*pi*y[j]);
        // prim[6]=std::sin(4*pi*x[i]);
        // prim[7]=0;
        // prim[8]=0;

        // //Kelvin-Helmholtz
        // prim[0] = 1.0; // Density
        // prim[1] = 0.5*std::tanh(20*y[j]); // Velocity
        // prim[2] = 0.01*std::sin(2.0*pi*x[i])*std::exp(-y[j]*y[j]/(0.01));
        // prim[3] = 0.0;
        // prim[4] = 1.0 / gamma; // pressure
        // prim[5] = 0.1*std::cos(pi/3.0); // magnetic field
        // prim[6] = 0.0;
        // prim[7] = 0.1*std::sin(pi/3.0);
        // prim[8] = 0;



        u0[i][j]=primitiveToConservative(prim,gamma);
        }}
    update_bcs(u0,num_xcells,num_ycells);
       
    return u0;
    }
    
StateVector MHD_xflux(const StateVector& u0,double gamma,double c_h){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    StateVector flux;
    StateVector prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[1];
    flux[1]=prim[0]*prim[1]*prim[1]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[5]*prim[5];
    flux[2]= prim[0]*prim[1]*prim[2]-prim[5]*prim[6];
    flux[3]= prim[0]*prim[1]*prim[3]-prim[5]*prim[7];
    flux[4]=(u0[4]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[1]-(prim[1]*prim[5]+prim[2]*prim[6]+prim[3]*prim[7])*prim[5];

    flux[5]=u0[8];
    flux[8]=u0[5]*c_h*c_h;
    // flux[5]=0;
    // flux[8]=0;
    
    
    
    flux[6]=prim[6]*prim[1]-prim[5]*prim[2];
    flux[7]=prim[7]*prim[1]-prim[5]*prim[3];

    
    return flux;

}

StateVector MHD_yflux(const StateVector& u0,double gamma,double c_h){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    StateVector flux;
    StateVector prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[2];
    flux[1]= prim[0]*prim[2]*prim[1]-prim[6]*prim[5];
    flux[2]=prim[0]*prim[2]*prim[2]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[6]*prim[6];
    flux[3]= prim[0]*prim[2]*prim[3]-prim[6]*prim[7];
    flux[4]=(u0[4]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[2]-(prim[1]*prim[5]+prim[2]*prim[6]+prim[3]*prim[7])*prim[6];
    flux[5]=-(prim[6]*prim[1]-prim[2]*prim[5]);
    flux[6]=u0[8];
    // flux[6]=0;
    flux[7]=-(prim[6]*prim[3]-prim[2]*prim[7]);
    
    flux[8]=u0[6]*c_h*c_h;
    // flux[8]=0;
    
    return flux;

}


//HLL METHOD
double get_cfx(double gamma, const StateVector& u){
    
    // u is the vector of parameters at lattice site
    StateVector prim =conservativeToPrimitive(u,gamma);

    
    double factor= (gamma*prim[4]+u[5]*u[5]+u[6]*u[6]+u[7]*u[7])/u[0];

    double c_fx=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*u[5]*u[5]/(u[0]*u[0]))));
    return c_fx;

}
double get_cfy(double gamma, const StateVector& u){
    
    // u is the vector of parameters at lattice site
    StateVector prim =conservativeToPrimitive(u,gamma);
    
    double factor= (gamma*prim[4]+u[5]*u[5]+u[6]*u[6]+u[7]*u[7])/u[0];
    double c_fy=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*u[6]*u[6]/(u[0]*u[0]))));
    return c_fy;
}
double get_cfz(double gamma, const StateVector& u){
    
    // u is the vector of parameters at lattice site
    StateVector prim =conservativeToPrimitive(u,gamma);
    
    double factor= (gamma*prim[4]+u[5]*u[5]+u[6]*u[6]+u[7]*u[7])/u[0];
    double c_fz=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*u[7]*u[7]/(u[0]*u[0]))));
    return c_fz;
}


std::tuple< double, double > get_S_LR_x(const StateVector& u_L,const StateVector& u_R, double gamma){
    StateVector prim_L=conservativeToPrimitive(u_L,gamma);
    StateVector prim_R=conservativeToPrimitive(u_R,gamma);
    double v_xL = prim_L[1];
    double v_xR = prim_R[1];
    double c_fL=get_cfx(gamma,u_L);
    double c_fR=get_cfx(gamma,u_R);
    double S_L = std::min(v_xL,v_xR)-std::max(c_fL,c_fR);
    double S_R = std::max(v_xL,v_xR)+std::max(c_fL,c_fR);


    return {S_L,S_R};

}
std::tuple< double, double > get_S_LR_y(const StateVector& u_L,const StateVector& u_R, double gamma){
    StateVector prim_L=conservativeToPrimitive(u_L,gamma);
    StateVector prim_R=conservativeToPrimitive(u_R,gamma);
    double v_yL = prim_L[2];
    double v_yR = prim_R[2];
    double c_fL=get_cfy(gamma,u_L);
    double c_fR=get_cfy(gamma,u_R);
    double S_L = std::min(v_yL,v_yR)-std::max(c_fL,c_fR);
    double S_R = std::max(v_yL,v_yR)+std::max(c_fL,c_fR);



    return {S_L,S_R};

}

StateVector get_u_HLL_x(const StateVector& u_L,const StateVector& u_R, double gamma,double c_h){
    //uHLL states are an intermediate for solving the 7-wave part of the MHD equations
    //Bx and psi are fixed at their tilde values at every step
    //In theory the for loops would only update 7 of the values not all 9
    //uHLL for Bx and psi would just also be the tilde values
    
    auto [S_L,S_R]= get_S_LR_x(u_L,u_R,gamma);

    StateVector f_L=MHD_xflux(u_L,gamma, c_h);
    StateVector f_R=MHD_xflux(u_R,gamma,c_h);
    StateVector u_HLL;

    for(size_t i=0;i<9;i++){
        //This performs an intermediate update for 7 parameters
        //For Bx and psi since L and R vals are identical this just returns psi_HLL=psi_L=psi_R=psi_tilde
        u_HLL[i]=1/(S_R-S_L)*(S_R*u_R[i]-S_L*u_L[i]+f_L[i]-f_R[i]);
        
    }

    if(S_L>=0){
        return u_L;
    }else if(S_L<0 && S_R>0){
        return u_HLL;
    } else{
        return u_R;
    }
}

StateVector get_u_HLL_y(const StateVector& u_L,const StateVector& u_R, double gamma,double c_h){
    auto [S_L,S_R]= get_S_LR_y(u_L,u_R,gamma);

    StateVector f_L=MHD_yflux(u_L,gamma, c_h);
    StateVector f_R=MHD_yflux(u_R,gamma, c_h);
    StateVector u_HLL;

    for(size_t i=0;i<9;i++){
        u_HLL[i]=1/(S_R-S_L)*(S_R*u_R[i]-S_L*u_L[i]+f_L[i]-f_R[i]);
    }

    if(S_L>=0){
        return u_L;
    }else if(S_L<0 && S_R>0){
        return u_HLL;
    }else{
        return u_R;
    }
}

StateVector get_u_star_x(const StateVector& u_L,const StateVector& u_R, double gamma,bool L,double B_xtilde,double psi_tilde, double S_L, double S_R,double c_h){
    //L bool dictates whether you are calculating u_starL or u_starR 
    //Everything here should have BxL=BxR= Bxtilde apart from the S_L and S_R calcs
    
    StateVector u_HLL=get_u_HLL_x(u_L,u_R,gamma,c_h);
    StateVector prim_HLL=conservativeToPrimitive(u_HLL,gamma);
    
    
    StateVector prim_L= conservativeToPrimitive(u_L,gamma);
    StateVector prim_R= conservativeToPrimitive(u_R,gamma);
    
    double S; // Will be used as S_L for L case and S_R for R_case
    StateVector prim; // Will be used as L for L case and R for R_case
    StateVector u;

    if(L){
        S=S_L;
        prim=prim_L;
        u=u_L;

    }else{
        S=S_R;
        prim=prim_R;
        u=u_R;

    }
    //universal quantity used for both L and R
    double q_star= (prim_R[0]*prim_R[1]*(S_R-prim_R[1])-prim_L[0]*prim_L[1]*(S_L-prim_L[1])+prim_L[4]+0.5*(prim_L[5]*prim_L[5]+prim_L[6]*prim_L[6]+prim_L[7]*prim_L[7])-prim_R[4]-0.5*(prim_R[5]*prim_R[5]+prim_R[6]*prim_R[6]+prim_R[7]*prim_R[7])-prim_L[5]*prim_L[5]+prim_R[5]*prim_R[5])/(prim_R[0]*(S_R-prim_R[1])-prim_L[0]*(S_L-prim_L[1]));
    
    StateVector u_star;

    // Bx and psi decouple and can be handled exactly with a two wave solution with speeds, +-c_h 
    // Bx and psi are always in the central region of the two wave system
    //These solutions are given by the tilde values

    u_star[5]=u_HLL[5];
    u_star[8]=psi_tilde;


    u_star[6]=u_HLL[6];
    u_star[7]=u_HLL[7];

    //Calc. the momentum and density terms
    u_star[0]=prim[0]*(S-prim[1])/(S-q_star);
    u_star[1]=u_star[0]*q_star;
    u_star[2]=u_star[0]*prim[2]-(u_star[5]*u_star[6]-prim[5]*prim[6])/(S-q_star);
    u_star[3]=u_star[0]*prim[3]-(u_star[5]*u_star[7]-prim[5]*prim[7])/(S-q_star);


    //Calc. the E_star term
    double p_star= prim[0]*(S-prim[1])*(q_star-prim[1])+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[5]*prim[5]+u_star[5]*u_star[5];
    
    u_star[4]=u[4]*(S-prim[1])/(S-q_star)+((p_star*q_star-(prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[1])-(u_star[5]*(prim_HLL[5]*prim_HLL[1]+prim_HLL[6]*prim_HLL[2]+prim_HLL[7]*prim_HLL[3])-prim[5]*(prim[5]*prim[1]+prim[6]*prim[2]+prim[7]*prim[3])))/(S-q_star);

    
    
    return u_star;
}

StateVector get_u_star_y(const StateVector& u_L,const StateVector& u_R, double gamma,bool L,double B_ytilde,double psi_tilde, double S_L, double S_R,double c_h){
    //L bool dictates whether you are calculating u_starL or u_starR 
    
    StateVector u_HLL=get_u_HLL_y(u_L,u_R,gamma,c_h);
    StateVector prim_HLL=conservativeToPrimitive(u_HLL,gamma);
   
    StateVector prim_L= conservativeToPrimitive(u_L,gamma);
    StateVector prim_R= conservativeToPrimitive(u_R,gamma);
    
    double S; // Will be used as S_L for L case and S_R for R_case
    StateVector prim; // Will be used as L for L case and R for R_case
    StateVector u;

    if(L){
        S=S_L;
        prim=prim_L;
        u=u_L;

    }else{
        S=S_R;
        prim=prim_R;
        u=u_R;

    }
    //universal quantity used for both L and R
    double q_star= (prim_R[0]*prim_R[2]*(S_R-prim_R[2])-prim_L[0]*prim_L[2]*(S_L-prim_L[2])+prim_L[4]+0.5*(prim_L[5]*prim_L[5]+prim_L[6]*prim_L[6]+prim_L[7]*prim_L[7])-prim_R[4]-0.5*(prim_R[5]*prim_R[5]+prim_R[6]*prim_R[6]+prim_R[7]*prim_R[7])-prim_L[6]*prim_L[6]+prim_R[6]*prim_R[6])/(prim_R[0]*(S_R-prim_R[2])-prim_L[0]*(S_L-prim_L[2]));
    
    StateVector u_star;
    

    u_star[6]=u_HLL[6];
    u_star[8]=psi_tilde;
    
    u_star[5]=u_HLL[5];
    u_star[7]=u_HLL[7];

    //Calc. the momentum and density terms
    u_star[0]=prim[0]*(S-prim[2])/(S-q_star);
    u_star[1]=u_star[0]*prim[1]-(u_star[5]*u_star[6]-prim[5]*prim[6])/(S-q_star);
    u_star[2]=u_star[0]*q_star;
    u_star[3]=u_star[0]*prim[3]-(u_star[6]*u_star[7]-prim[6]*prim[7])/(S-q_star);

    //Calc. the E_star term
    double p_star= prim[0]*(S-prim[2])*(q_star-prim[2])+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[6]*prim[6]+u_star[6]*u_star[6];
    
    u_star[4]=u[4]*(S-prim[2])/(S-q_star)+((p_star*q_star-(prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[2])-(u_star[6]*(prim_HLL[5]*prim_HLL[1]+prim_HLL[6]*prim_HLL[2]+prim_HLL[7]*prim_HLL[3])-prim[6]*(prim[5]*prim[1]+prim[6]*prim[2]+prim[7]*prim[3])))/(S-q_star);

    return u_star;
}


double get_c_h(const std::vector<std::vector<StateVector > >& u,double num_xcells,double num_ycells,double gamma){

    double max_ch=0.;
    StateVector prim;
    
    for(size_t i=2; i<num_xcells+2;i++){
        for(size_t j=2; j<num_ycells+2;j++){
        prim=conservativeToPrimitive(u[i][j],gamma);

        double c_fx=get_cfx(gamma,u[i][j]);
        double c_fy = get_cfy(gamma,u[i][j]);
        double c_fz = get_cfz(gamma,u[i][j]);

        double new_max=std::max({std::abs(prim[1])+c_fx,std::abs(prim[2])+c_fy,std::abs(prim[3])+c_fz});

        if(new_max> max_ch){
            max_ch=new_max;
        }   
    }
}

return max_ch;
}

double computeTimeStep(double dx,double dy, double C,double c_h){
 //Making an adjustment here
return C* std::min(dx,dy)/(c_h);
}


std::tuple< double,double > get_xtilde_vals(const StateVector& uL, const StateVector& uR,double gamma,double c_h){
    
    double Bx_tilde = 0.5*(uL[5]+uR[5])-0.5*1./c_h*(uR[8]-uL[8]);
    double psi_tilde= 0.5*(uL[8]+uR[8])-c_h/2.*(uR[5]-uL[5]);
    //These are the decoupled solutions to the 2D Riemann problem
    return {Bx_tilde,psi_tilde};

}

std::tuple< double,double > get_ytilde_vals(const StateVector& uL, const StateVector& uR,double gamma,double c_h){
    
    double By_tilde = 0.5*(uL[6]+uR[6])-0.5*1./c_h*(uR[8]-uL[8]);
    double psi_tilde= 0.5*(uL[8]+uR[8])-c_h/2.*(uR[6]-uL[6]);
    return {By_tilde,psi_tilde};

}

StateVector get_HLLC_flux_x(StateVector u_L,StateVector u_R, double gamma,double c_h){
    //gets flux between L and R position
    //Solve decoupled 2 wave system first
    auto[B_tilde,psi_tilde]=get_xtilde_vals(u_L,u_R,gamma,c_h);

    //Now use B_tilde and x_tilde values as fixed parameters in solving the remaining 7 wave system
    
    u_L[5]=B_tilde;
    u_R[5]=B_tilde;
    u_L[8] = psi_tilde;
	u_R[8] = psi_tilde;
    
    auto [S_L,S_R]= get_S_LR_x(u_L,u_R,gamma);

    
    
    
    StateVector f_L=MHD_xflux(u_L,gamma,c_h);
    StateVector f_R=MHD_xflux(u_R,gamma,c_h);
    
    StateVector u_starL=get_u_star_x(u_L,u_R,gamma,true,B_tilde,psi_tilde, S_L,S_R,c_h);
    StateVector u_starR=get_u_star_x(u_L,u_R,gamma,false,B_tilde,psi_tilde,S_L,S_R,c_h);

    double q_star= u_starL[1]/u_starL[0];
    StateVector flux_out;

    //Deciding what region/flux is being used 
    if(S_L>=0){
        flux_out=f_L;

    }else if(S_L<0 && q_star>=0){
        for(size_t i=0;i<9;i++){
        flux_out[i]=f_L[i]+S_L*(u_starL[i]-u_L[i]);
    }}else if(q_star<=0 && 0<=S_R){
        for(size_t i=0;i<9;i++){
        flux_out[i]=f_R[i]+S_R*(u_starR[i]-u_R[i]);

    }}else {
        flux_out=f_R; 
    }
    
    

    // flux_out[5]=0;
    // flux_out[8]=0;

    

    return flux_out;


    // std::vector<double> fluxHLL(9);
    // for(size_t i=0;i<9;i++){
    //     fluxHLL[i]=(S_R*f_L[i]-S_L*f_R[i] +S_L*S_R*(u_R[i]-u_L[i]))/(S_R-S_L);
    // }
    // fluxHLL[5]=psi_tilde;
    // fluxHLL[8]=B_tilde*c_h*c_h;
    // return fluxHLL;
}

StateVector get_HLLC_flux_y(StateVector u_L,StateVector u_R, double gamma,double c_h){
    //gets flux between L and R position
    auto[B_tilde,psi_tilde]=get_ytilde_vals(u_L,u_R,gamma,c_h);
    auto [S_L,S_R]= get_S_LR_y(u_L,u_R,gamma);

    u_L[6]=B_tilde;
    u_R[6]=B_tilde;
    u_L[8] = psi_tilde;
	u_R[8] = psi_tilde;


    StateVector f_L=MHD_yflux(u_L,gamma,c_h);
    StateVector f_R=MHD_yflux(u_R,gamma,c_h);
    StateVector u_starL=get_u_star_y(u_L,u_R,gamma,true,B_tilde,psi_tilde,S_L,S_R,c_h);
    StateVector u_starR=get_u_star_y(u_L,u_R,gamma,false,B_tilde,psi_tilde, S_L,S_R,c_h);

    double q_star= u_starL[2]/u_starL[0];
    StateVector flux_out;

    //Deciding what region/flux is being used
    if(S_L>=0){
        flux_out=f_L;
    }else if(S_L<0 && q_star>=0){
        for(size_t i=0;i<9;i++){
        flux_out[i]=f_L[i]+S_L*(u_starL[i]-u_L[i]);
    }}else if(q_star<=0 && 0<=S_R){
        for(size_t i=0;i<9;i++){
        flux_out[i]=f_R[i]+S_R*(u_starR[i]-u_R[i]);
    }}else {
        flux_out=f_R; 
    }




    // flux_out[6]=0;
    // flux_out[8]=0;

    return flux_out;

    // std::vector<double> fluxHLL(9);
    // for(size_t i=0;i<9;i++){
    //     fluxHLL[i]=(S_R*f_L[i]-S_L*f_R[i] +S_L*S_R*(u_R[i]-u_L[i]))/(S_R-S_L);
    // }
    // fluxHLL[6]=psi_tilde;
    // fluxHLL[8]=B_tilde*c_h*c_h;
    // return fluxHLL;
}


//Adding Slope Limiting

StateVector getDelta_plus(const StateVector& u_0, const StateVector& u_plus){

    StateVector delta_plus;
    for(size_t i=0;i<9;i++){
        delta_plus[i]=u_plus[i]-u_0[i];
    }
    return delta_plus;
}
StateVector getDelta_minus(const StateVector& u_minus, const StateVector& u_0){
    StateVector delta_minus;
    for(size_t i=0;i<9;i++){
        delta_minus[i]=u_0[i]-u_minus[i];
    }
    return delta_minus;
}
StateVector getDelta(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w=0){
StateVector delta_minus = getDelta_minus(u_minus,u_0);
StateVector delta_plus = getDelta_plus(u_0,u_plus);
StateVector delta;
    for(size_t i=0;i<9;i++){
        delta[i]=0.5*(1+w)*delta_minus[i] + 0.5*(1-w)*delta_plus[i];
    }
    return delta;
}

StateVector get_r(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double gamma){
    StateVector delta_minus = getDelta_minus(u_minus,u_0);
    StateVector delta_plus = getDelta_plus(u_0,u_plus);
    StateVector r;
    
    
    double min_r=2;   
    
    for(size_t i=0;i<9;i++){
    
        r[i]=delta_minus[i]/(delta_plus[i]+1e-12);
        // if(r[i]<min_r && r[i] != 0.){
        //     min_r=r[i];
        // }
        
    }
    
    // for(size_t i=0;i<9;i++){
    //     r[i]=min_r;
    //     }
        return r;
}
StateVector getXi_VL(StateVector r){
    //Van Leer slope limiting
    StateVector xi;

    for(size_t i=0;i<9;i++){

        if(r[i]<=0){
            xi[i]=0;
        } else{
            xi[i]=std::min(2*r[i]/(1+r[i]),2/(1+r[i]));
        } 
    }
    return xi;
}

std::tuple< StateVector,StateVector  > calc_ubar(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w,double gamma){
    //returns {ubar_L,ubar_R}
    StateVector prim_minus = conservativeToPrimitive(u_minus,gamma);
    StateVector prim_0 = conservativeToPrimitive(u_0,gamma);
    StateVector prim_plus = conservativeToPrimitive(u_plus,gamma);
    
    StateVector ubar_L;
    StateVector ubar_R;

    StateVector delta=getDelta(u_minus,u_0,u_plus,w);
    StateVector r= get_r(prim_minus,prim_0,prim_plus,gamma);
    StateVector xi= getXi_VL(r);
    
    for(size_t j=0;j<9;j++){
        //Looping over variables
            ubar_L[j]=u_0[j]-0.5*xi[j]*delta[j];
            ubar_R[j]=u_0[j]+0.5*xi[j]*delta[j];
        }  

        return {ubar_L,ubar_R};
    }
std::tuple< StateVector,StateVector > calc_ubarx_plus(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w,double dt, double dx,double gamma,double c_h){
    //returns {ubar_L^(n+1/2),ubar_R^(n+1/2) }
    auto [ubar_L, ubar_R]= calc_ubar(u_minus,u_0,u_plus,w,gamma);
    StateVector ubar_L_plus;
    StateVector ubar_R_plus;
    StateVector f_R= MHD_xflux(ubar_R,gamma,c_h);
    StateVector f_L= MHD_xflux(ubar_L,gamma,c_h);
    for(size_t j=0;j<9;j++){//Looping over variables
            ubar_L_plus[j]= ubar_L[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            ubar_R_plus[j]= ubar_R[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            
        }  
        return {ubar_L_plus,ubar_R_plus};

    }
std::tuple< StateVector,StateVector > calc_ubary_plus(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w,double dt, double dx,double gamma,double c_h){
    //returns {ubar_L^(n+1/2),ubar_R^(n+1/2) }
    auto [ubar_L, ubar_R]= calc_ubar(u_minus,u_0,u_plus,w,gamma);
    StateVector ubar_L_plus;
    StateVector ubar_R_plus;
    StateVector f_R= MHD_yflux(ubar_R,gamma,c_h);
    StateVector f_L= MHD_yflux(ubar_L,gamma,c_h);
    for(size_t j=0;j<9;j++){//Looping over variables
            ubar_L_plus[j]= ubar_L[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            ubar_R_plus[j]= ubar_R[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            
        }  
        return {ubar_L_plus,ubar_R_plus};

    }

std::tuple< std::vector<std::vector<StateVector > >,std::vector<std::vector<StateVector > > > get_ubar_plus_xarrays(const std::vector<std::vector<StateVector > >& u, double w,double dt, double dx,double gamma,double num_xcells,double num_ycells,double c_h){
    std::vector<std::vector<StateVector > > ubar_Rplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > ubar_Lplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    
    for(size_t i=1;i<num_xcells+3;i++){
        for(size_t j=1;j<num_ycells+3;j++){
        auto [L, R] = calc_ubarx_plus(u[i-1][j],u[i][j],u[i+1][j],w,dt,dx,gamma,c_h);
        ubar_Lplus[i][j]=L;
        ubar_Rplus[i][j]=R;  
    }}


    update_bcs(ubar_Lplus,num_xcells,num_ycells);
    update_bcs(ubar_Rplus,num_xcells,num_ycells);
    
    return {ubar_Lplus,ubar_Rplus};
}
std::tuple< std::vector<std::vector<StateVector > >,std::vector<std::vector<StateVector > > > get_ubar_plus_yarrays(const std::vector<std::vector<StateVector > >& u, double w,double dt, double dx,double gamma,double num_xcells,double num_ycells,double c_h){
    std::vector<std::vector<StateVector > > ubar_Rplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > ubar_Lplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    
    for(size_t i=1;i<num_xcells+3;i++){
        for(size_t j=1;j<num_ycells+3;j++){
        auto [L, R] = calc_ubary_plus(u[i][j-1],u[i][j],u[i][j+1],w,dt,dx,gamma,c_h);
        ubar_Lplus[i][j]=L;
        ubar_Rplus[i][j]=R;  
    }}

    
    update_bcs(ubar_Lplus,num_xcells,num_ycells);
    update_bcs(ubar_Rplus,num_xcells,num_ycells);
    
    return {ubar_Lplus,ubar_Rplus};
}


int main(void){
    
    //Orszang-Tang vortex test
    double C=0.8;
    double tEnd=1.1;
    double tStart=0;
    double x0=0;
    double xf=1;
    double y0=0;
    double yf=1;
    double num_xcells=256;
    double num_ycells=256;
    double gamma=5./3. ;
    double w=0;
    int save_interval=10;

    // // Shock Tube Test
    // double C=0.8;
    // double tEnd=80;
    // double tStart=0;
    // double x0=0;
    // double xf=800;
    // double y0=0;
    // double yf=800;
    // double num_xcells=250;
    // double num_ycells=250;
    // double gamma=2. ;
    // double w=0;
    // int save_interval=10;


    // //Kelvin-Helmholtz test
    // double C=0.8;
    // double tEnd=20;
    // double tStart=0;
    // double x0=0;
    // double xf=1;
    // double y0=-1;
    // double yf=1;
    // double num_xcells=128;
    // double num_ycells=256;
    // double gamma=5./3. ;
    // double w=0;
    // int save_interval=50;


    


    double dx=(xf-x0)/num_xcells;
    double dy=(yf-y0)/num_ycells;

    //Initialise x and y vectors which mark positions of cell centres
    std::vector<double> x(num_xcells+4);
    for(int i=0;i<x.size();i++){
        x[i]=x0+(i-0.5)*dx;

    };
    
    std::vector<double> y(num_ycells+4);
    for(int i=0;i<y.size();i++){
        y[i]=y0+(i-0.5)*dy;

    };



    //Initialise u(x) vector
    std::vector<std::vector<StateVector > > u(num_xcells+4, std::vector<StateVector >(num_ycells+4));

    // store flux: flux[i] is flux at x(i+1/2)
    std::vector<std::vector<StateVector > > flux_x(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > flux_y(num_xcells+4, std::vector<StateVector >(num_ycells+4));

    u=set_u0(x,y,gamma,num_xcells,num_ycells);

    double t=tStart;
    
    int count=0; //counts number of update steps performed

    do {


        if(count%save_interval==0){
            std::string name= "Plot_" + std::to_string(t);
             //Output the data
            std::string dir = "Test/";
            save_to_file(u,x,num_xcells,y,num_ycells,gamma,dir,name);
        }

        update_bcs(u,num_xcells,num_ycells);
        
        double c_h= get_c_h(u,num_xcells,num_ycells,gamma); //Max across entirety of u, should be updated every time u is actually updated 
        double c_p_squared=c_h*0.18;
        double dt=computeTimeStep(dx,dy,C,c_h);

        t+=dt;
        count++;
        
        std::cout<<"Starting iteration: "<<count<<std::endl;

        //Do x updates first
        
        //Do half time step source update for all psi
        
        for(size_t i=0; i<num_xcells+4; i++){
            for(size_t j=0; j<num_ycells+4; j++){

                u[i][j][8]*=std::exp((-dt/2)*(c_h*c_h)/c_p_squared);
                }}
        
        //First update the u_bar arrays to correspond to (x) flux for slope limiting
        auto [ubar_L_plusx,ubar_R_plusx] = get_ubar_plus_xarrays(u,w,dt,dx,gamma,num_xcells,num_ycells,c_h);

         // update x flux array
        

        for(size_t i=1; i<num_xcells+2; i++){
            for(size_t j=1; j<num_ycells+2; j++){

                flux_x[i][j]=get_HLLC_flux_x(ubar_R_plusx[i][j],ubar_L_plusx[i+1][j],gamma,c_h);
            }}

        //Update all real cells using x flux
        for(size_t i=2; i<num_xcells+2; i++){
            for(size_t j=2; j<num_ycells+2; j++){
                for(size_t k=0;k<9;k++){
                u[i][j][k]=u[i][j][k] - (dt/dx) * (flux_x[i][j][k]-flux_x[i-1][j][k]);
                }
            }
        }
  
        update_bcs(u,num_xcells,num_ycells);

        //Now perform y updates with the same timestep

        auto [ubar_L_plusy,ubar_R_plusy] = get_ubar_plus_yarrays(u,w,dt,dy,gamma,num_xcells,num_ycells,c_h);

    
         // update y flux array
        
        for(size_t i=1; i<num_xcells+2; i++){
            for(size_t j=1; j<num_ycells+2; j++){

                flux_y[i][j]=get_HLLC_flux_y(ubar_R_plusy[i][j],ubar_L_plusy[i][j+1],gamma,c_h);
            }}

        //Update all real cells using y flux
        for(size_t i=2; i<num_xcells+2; i++){
            for(size_t j=2; j<num_ycells+2; j++){
                for(size_t k=0;k<9;k++){
                u[i][j][k]=u[i][j][k] - (dt/dy) * (flux_y[i][j][k]-flux_y[i][j-1][k]);
                }
                
                //Now do the second partial source update for psi
               
                u[i][j][8]*=std::exp(-dt*(c_h*c_h)/c_p_squared);
         
            }
        }
        
        update_bcs(u,num_xcells,num_ycells);

      
        

   

    }while (t<tEnd);

}

    



