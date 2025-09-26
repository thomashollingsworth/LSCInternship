#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <filesystem>  // C++17
#include <tuple>

typedef std::array< double ,4> StateVector; //Stores the 4 conserved variables defined at each grid point

StateVector primitiveToConservative(const StateVector& prim,double gamma){
    //Convert a 1D vector of density,vx,vy,p to the conserved quantities density,momentum_x,momentum_y,energy
    StateVector conserved;
    conserved[0]= prim[0];
    conserved[1]=prim[0]*prim[1];
    conserved[2]=prim[0]*prim[2];
    conserved[3] = prim[3]/((gamma-1)) +0.5*prim[0]*(prim[1]*prim[1]+prim[2]*prim[2]);

    return conserved;

}

StateVector conservativeToPrimtive(const StateVector& conserved,double gamma){
    //Convert a 1D vector of density,mom_x,mom_y, energy to density,vx,vy,pressure
    StateVector prim;
    prim[0]= conserved[0];
    prim[1]=conserved[1]/conserved[0];
    prim[2]=conserved[2]/conserved[0];
    prim[3] = (gamma-1)*(conserved[3]-0.5*(conserved[1]*conserved[1]+conserved[2]*conserved[2])/conserved[0]);

    return prim;

}
    
void save_to_file(const std::vector<std::vector<StateVector > >& u,const std::vector<double >& x, double num_xcells,const std::vector<double>& y, double num_ycells, double gamma, std::string filename){
    //Convert to primitive variables
    std::vector<std::vector<StateVector > > output(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    for(size_t i=0;i<num_xcells+4;i++){
        for(size_t j=0;j<num_ycells+4;j++){
        output[i][j]=conservativeToPrimtive(u[i][j],gamma);
    }}

    //Output the data
    std::string dir = "2D_output_data/";
    std::filesystem::create_directories(dir);
    
    filename=dir+filename;

    std::ofstream outfile(filename);

    for(size_t j=2; j<num_ycells+2;j++){
        for(size_t i=2; i<num_xcells+2;i++){
            outfile << x[i] << " " << y[j] << " " << output[i][j][0] << " "<< output[i][j][1]<< " "<< output[i][j][2]<< " "<< output[i][j][3]<< "\n";
}
outfile<<"\n"; //Blank line between y rows for correct gnuplot formatting


}
    std::string gnuplot_cmd = "gnuplot -e \"filename='" + filename + "'\" plot_script.plt";
    std::system(gnuplot_cmd.c_str());
    std::cout << "Plot saved to " <<filename<< ".png" << std::endl;
    //std::filesystem::remove(filename);
}

void update_bcs_trans(std::vector<std::vector<StateVector > >& u){
    u[0]=u[2];
    u[1]=u[2];
    u[u.size()-1]= u[u.size()-3];
    u[u.size()-2]= u[u.size()-3];
    for(size_t i=0;i<u.size();i++){
        //For every row set the four ghost cells
        u[i][0]=u[i][2];
        u[i][1]=u[i][2];
        u[i][u[i].size()-1]=u[i][u[i].size()-3];
        u[i][u[i].size()-2] = u[i][u[i].size()-3];


    } 

}

double computeTimeStep(const std::vector<std::vector<StateVector > >& u,double dx,double dy, double C,double gamma,double num_xcells,double num_ycells){
    StateVector prim0;
    prim0=conservativeToPrimtive(u[0][0],gamma);
    double max=std::abs(std::sqrt(prim0[1]*prim0[1]+prim0[2]*prim0[2]))+std::sqrt(gamma*prim0[3]/prim0[0]);
    
    StateVector prim;
    //Calculate max wave speed for all real cells
    for(size_t i=2; i<num_xcells+2;i++){
        for(size_t j=2; j<num_ycells+2;j++) 
        prim=conservativeToPrimtive(u[i][j],gamma);
        double new_max=std::abs(std::sqrt(prim[1]*prim[1]+prim[2]*prim[2]))+std::sqrt(gamma*prim[3]/prim[0]);
        if(new_max> max){
            max=new_max;
        }
        
    }
 return C*std::min(dx,dy)/max;
}

//Set the initial t=0 condition from a given xarray

std::vector<std::vector<StateVector > > set_u0(const std::vector<double>& x,const std::vector<double>& y,double gamma){
    StateVector prim;
    std::vector<std::vector<StateVector > > u0(x.size(),std::vector<StateVector >(y.size()));

    //Define intiial conditions here!
    for(size_t i=0; i<x.size();i++){
        for(size_t j=0;j<y.size();j++){
            //Set initial condition in primitive coords
        if(x[i]<=0.5 && y[j]<=0.5){
            //bottom left corner
            prim[0]=1.0;
            prim[1]=0.1;
            prim[2]=0.1;
            prim[3]=1.0;
        }
        else if(x[i]<=0.5 && y[j]>0.5){
            //top left corner
            prim[0]=0.125;
            prim[1]=0.0;
            prim[2]=1.0;
            prim[3]=1.0;
        }else if(x[i]>0.5 && y[j]<0.5){
            //bottom right corner
            prim[0]=1.0;
            prim[1]=-2.0;
            prim[2]=1.0;
            prim[3]=0.01;
        }else if(x[i]>0.5 && y[j]>0.5){
            //top right corner
            prim[0]=2.0;
            prim[1]=2.0;
            prim[2]=0.0;
            prim[3]=0.125;
        }
        u0[i][j]=primitiveToConservative(prim,gamma);
        }}
    return u0;
    }
    
StateVector euler_flux(const StateVector& u0,double gamma){
    StateVector flux;
    StateVector prim=conservativeToPrimtive(u0,gamma);
    flux[0]=prim[0]*prim[1];
    flux[1]=prim[0]*prim[1]*prim[1] +prim[3];
    flux[2]=prim[0]*prim[1]*prim[2];
    flux[3]= (u0[3]+prim[3])*prim[1];
    return flux;

}

StateVector getFluxRI(const StateVector& u_0,const StateVector& u_1,const double& dx,const double& dt,double gamma){
    //returns the flux f_(i+1/2) using Richtymer method
    StateVector f_plus=euler_flux(u_1,gamma);
    StateVector f_minus=euler_flux(u_0,gamma);
    StateVector u_plus_half;
    for(size_t i=0;i<4;i++){
        u_plus_half[i]=0.5*(u_0[i]+u_1[i]) -0.5*dt/dx*(f_plus[i]-f_minus[i]);

    }
  
    return euler_flux(u_plus_half,gamma);


}

StateVector getFluxLF(const StateVector& u_0,const StateVector& u_1,const double& dx,const double& dt,double gamma){
    //Lex Friedrichs flux
    //returns the flux f_(i+1/2)
    StateVector f_plus=euler_flux(u_1,gamma);
    StateVector f_minus=euler_flux(u_0,gamma);
    StateVector flux;
    for(size_t i=0;i<4;i++){
        flux[i]=0.5*dx/dt*(u_0[i]-u_1[i])+0.5*(f_minus[i]+f_plus[i]);
    }

    return flux;

}

StateVector getFluxFORCE(const StateVector& u_0, const StateVector& u_1,const double& dx,const double& dt,double gamma){
    //Force flux (avg. of RI and LF flux)
    StateVector flux;
    StateVector fluxLF=getFluxLF(u_0,u_1,dx,dt,gamma);
    StateVector fluxRI=getFluxRI(u_0,u_1,dx,dt,gamma);
    for(size_t i=0;i<4;i++){
        flux[i]=0.5*(fluxLF[i]+fluxRI[i]);
    }
    return flux;
}


//USED FOR SLOPE LIMITING METHOD
StateVector getDelta_plus(const StateVector& u_0, const StateVector& u_plus){

    StateVector delta_plus;
    for(size_t i=0;i<4;i++){
        delta_plus[i]=u_plus[i]-u_0[i];
    }
    return delta_plus;
}
StateVector getDelta_minus(const StateVector& u_minus, const StateVector& u_0){
    StateVector delta_minus;
    for(size_t i=0;i<4;i++){
        delta_minus[i]=u_0[i]-u_minus[i];
    }
    return delta_minus;
}
StateVector getDelta(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w){
StateVector delta_minus = getDelta_minus(u_minus,u_0);
StateVector delta_plus = getDelta_plus(u_0,u_plus);
StateVector delta;
    for(size_t i=0;i<4;i++){
        delta[i]=0.5*(1+w)*delta_minus[i] + 0.5*(1-w)*delta_plus[i];
    }
    return delta;
}

StateVector get_r(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus){
    StateVector delta_minus = getDelta_minus(u_minus,u_0);
    StateVector delta_plus = getDelta_plus(u_0,u_plus);
    StateVector r;
        for(size_t i=0;i<4;i++){
        r[i]=delta_minus[i]/delta_plus[i];
        }
        return r;
}
StateVector getXi_Minbee(const StateVector& r){
    StateVector xi;
    for(size_t i=0;i<4;i++){
        if(r[i]<=0){
            xi[i]=0;
        } else if(r[i]<=1){
            xi[i]=r[i];

        } else{
            double xi_R = 2/(1+r[i]);
            xi[i]=std::min(1.0,xi_R);
        }
    }
    return xi;

}


std::tuple< StateVector,StateVector  > calc_ubar(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w){
    //returns {ubar_L,ubar_R}
    StateVector ubar_L;
    StateVector ubar_R;
    StateVector delta=getDelta(u_minus,u_0,u_plus,w);
    StateVector r= get_r(u_minus,u_0,u_plus);
    StateVector xi= getXi_Minbee(r);
    for(size_t j=0;j<4;j++){
        //Looping over variables
            ubar_L[j]=u_0[j]-0.5*xi[j]*delta[j];
            ubar_R[j]=u_0[j]+0.5*xi[j]*delta[j];
        }  
        return {ubar_L,ubar_R};
    }
std::tuple< StateVector,StateVector > calc_ubar_plus(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w,double dt, double dx,double gamma){
    //returns {ubar_L^(n+1/2),ubar_R^(n+1/2) }
    auto [ubar_L, ubar_R]= calc_ubar(u_minus,u_0,u_plus,w);
    StateVector ubar_L_plus;
    StateVector ubar_R_plus;
    StateVector f_R= euler_flux(ubar_R,gamma);
    StateVector f_L= euler_flux(ubar_L,gamma);
    for(size_t j=0;j<4;j++){//Looping over variables
            ubar_L_plus[j]= ubar_L[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            ubar_R_plus[j]= ubar_R[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            
        }  
        return {ubar_L_plus,ubar_R_plus};

    }
std::tuple< std::vector<std::vector<StateVector > >,std::vector<std::vector<StateVector > > > get_ubar_plus_xarrays(const std::vector<std::vector<StateVector > >& u, double w,double dt, double dx,double gamma,double num_xcells,double num_ycells){
    std::vector<std::vector<StateVector > > ubar_Rplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > ubar_Lplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    
    for(size_t i=1;i<num_xcells+3;i++){
        for(size_t j=1;j<num_ycells+3;j++){
        auto [L, R] = calc_ubar_plus(u[i-1][j],u[i][j],u[i+1][j],w,dt,dx,gamma);
        ubar_Lplus[i][j]=L;
        ubar_Rplus[i][j]=R;
        }
    }
    return {ubar_Lplus,ubar_Rplus};
}
std::tuple< std::vector<std::vector<StateVector > >,std::vector<std::vector<StateVector > > > get_ubar_plus_yarrays(const std::vector<std::vector<StateVector > >& u, double w,double dt, double dy,double gamma,double num_xcells,double num_ycells){
    std::vector<std::vector<StateVector > > ubar_Rplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > ubar_Lplus(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    
    for(size_t i=1;i<num_xcells+3;i++){
        for(size_t j=1;j<num_ycells+3;j++){
        auto [L, R] = calc_ubar_plus(u[i][j-1],u[i][j],u[i][j+1],w,dt,dy,gamma);
        ubar_Lplus[i][j]=L;
        ubar_Rplus[i][j]=R;
        }
    }
    return {ubar_Lplus,ubar_Rplus};
}

StateVector calcFluxSLIC(const StateVector& ubar_R_plus0,const StateVector& ubar_L_plus1,const double& d_space,const double& dt,double gamma, double w){
    //Calculates the flux for a single point flux[i,j] refers to f_(i+1/2,j+1/2) 
    StateVector flux;
    flux = getFluxFORCE(ubar_R_plus0,ubar_L_plus1,d_space,dt,gamma);
         
    return flux;   
}
std::vector<std::vector<StateVector > > get_FluxSLIC_xarray(const std::vector<std::vector<StateVector > >& ubar_L_plusx,const std::vector<std::vector<StateVector > >& ubar_R_plusx, double w,double dt, double dx,double gamma,double num_xcells,double num_ycells){
    //Fills an array with the x fluxes at every cell where ubar[i] and ubar[i+1] exist
    // f[i]=f(x_(i+1/2))
    //ubar is only defined on the first ring of ghost cells, therefore f can be defined along the [1] edge and not the [N-1] edge
    std::vector<std::vector<StateVector > > x_flux(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    for(size_t i=1;i<num_xcells+2;i++){
        for(size_t j=1;j<num_ycells+2;j++){
        x_flux[i][j] = calcFluxSLIC(ubar_R_plusx[i][j],ubar_L_plusx[i+1][j],dx,dt,gamma,w);
        }
    }
    return x_flux;
}
std::vector<std::vector<StateVector > > get_FluxSLIC_yarray(const std::vector<std::vector<StateVector > >& ubar_L_plusy,const std::vector<std::vector<StateVector > >& ubar_R_plusy, double w,double dt, double dy,double gamma,double num_xcells,double num_ycells){
    //Fills an array with the x fluxes at every cell where ubar[i] and ubar[i+1] exist
    // f[i]=f(x_(i+1/2))
    //ubar is only defined on the first ring of ghost cells, therefore f can be defined along the [1] edge and not the [N-1] edge
    std::vector<std::vector<StateVector > > y_flux(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    for(size_t i=1;i<num_xcells+2;i++){
        for(size_t j=1;j<num_ycells+2;j++){
        y_flux[i][j] = calcFluxSLIC(ubar_R_plusy[i][j],ubar_L_plusy[i][j+1],dy,dt,gamma,w);
        }
    }
    return y_flux;
}




int main(void){
    double C=0.8;
    double tEnd=0.25;
    double tStart=0;
    double x0=0;
    double xf=1;
    double y0=0;
    double yf=1;
    double num_xcells=100;
    double num_ycells=100;
    double gamma=1.4;
    double w=0;

    double dx=(xf-x0)/num_xcells;
    double dy=(yf-y0)/num_ycells;

    //Initialise x vector which marks the x point of  cell centres (including boundary padding cells)
    std::vector<double> x(num_xcells+4);
    for(int i=0;i<x.size();i++){
        x[i]=x0+(i-0.5)*dx;

    };
    
    std::vector<double> y(num_ycells+4);
    for(int i=0;i<y.size();i++){
        y[i]=y0+(i-0.5)*dy;

    };


    //Initialise flux vector where the first value is the flux to the left of the first real cell i.e. flux[0] corresponds to x[0.5] where x[0] is the ghost cell
    std::vector<std::vector<StateVector > > x_flux(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    
    std::vector<std::vector<StateVector > > y_flux(num_xcells+4, std::vector<StateVector >(num_ycells+4));

    //Initialise u(x) vector
    std::vector<std::vector<StateVector > > u(num_xcells+4, std::vector<StateVector >(num_ycells+4));

    //Initialise vectors to store intermediate ubar_L_plus, ubar_R_plus values:
    std::vector<std::vector<StateVector > > ubar_L_plusx(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > ubar_R_plusx(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > ubar_L_plusy(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    std::vector<std::vector<StateVector > > ubar_R_plusy(num_xcells+4, std::vector<StateVector >(num_ycells+4));
    


    //Array to store updated values
    std::vector<std::vector<StateVector > > uPlus1(num_xcells+4, std::vector<StateVector >(num_ycells+4));

    //Set initial state

    u=set_u0(x,y,gamma);
 
    
    double t=tStart;
    int save_interval=10; //Saves per n iterations
    int count=0; //counts number of update steps performed

    do {
        double dt=computeTimeStep(u,dx,dy,C,gamma,num_xcells,num_ycells);
        t+=dt;
        count++;

        //Perform dimensionally split update, x then y
        
        //First update the u_bar arrays to correspond to x flux
        auto [ubar_L_plusx,ubar_R_plusx] = get_ubar_plus_xarrays(u,w,dt,dx,gamma,num_xcells,num_ycells);
  

        //Then calculate the xflux array

        x_flux=get_FluxSLIC_xarray(ubar_L_plusx,ubar_R_plusx,w,dt,dx,gamma,num_xcells,num_ycells);

        //Loop through x and y positions and update all real cells using x flux

        for(size_t i=2; i<num_xcells+2; i++){
            for(size_t j=2; j<num_ycells+2; j++){
            for(size_t k=0;k<4;k++){
                uPlus1[i][j][k]=u[i][j][k] + (dt/dx) * (x_flux[i-1][j][k]-x_flux[i][j][k]);
            }
        }
    }
  
        u=uPlus1;
        update_bcs_trans(u);

        //Now repeat process for y flux updates

        auto [ubar_L_plusy,ubar_R_plusy] = get_ubar_plus_yarrays(u,w,dt,dy,gamma,num_xcells,num_ycells);


        y_flux=get_FluxSLIC_yarray(ubar_L_plusy,ubar_R_plusy,w,dt,dy,gamma,num_xcells,num_ycells);


        for(size_t i=2; i<num_xcells+2; i++){
            for(size_t j=2; j<num_ycells+2; j++){
            for(size_t k=0;k<4;k++){
                uPlus1[i][j][k]=u[i][j][k] + (dt/dy) * (y_flux[i][j-1][k]-y_flux[i][j][k]);
            }
        }
    }

        u=uPlus1;
        update_bcs_trans(u);

        if(count%save_interval==0){
            std::string name= "2D_" + std::to_string(t);
            save_to_file(u,x,num_xcells,y,num_ycells,gamma,name);
        }
        std::cout<<"Completed Iteration"<<count<<std::endl;

      
        

    }while (t<tEnd);

}

    
    


