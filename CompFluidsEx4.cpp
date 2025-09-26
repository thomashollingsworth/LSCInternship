#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <filesystem>  // C++17



std::vector<double> primitiveToConservative(std::vector<double> prim,double gamma){
    //Convert a 1D vector of density,v,p to the conserved quantities density,momentum,energy
    std::vector<double> conserved(3);
    conserved[0]= prim[0];
    conserved[1]=prim[0]*prim[1];
    conserved[2] = prim[2]/((gamma-1)) +0.5*prim[0]*prim[1]*prim[1];

    return conserved;

}

std::vector<double> conservativeToPrimtive(std::vector<double> conserved,double gamma){
    //Convert a 1D vector of density,mom., energy to density,velocity,pressure
    std::vector<double> prim(3);
    prim[0]= conserved[0];
    prim[1]=conserved[1]/conserved[0];
    prim[2] = (gamma-1)*(conserved[2]-0.5*(conserved[1]*conserved[1])/conserved[0]);

    return prim;

}
    
void save_to_file(std::vector<std::vector<double> >& u,std::vector<double> x, double num_xcells, double gamma, std::string filename,std::string dir){
    //Convert to primitive variables
    std::vector<std::vector<double> > output(num_xcells+2, std::vector<double>(3, 0.0));
    for(size_t i=0;i<output.size();i++){
        output[i]=conservativeToPrimtive(u[i],gamma);
    }

    try {
        std::filesystem::create_directories(dir); 
        
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    filename=dir+filename;

    std::ofstream outfile(filename);

    for(size_t i=1; i<u.size()-2;i++){
            outfile << x[i] << " " << output[i][0] << " "<< output[i][1]<< " "<< output[i][2]<<std::endl;

}
    std::string gnuplot_cmd = "gnuplot -e \"filename='" + filename + "'\" 1Dplot_script.plt";
    std::system(gnuplot_cmd.c_str());
    std::cout << "Plot saved to " <<filename<< ".png" << std::endl;
    //std::filesystem::remove(filename);
}

void update_bcs_periodic(std::vector<std::vector<double> >& u){
    u[0]=u[u.size()-2];
    u[u.size()-1]= u[1];

    
    
}

void update_bcs_trans(std::vector<std::vector<double> >& u){
    u[0]=u[1];
    u[u.size()-1]= u[u.size()-2];

}

double computeTimeStep(const std::vector<std::vector<double> >& u,double dx,double C,double gamma){
    std::vector<double> prim0(3);
    prim0=conservativeToPrimtive(u[0],gamma);
    double max=std::abs(prim0[1])+std::sqrt(gamma*prim0[2]/prim0[0]);
    
    std::vector<double> prim(3);
    
    for(size_t i=1; i<u.size();i++){
        prim=conservativeToPrimtive(u[i],gamma);
        double new_max=std::abs(prim[1])+std::sqrt(gamma*prim[2]/prim[0]);
        if(new_max> max){
            max=new_max;
        }
        
    }
 return C*dx/max;
}

//Set the initial t=0 condition from a given xarray

std::vector<std::vector<double> > set_u0(std::vector<double> x,double gamma){
    std::vector<double> prim(3);
    std::vector<std::vector<double> > u0(x.size(),std::vector<double>(3));

    
    for(size_t i=0; i<x.size();i++){

    //Set initial condition in primitive coords
        if(x[i]<=0.5){
            prim[0]=1.;
            prim[1]=0.;
            prim[2]=0.01;
        }
        else{
            prim[0]=1.0;
            prim[1]=0.;
            prim[2]=100.;
        }
       

        u0[i]=primitiveToConservative(prim,gamma);

    }
    return u0;
    }
    
std::vector<double> euler_flux(std::vector<double> u0,double gamma){
    std::vector<double> flux(3);
    std::vector<double> prim=conservativeToPrimtive(u0,gamma);
    flux[0]=prim[0]*prim[1];
    flux[1]=prim[0]*prim[1]*prim[1] +prim[2];
    flux[2]= (u0[2]+prim[2])*prim[1];
    return flux;

}

double burger_flux(const double& u_0){
    return 0.5*u_0*u_0;
}

std::vector<double> getFluxRI(const std::vector<double>& u_0,const std::vector<double>& u_1,const double& dx,const double& dt,double gamma){
    //returns the flux f_(i+1/2)
    std::vector<double> f_plus=euler_flux(u_1,gamma);
    std::vector<double> f_minus=euler_flux(u_0,gamma);
    std::vector<double> u_plus_half(3);
    for(size_t i=0;i<3;i++){
        u_plus_half[i]=0.5*(u_0[i]+u_1[i]) -0.5*dt/dx*(f_plus[i]-f_minus[i]);

    }
  
    return euler_flux(u_plus_half,gamma);


}

std::vector<double> getFluxLF(const std::vector<double>& u_0,const std::vector<double>& u_1,const double& dx,const double& dt,double gamma){

    //returns the flux f_(i+1/2)
    std::vector<double> f_plus=euler_flux(u_1,gamma);
    std::vector<double> f_minus=euler_flux(u_0,gamma);
    std::vector<double> flux(3);
    for(size_t i=0;i<3;i++){
        flux[i]=0.5*dx/dt*(u_0[i]-u_1[i])+0.5*(f_minus[i]+f_plus[i]);
    }

    return flux;

}

std::vector<double> getFluxFORCE(const std::vector<double>& u_0, const std::vector<double>& u_1,const double& dx,const double& dt,double gamma){
    std::vector<double> flux(3);
    std::vector<double> fluxLF=getFluxLF(u_0,u_1,dx,dt,gamma);
    std::vector<double> fluxRI=getFluxRI(u_0,u_1,dx,dt,gamma);
    for(size_t i=0;i<3;i++){
        flux[i]=0.5*(fluxLF[i]+fluxRI[i]);
    }
    return flux;
}

double getBurgerFluxGOD(const double& u_0, const double& u_1,const double& dx,const double& dt,const double& a){
    double u_plus_half;
    bool left_larger=(u_0>u_1);
    
    if(left_larger){
        double s= 0.5*(u_0+u_1);
        if(s>0){
            u_plus_half=u_0;
        } else{
            u_plus_half=u_1;

        }
    }else{
        if(u_0>0){
            u_plus_half=u_0;
        }else if(u_1<0){
            u_plus_half=u_1;
        }else{
            u_plus_half=0;
        }
    }
    return burger_flux(u_plus_half);




}


double get_p0_PV(const std::vector<double>& u_0, const std::vector<double>& u_1,double gamma,const double TOL=1e-6){
    //Use PV approach to guess an initial p value for Newton-Raphson Solver
    std::vector<double> prim_L=conservativeToPrimtive(u_0,gamma);
    std::vector<double> prim_R=conservativeToPrimtive(u_1,gamma);
    
    double p_PV= 0.5*(prim_L[2]+prim_R[2])-0.125*(prim_R[1]-prim_L[1])*(prim_L[0]+prim_R[0])*(std::sqrt(gamma*prim_L[2]/prim_L[0])+std::sqrt(gamma*prim_R[2]/prim_R[0]));

    double output =std::max(p_PV,TOL);
    return output;
}

double  get_f_K(const std::vector<double>& u_K,double gamma,double p_star){
    //Can replace K with (left or right)/(u_0,u_1) 
    std::vector<double> prim_K=conservativeToPrimtive(u_K,gamma);
    double f_K;
    if(p_star>prim_K[2]){
        double A_K=2/((gamma+1)*prim_K[0]);
        double B_K= (gamma-1)/(gamma+1)*prim_K[2];
        f_K= (p_star-prim_K[2])*std::sqrt(A_K/(p_star+B_K));
    }else{
        double power= (gamma-1)/(2*gamma);
        double c_K=std::sqrt(gamma*prim_K[2]/prim_K[0]);
        f_K=(2*c_K)/(gamma-1)*(std::pow(p_star/prim_K[2],power)-1);
    }
    return f_K;
}

double  get_f_derivative_K(const std::vector<double>& u_K,double gamma,double p_star){
    //Can replace K with (left or right)/(u_0,u_1) 
    std::vector<double> prim_K=conservativeToPrimtive(u_K,gamma);
    double f_derivative_K;
    if(p_star>prim_K[2]){
        double A_K=2/((gamma+1)*prim_K[0]);
        double B_K= (gamma-1)/(gamma+1)*prim_K[2];
        f_derivative_K= std::sqrt(A_K/(p_star+B_K))*(1-0.5*(p_star-prim_K[2])/(p_star+B_K));
    }else{
        double power= (gamma+1)/(2*gamma);
        double c_K=std::sqrt(gamma*prim_K[2]/prim_K[0]);
        f_derivative_K=1.0/(c_K*prim_K[0])*std::pow(prim_K[2]/p_star,power);
    }
    return f_derivative_K;
}

double get_p_star(const std::vector<double>& u_0, const std::vector<double>& u_1,double gamma,const double TOL=1e-6){
    //   0/1 and L/R are interchangeable
    //Use PV method to estimate p_0
    double p_0=get_p0_PV(u_0,u_1,gamma,TOL);
    double p_old=p_0;
    std::vector<double> prim_0=conservativeToPrimtive(u_0,gamma);
    std::vector<double> prim_1=conservativeToPrimtive(u_1,gamma);
    double delta_v=prim_1[1]-prim_0[1];
    double change;
    
    do{
        double p_new=p_old-(get_f_K(u_0,gamma,p_old)+get_f_K(u_1,gamma,p_old)+delta_v)/(get_f_derivative_K(u_1,gamma,p_old)+get_f_derivative_K(u_0,gamma,p_old));
        change=2*std::abs(p_new-p_old)/(p_new+p_old);
        p_old=p_new;
    }while(change>TOL);
return p_old;
}

double get_v_star(const std::vector<double>& u_0, const std::vector<double>& u_1,double gamma,const double p_star){
    std::vector<double> prim_L=conservativeToPrimtive(u_0,gamma);
    std::vector<double> prim_R=conservativeToPrimtive(u_1,gamma);
    return 0.5*(prim_L[1]+prim_R[1])+0.5*(get_f_K(u_1,gamma,p_star)-get_f_K(u_0,gamma,p_star));
}

std::vector<double> get_exact_u_half(const std::vector<double>& u_L, const std::vector<double>& u_R,double v_star,double p_star,double gamma,double S=0){
    //Uses p_star and v_star and the large flow chart to analytically evaluate W(rho,u,p)
    //Then convert to u(rho,momentum,Energy) for the x=0 i.e. cell boundary point
    //S=x/t therefore is 0

    std::vector<double> prim_L=conservativeToPrimtive(u_L,gamma);
    std::vector<double> prim_R=conservativeToPrimtive(u_R,gamma);
    double c_L=std::sqrt(gamma*prim_L[2]/prim_L[0]);
    double c_R=std::sqrt(gamma*prim_R[2]/prim_R[0]);

    std::vector<double> W(3);
 
    //Defining booleans for the flow chart:
    bool left_of_contact_discontinuity = v_star>0;
    bool left_wave_shock =p_star>prim_L[2];
    bool right_wave_shock=p_star>prim_R[2];

    if(left_of_contact_discontinuity){
        if(left_wave_shock){
            double S_L = prim_L[1]-c_L*std::sqrt((gamma+1)/(2*gamma)*p_star/prim_L[2]+(gamma-1)/(2*gamma));
            if(S_L>0){
                for(size_t i=0;i<4;i++){
                    W[i]=prim_L[i];
                }} else {
                    W[0]=prim_L[0]*std::pow(p_star/prim_L[2],1/gamma);
                    W[1]=v_star;
                    W[2]=p_star;

                }}else{
                    //left wave is rarefaction
                    double S_HL = prim_L[1]-c_L;
                    if(S_HL>0){
                        for(size_t i=0;i<4;i++){
                        W[i]=prim_L[i];
                }}else{
                    double S_TL=v_star-c_L*std::pow(p_star/prim_L[2],(gamma-1)/(2*gamma));
                    if(S_TL<0){
                    W[0]=prim_L[0]*std::pow(p_star/prim_L[2],1/gamma);
                    W[1]=v_star;
                    W[2]=p_star;

                    }else{//W_LFAN
                        W[0]=prim_L[0]*std::pow(2/(gamma+1)+(gamma-1)/(c_L*(gamma+1))*prim_L[1],2/(gamma-1));
                        W[1]=2/(gamma+1)*(c_L+(gamma-1)/2*prim_L[1]);
                        W[2]=prim_L[2]*std::pow(2/(gamma+1)+(gamma-1)/(c_L*(gamma+1))*prim_L[1],(2*gamma)/(gamma-1));
                    }
                }

                    }
                }else{
                    if(right_wave_shock){
                        double S_R=prim_R[1]+c_R*std::sqrt((gamma+1)/(2*gamma)*p_star/prim_R[2]+(gamma-1)/(2*gamma));
                        if(S_R>0){
                        W[0]=prim_R[0]*std::pow(p_star/prim_R[2],1/gamma);
                        W[1]=v_star;
                        W[2]=p_star;

                        }else{
                            for(size_t i=0;i<4;i++){
                            W[i]=prim_R[i];
                }

                        }
                    }else{//Right wave is a rarefaction
                        double S_TR=v_star+c_R*std::pow(p_star/prim_R[2],(gamma-1)/(2*gamma));
                        if(S_TR>0){
                        W[0]=prim_R[0]*std::pow(p_star/prim_R[2],1/gamma);
                        W[1]=v_star;
                        W[2]=p_star;
                        }else{
                        double S_HR= prim_R[1]+c_R;
                        if(S_HR>0){//W_RFAN
                        W[0]=prim_R[0]*std::pow(2/(gamma+1)-(gamma-1)/(c_R*(gamma+1))*prim_R[1],2/(gamma-1));
                        W[1]=2/(gamma+1)*(-c_R+(gamma-1)/2*prim_R[1]);
                        W[2]=prim_R[2]*std::pow(2/(gamma+1)-(gamma-1)/(c_R*(gamma+1))*prim_R[1],(2*gamma)/(gamma-1));
                        }else{
                            for(size_t i=0;i<4;i++){
                            W[i]=prim_R[i];
                }   } } }}

                std::vector<double> u_out=primitiveToConservative(W,gamma);
                return u_out;
            }


    


int main(void){
    double C=0.8;
    double tEnd=0.04;
    double tStart=0;
    double x0=0;
    double xf=1;
    double num_xcells=250;
    double gamma=1.4;

    double dx=(xf-x0)/num_xcells;

    //Initialise x vector which marks the x point of  cell centres (including boundary padding cells)
    std::vector<double> x(num_xcells+2);
    for(int i=0;i<x.size();i++){
        x[i]=x0+(i-0.5)*dx;

    };

    //Initialise flux vector where the first value is the flux to the left of the first real cell i.e. flux[0] corresponds to x[0.5] where x[0] is the ghost cell
    std::vector<std::vector<double> > flux(num_xcells+2,std::vector<double>(3, 0.0));

    //Initialise u(x) vector
    std::vector<std::vector<double> > u(num_xcells+2, std::vector<double>(3, 0.0));
    
    u=set_u0(x,gamma);
    std::vector<std::vector<double> > uPlus1(num_xcells+2, std::vector<double>(3, 0.0));
    
    double t=tStart;
    int save_interval=10;
    int count=0; //counts number of update steps performed

    do {
        update_bcs_trans(u);
        double dt=computeTimeStep(u,dx,C,gamma);
        t+=dt;
        count++;
        
        //Update flux array

        for(size_t i=0; i<num_xcells+1;i++){
            flux[i]=getFluxFORCE(u[i],u[i+1],dx,dt,gamma);
        }

        //Use flux arrays to update u

        for(size_t i=1; i<num_xcells+1; i++){
            for(size_t v=0;v<3;v++){
                uPlus1[i][v]=u[i][v] - (dt/dx) * (flux[i][v]-flux[i-1][v]);
            }
        }
       
        u=uPlus1;
      
        if(count%save_interval==0){
            std::string name= "EulerApprox_" + std::to_string(t);
             //Output the data
            std::string dir = "output_approx/";
            save_to_file(u,x,num_xcells,gamma,name,dir);
        }

    }while (t<tEnd);

    //Repeat with exact Solver
    u=set_u0(x,gamma);
    t=tStart;
    save_interval=10;
    count=0; //counts number of update steps performed

    do {
        update_bcs_trans(u);
        double dt=computeTimeStep(u,dx,C,gamma);
        t+=dt;
        count++;

        

        //Update flux array

        for(size_t i=0; i<num_xcells+1;i++){
            double p_star=get_p_star(u[i],u[i+1],gamma);
            double v_star=get_v_star(u[i],u[i+1],gamma,p_star);
            
            flux[i]=euler_flux(get_exact_u_half(u[i],u[i+1],v_star,p_star,gamma),gamma);
        }

        //Use flux arrays to update u

        for(size_t i=1; i<num_xcells+1; i++){
            for(size_t v=0;v<3;v++){
                uPlus1[i][v]=u[i][v] - (dt/dx) * (flux[i][v]-flux[i-1][v]);
            }
        }
       
        u=uPlus1;
      
        if(count%save_interval==0){
            std::string name= "EulerExact_" + std::to_string(t);
            std::string dir = "output_exact/";
            save_to_file(u,x,num_xcells,gamma,name,dir);
        }

    }while (t<tEnd);

}

    



