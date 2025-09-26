#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <filesystem>  // C++17



std::vector<double> primitiveToConservative(std::vector<double> prim,double gamma){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    std::vector<double> conserved(8);
    conserved[0]= prim[0];
    conserved[1]=prim[0]*prim[1];
    conserved[2]=prim[0]*prim[2];
    conserved[3]=prim[0]*prim[3];
    conserved[4]=prim[4]/((gamma-1)) +0.5*prim[0]*(prim[1]*prim[1]+prim[2]*prim[2]+prim[3]*prim[3])+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]);
    conserved[5]=prim[5];
    conserved[6]=prim[6];
    conserved[7]=prim[7];

    return conserved;

}

std::vector<double> conservativeToPrimitive(std::vector<double> conserved,double gamma){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    std::vector<double> prim(8);
    prim[0]= conserved[0];
    prim[1]=conserved[1]/conserved[0];
    prim[2]=conserved[2]/conserved[0];
    prim[3]=conserved[3]/conserved[0];
    prim[4] = (gamma-1)*(conserved[4]-0.5*(conserved[1]*conserved[1]+conserved[2]*conserved[2]+conserved[3]*conserved[3])/conserved[0]-0.5*(conserved[5]*conserved[5]+conserved[6]*conserved[6]+conserved[7]*conserved[7]));
    prim[5]=conserved[5];
    prim[6]=conserved[6];
    prim[7]=conserved[7];

    return prim;

}
    
void save_to_file(std::vector<std::vector<double> >& u,std::vector<double> x, double num_xcells, double gamma, std::string filename,std::string dir){
    //Convert to primitive variables
    std::vector<std::vector<double> > output(num_xcells+4, std::vector<double>(8, 0.0));
    for(size_t i=0;i<num_xcells+4;i++){
        output[i]=conservativeToPrimitive(u[i],gamma);
    }

    try {
        std::filesystem::create_directories(dir); 
        
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    filename=dir+filename;

    std::ofstream outfile(filename);

    for(size_t i=2; i<num_xcells+2;i++){
            outfile << x[i] << " " << output[i][0] << " "<< output[i][1]<< " "<< output[i][2]<< " "<< output[i][3]<< " "<< output[i][4]<< " "<< output[i][5]<< " "<< output[i][6]<< " "<< output[i][7]<<std::endl;

}
    std::string gnuplot_cmd = "gnuplot -e \"filename='" + filename + "'\" MHDplot_script.plt";
    std::system(gnuplot_cmd.c_str());
    std::cout << "Plot saved to " <<filename<< ".png" << std::endl;
    //std::filesystem::remove(filename);
}

void update_bcs_trans(std::vector<std::vector<double> >& u,double num_xcells){
    u[0]=u[2];
    u[1]=u[2];
    u[num_xcells+3]= u[num_xcells+1];
    u[num_xcells+2]= u[num_xcells+1];

}

double computeTimeStep(const std::vector<std::vector<double> >& u,double dx,double C,double gamma,double num_xcells){
    //Prim: rho,v_x,v_y,v_z,p,B_y,B_x,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    std::vector<double> prim0(8);
    prim0=conservativeToPrimitive(u[0],gamma);
 
    //Breaking calc into steps for readibility
    double factor0= (gamma*prim0[4]+prim0[5]*prim0[5]+prim0[6]*prim0[6]+prim0[7]*prim0[7])/prim0[0];
    double c_f0=std::sqrt(0.5*(factor0+std::sqrt(factor0*factor0-4.0*gamma*prim0[4]*prim0[5]*prim0[5]/(prim0[0]*prim0[0]))));
    double max=std::sqrt(prim0[1]*prim0[1]+prim0[2]*prim0[2]+prim0[3]*prim0[3])+c_f0;
    

    
    std::vector<double> prim(8);
    
    for(size_t i=2; i<num_xcells+2;i++){
        prim=conservativeToPrimitive(u[i],gamma);
        double factor= (gamma*prim[4]+prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])/prim[0];
        double c_f=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*prim[5]*prim[5]/(prim[0]*prim[0]))));
        double new_max=std::sqrt(prim[1]*prim[1]+prim[2]*prim[2]+prim[3]*prim[3])+c_f;

        if(new_max> max){
            max=new_max;
        }
        
    }
 return C*dx/max;
}

//Set the initial t=0 condition from a given xarray

std::vector<std::vector<double> > set_u0(std::vector<double> x,double gamma){
    std::vector<double> prim(8);
    std::vector<std::vector<double> > u0(x.size(),std::vector<double>(8));
    double pi=std::atan(1)*4;
    for(size_t i=0; i<x.size();i++){

        //1.08,1.2,0.01,0.5,3.6/√4π,2/√4π,2/√4π,0.95
        

    //Set initial condition in primitive coords
        if(x[i]<=400){
            prim[0]=1.;
            prim[1]=0.;
            prim[2]=0.;
            prim[3]=0.;
            prim[4]=1.;
            prim[5]=0.75;
            prim[6]=1.;
            prim[7]=0.;
        }
        //1,0,0,0,4/√4π,2/√4π,2/√4π,1
        else{
            prim[0]=0.125;
            prim[1]=0.;
            prim[2]=0.;
            prim[3]=0.;
            prim[4]=0.1;
            prim[5]=0.75;
            prim[6]=-1.0;
            prim[7]=0.0;
        }
        u0[i]=primitiveToConservative(prim,gamma);
    }
   
    return u0;
    }
    
std::vector<double> MHD_flux(const std::vector<double>& u0,double gamma){
    //Prim: rho,v_x,v_y,v_z,p,B_x,B_y,B_z (B_x is const.)
    //Conserved: rho, mom_x, mom_y, mom_z, U, B_x, B_y, B_z
    std::vector<double> flux(8);
    std::vector<double> prim=conservativeToPrimitive(u0,gamma);
    flux[0]=prim[0]*prim[1];
    flux[1]=prim[0]*prim[1]*prim[1]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[5]*prim[5];
    flux[2]= prim[0]*prim[1]*prim[2]-prim[5]*prim[6];
    flux[3]= prim[0]*prim[1]*prim[3]-prim[5]*prim[7];
    flux[4]=(u0[4]+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[1]-(prim[1]*prim[5]+prim[2]*prim[6]+prim[3]*prim[7])*prim[5];
    flux[5]=0.;//B_x is constant for 1D problem
    flux[6]=prim[6]*prim[1]-prim[5]*prim[2];
    flux[7]=prim[7]*prim[1]-prim[5]*prim[3];
    
    return flux;

}

std::vector<double> getFluxRI(const std::vector<double>& u_0,const std::vector<double>& u_1,const double& dx,const double& dt,double gamma){
    //returns the flux f_(i+1/2)
    std::vector<double> f_plus=MHD_flux(u_1,gamma);
    std::vector<double> f_minus=MHD_flux(u_0,gamma);
    std::vector<double> u_plus_half(8);
    for(size_t i=0;i<8;i++){
        u_plus_half[i]=0.5*(u_0[i]+u_1[i]) -0.5*dt/dx*(f_plus[i]-f_minus[i]);

    }
  
    return MHD_flux(u_plus_half,gamma);


}

std::vector<double> getFluxLF(const std::vector<double>& u_0,const std::vector<double>& u_1,const double& dx,const double& dt,double gamma){

    //returns the flux f_(i+1/2)
    std::vector<double> f_plus=MHD_flux(u_1,gamma);
    std::vector<double> f_minus=MHD_flux(u_0,gamma);
    std::vector<double> flux(8);
    for(size_t i=0;i<8;i++){
        flux[i]=0.5*dx/dt*(u_0[i]-u_1[i])+0.5*(f_minus[i]+f_plus[i]);
    }

    return flux;

}

std::vector<double> getFluxFORCE(const std::vector<double>& u_0, const std::vector<double>& u_1,const double& dx,const double& dt,double gamma){
    std::vector<double> flux(8);
    std::vector<double> fluxLF=getFluxLF(u_0,u_1,dx,dt,gamma);
    std::vector<double> fluxRI=getFluxRI(u_0,u_1,dx,dt,gamma);
    for(size_t i=0;i<8;i++){
        flux[i]=0.5*(fluxLF[i]+fluxRI[i]);
    }
    return flux;
}

std::vector<double> getDelta_plus(const std::vector<double>& u_0, const std::vector<double>& u_plus){

    std::vector<double> delta_plus(8);
    for(size_t i=0;i<8;i++){
        delta_plus[i]=u_plus[i]-u_0[i];
    }
    return delta_plus;
}
std::vector<double> getDelta_minus(const std::vector<double>& u_minus, const std::vector<double>& u_0){
    std::vector<double> delta_minus(8);
    for(size_t i=0;i<8;i++){
        delta_minus[i]=u_0[i]-u_minus[i];
    }
    return delta_minus;
}
std::vector<double> getDelta(const std::vector<double>& u_minus,const std::vector<double>& u_0, const std::vector<double>& u_plus,double w){
std::vector<double> delta_minus = getDelta_minus(u_minus,u_0);
std::vector<double> delta_plus = getDelta_plus(u_0,u_plus);
std::vector<double> delta(8);
    for(size_t i=0;i<8;i++){
        delta[i]=0.5*(1+w)*delta_minus[i] + 0.5*(1-w)*delta_plus[i];
    }
    return delta;
}

std::vector<double> get_r(const std::vector<double>& u_minus,const std::vector<double>& u_0, const std::vector<double>& u_plus){
    std::vector<double> delta_minus = getDelta_minus(u_minus,u_0);
    std::vector<double> delta_plus = getDelta_plus(u_0,u_plus);
    std::vector<double> r(8);
        for(size_t i=0;i<8;i++){
    
        r[i]=delta_minus[i]/(delta_plus[i]+1e-8);
        }
        return r;
}
std::vector<double> getXi_Minbee(std::vector<double> r){
    std::vector<double> xi(8);
    for(size_t i=0;i<8;i++){
        if(r[i]<=0){
            xi[i]=0;
        } else if(r[i]<=1){
            xi[i]=r[i];

        } else{
            double xi_R = 2/(1+r[i]);
            xi[i]=std::fmin(1.0,xi_R);
        }
    }
    return xi;

}
std::tuple< std::vector<double>,std::vector<double>  > calc_ubar(const std::vector<double>& u_minus,const std::vector<double>& u_0, const std::vector<double>& u_plus,double w){
    //returns {ubar_L,ubar_R}
    std::vector<double> ubar_L(8);
    std::vector<double> ubar_R(8);
    std::vector<double> delta=getDelta(u_minus,u_0,u_plus,w);
    std::vector<double> r= get_r(u_minus,u_0,u_plus);
    std::vector<double> xi= getXi_Minbee(r);
    for(size_t j=0;j<8;j++){
        //Looping over variables
            ubar_L[j]=u_0[j]-0.5*xi[j]*delta[j];
            ubar_R[j]=u_0[j]+0.5*xi[j]*delta[j];
        }  
        return {ubar_L,ubar_R};
    }

std::tuple< std::vector<double>,std::vector<double> > calc_ubar_plus(const std::vector<double>& u_minus,const std::vector<double>& u_0, const std::vector<double>& u_plus,double w,double dt, double dx,double gamma){
    //returns {ubar_L^(n+1/2),ubar_R^(n+1/2) }
    auto [ubar_L, ubar_R]= calc_ubar(u_minus,u_0,u_plus,w);
    std::vector<double> ubar_L_plus(8);
    std::vector<double> ubar_R_plus(8);
    std::vector<double> f_R= MHD_flux(ubar_R,gamma);
    std::vector<double> f_L= MHD_flux(ubar_L,gamma);
    for(size_t j=0;j<8;j++){//Looping over variables
            ubar_L_plus[j]= ubar_L[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            ubar_R_plus[j]= ubar_R[j]-0.5*dt/dx*(f_R[j]-f_L[j]);
            
        }  
        return {ubar_L_plus,ubar_R_plus};

    }
std::tuple< std::vector<std::vector<double> >,std::vector<std::vector<double > > > get_ubar_plus_xarrays(const std::vector<std::vector<double> >& u, double w,double dt, double dx,double gamma,double num_xcells){
    std::vector<std::vector<double> > ubar_Rplus(num_xcells+4, std::vector<double>(8,0.0));
    std::vector<std::vector<double> > ubar_Lplus(num_xcells+4, std::vector<double>(8,0.0));
    
    for(size_t i=1;i<num_xcells+3;i++){
        auto [L, R] = calc_ubar_plus(u[i-1],u[i],u[i+1],w,dt,dx,gamma);
        ubar_Lplus[i]=L;
        ubar_Rplus[i]=R;  
    }
    update_bcs_trans(ubar_Lplus,num_xcells);
    update_bcs_trans(ubar_Rplus,num_xcells);
    
    return {ubar_Lplus,ubar_Rplus};
}

std::vector<double> calcFluxSLIC(const std::vector<double>& ubar_R_plus0,const std::vector<double>& ubar_L_plus1,const double& dx,const double& dt,double gamma, double w){
    //Calculates the flux for a single point flux[i,j] refers to f_(i+1/2,j+1/2) 
    std::vector<double > flux(8);
    flux = getFluxFORCE(ubar_R_plus0,ubar_L_plus1,dx,dt,gamma);
         
    return flux;   
}
std::vector<std::vector<double> > get_FluxSLIC_xarray(const std::vector<std::vector<double> >& ubar_L_plusx,const std::vector<std::vector<double> >& ubar_R_plusx, double w,double dt, double dx,double gamma,double num_xcells){
    //Fills an array with the x fluxes at every cell where ubar[i] and ubar[i+1] exist
    // f[i]=f(x_(i+1/2))
    //ubar is only defined on the first ring of ghost cells, therefore f can be defined along the [1] edge and not the [N-1] edge
    std::vector<std::vector<double> > x_flux(num_xcells+4,std::vector<double>(8,0.0));
    for(size_t i=1;i<num_xcells+2;i++){
      
        x_flux[i] = calcFluxSLIC(ubar_R_plusx[i],ubar_L_plusx[i+1],dx,dt,gamma,w);
        
    }
    return x_flux;
}


//HLL METHOD
double get_cf(double gamma, const std::vector<double>& u){
    
    // u is the vector of parameters at lattice site
    std::vector<double> prim =conservativeToPrimitive(u,gamma);
    
    double factor= (gamma*prim[4]+u[5]*u[5]+u[6]*u[6]+u[7]*u[7])/u[0];
    double c_f=std::sqrt(0.5*(factor+std::sqrt(factor*factor-4.0*gamma*prim[4]*u[5]*u[5]/(u[0]*u[0]))));
    return c_f;

}


std::tuple< double, double > get_SL_R(const std::vector<double>& u_L,const std::vector<double>& u_R, double gamma){
    std::vector<double> prim_L=conservativeToPrimitive(u_L,gamma);
    std::vector<double> prim_R=conservativeToPrimitive(u_R,gamma);
    double v_xL = prim_L[1];
    double v_xR = prim_R[1];
    double c_fL=get_cf(gamma,u_L);
    double c_fR=get_cf(gamma,u_R);
    double S_L = std::min(v_xL,v_xR)-std::max(c_fL,c_fR);
    double S_R = std::max(v_xL,v_xR)+std::max(c_fL,c_fR);

    return {S_L,S_R};

}

std::vector<double> get_u_HLL(const std::vector<double>& u_L,const std::vector<double>& u_R, double gamma){
    auto [S_L,S_R]= get_SL_R(u_L,u_R,gamma);

    std::vector<double> f_L=MHD_flux(u_L,gamma);
    std::vector<double> f_R=MHD_flux(u_R,gamma);
    std::vector<double> u_HLL(8);

    for(size_t i=0;i<8;i++){
        u_HLL[i]=1/(S_R-S_L)*(S_R*u_R[i]-S_L*u_L[i]+f_L[i]-f_R[i]);

    }
    return u_HLL;
}



std::vector<double> get_u_star(const std::vector<double>& u_L,const std::vector<double>& u_R, double gamma,bool L=false){
    //L bool dictates whether you are calculating u_starL or u_starR 
    
    std::vector<double> u_HLL=get_u_HLL(u_L,u_R,gamma);
    std::vector<double> prim_HLL = conservativeToPrimitive(u_HLL,gamma);
    auto [S_L,S_R]= get_SL_R(u_L,u_R,gamma);
    std::vector<double> prim_L= conservativeToPrimitive(u_L,gamma);
    std::vector<double> prim_R= conservativeToPrimitive(u_R,gamma);
    
    double S; // Will be used as S_L for L case and S_R for R_case
    std::vector<double> prim(8); // Will be used as L for L case and R for R_case
    std::vector<double> u(8);


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
    
    std::vector<double> u_star(8);
    //Set B_star values to HLL values
    u_star[5]=u_HLL[5];
    u_star[6]=u_HLL[6];
    u_star[7]=u_HLL[7];

    //Calc. the E_star term
    double p_star= prim[0]*(S-prim[1])*(q_star-prim[1])+prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7])-prim[5]*prim[5]+u_star[5]*u_star[5];
    
    u_star[4]=u[4]*(S-prim[1])/(S-q_star)+((p_star*q_star-(prim[4]+0.5*(prim[5]*prim[5]+prim[6]*prim[6]+prim[7]*prim[7]))*prim[1])-(u_star[5]*(prim_HLL[5]*prim_HLL[1]+prim_HLL[6]*prim_HLL[2]+prim_HLL[7]*prim_HLL[3])-prim[5]*(prim[5]*prim[1]+prim[6]*prim[2]+prim[7]*prim[3])))/(S-q_star);

    //Calc. the momentum and density terms
    u_star[0]=prim[0]*(S-prim[1])/(S-q_star);
    u_star[1]=u_star[0]*q_star;
    u_star[2]=u_star[0]*prim[2]-(u_star[5]*u_star[6]-prim[5]*prim[6])/(S-q_star);
    u_star[3]=u_star[0]*prim[3]-(u_star[5]*u_star[7]-prim[5]*prim[7])/(S-q_star);
    
    return u_star;

}

std::vector<double> get_HLLC_flux(const std::vector<double>& u_L,const std::vector<double>& u_R, double gamma){
    //gets flux between L and R position
    auto [S_L,S_R]= get_SL_R(u_L,u_R,gamma);
    std::vector<double> f_L=MHD_flux(u_L,gamma);
    std::vector<double> f_R=MHD_flux(u_R,gamma);
    std::vector<double> u_starL=get_u_star(u_L,u_R,gamma,true);
    std::vector<double> u_starR=get_u_star(u_L,u_R,gamma,false);

    double q_star= u_starL[1]/u_starL[0];
    std::vector<double> flux_out(8);

    //Deciding what region/flux is being used
    if(S_L>=0){
        flux_out=f_L;

    }else if(S_L<q_star){
        for(size_t i=0;i<8;i++){
        flux_out[i]=f_L[i]+S_L*(u_starL[i]-u_L[i]);
    }}else if(q_star<=S_R){
        for(size_t i=0;i<8;i++){
        flux_out[i]=f_R[i]+S_R*(u_starR[i]-u_R[i]);

    }}else {
        flux_out=f_R; 
    }

    return flux_out;

}





int main(void){
    double C=1;
    double tEnd=80;
    double tStart=0;
    double x0=0;
    double xf=800;
    double num_xcells=800;
    double gamma=2 ;
    double w=1;

    double dx=(xf-x0)/num_xcells;

    //Initialise x vector which marks the x point of  cell centres (including boundary padding cells)
    std::vector<double> x(num_xcells+4);
    for(int i=0;i<x.size();i++){
        x[i]=x0+(i-0.5)*dx;

    };

    //Initialise flux vector where the first value is the flux to the left of the first real cell i.e. flux[0] corresponds to x[0.5] where x[0] is the ghost cell
    std::vector<std::vector<double> > flux(num_xcells+4,std::vector<double>(8, 0.0));
    std::vector<std::vector<double> > ubar_L_plus(num_xcells+4,std::vector<double>(8, 0.0));
    std::vector<std::vector<double> > ubar_R_plus(num_xcells+4,std::vector<double>(8, 0.0));

    //Initialise u(x) vector
    std::vector<std::vector<double> > u(num_xcells+4, std::vector<double>(8, 0.0));
    
    u=set_u0(x,gamma);
    
    std::vector<std::vector<double> > uPlus1(num_xcells+4, std::vector<double>(8, 0.0));
    
    double t=tStart;
    int save_interval=25;
    int count=0; //counts number of update steps performed

    do {
        
       
        double dt=computeTimeStep(u,dx,C,gamma,num_xcells);
   
        t+=dt;
        count++;
        update_bcs_trans(u,num_xcells);


        
        // //First update the u_bar arrays to correspond to (x) flux
        // auto [ubar_L_plusx,ubar_R_plusx] = get_ubar_plus_xarrays(u,w,dt,dx,gamma,num_xcells);
    

        // //Update flux array
        // flux=get_FluxSLIC_xarray(ubar_L_plusx,ubar_R_plusx,w,dt,dx,gamma,num_xcells);

        
        //Use flux to update all real x cells
        for(size_t i=2; i<num_xcells+2; i++){
            for(size_t j=0;j<8;j++){
                uPlus1[i][j]=u[i][j] - (dt/dx) * (get_HLLC_flux(u[i],u[i+1],gamma)[j]-get_HLLC_flux(u[i-1],u[i],gamma)[j]);
            
        }
    }


        u=uPlus1;
        update_bcs_trans(u,num_xcells);
      
        if(count%save_interval==0){
            std::string name= "Plot_" + std::to_string(t);
             //Output the data
            std::string dir = "HLLC/";
            save_to_file(u,x,num_xcells,gamma,name,dir);
        }

   

    }while (t<tEnd);

}

    



