#include "SLIC.h"


//Can be used for euler or MHD or any N-param system
template<typename StateVector>
StateVector getDelta_plus(const StateVector& u_0, const StateVector& u_plus){
    const size_t size= u_0.size();
    StateVector delta_plus;
    for(size_t i=0;i<size;i++){
        delta_plus[i]=u_plus[i]-u_0[i];
    }
    return delta_plus;
}
template<typename StateVector>
StateVector getDelta_minus(const StateVector& u_minus, const StateVector& u_0){
    const size_t size= u_0.size();
    StateVector delta_minus;
    for(size_t i=0;i<size;i++){
        delta_minus[i]=u_0[i]-u_minus[i];
    }
    return delta_minus;
}
template<typename StateVector>
StateVector getDelta(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w){
    const size_t size= u_0.size();
    StateVector delta_minus = getDelta_minus(u_minus,u_0);
    StateVector delta_plus = getDelta_plus(u_0,u_plus);
    StateVector delta;
    for(size_t i=0;i<size;i++){
        delta[i]=0.5*(1+w)*delta_minus[i] + 0.5*(1-w)*delta_plus[i];
    }
    return delta;
}
template<typename StateVector>
StateVector get_r(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus){
    const size_t size= u_0.size();
    StateVector delta_minus = getDelta_minus(u_minus,u_0);
    StateVector delta_plus = getDelta_plus(u_0,u_plus);
    StateVector r;
        for(size_t i=0;i<size;i++){
        r[i]=delta_minus[i]/delta_plus[i];
        }
        return r;
}
template<typename StateVector>
StateVector getXi_Minbee(const StateVector& r){
    const size_t size= r.size();
    StateVector xi;
    for(size_t i=0;i<size;i++){
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

template<typename StateVector>
std::tuple< StateVector,StateVector  > calc_ubar(const StateVector& u_minus,const StateVector& u_0, const StateVector& u_plus,double w){
    //returns {ubar_L,ubar_R}
    const size_t size= u_0.size();
    StateVector ubar_L;
    StateVector ubar_R;
    StateVector delta=getDelta(u_minus,u_0,u_plus,w);
    StateVector r= get_r(u_minus,u_0,u_plus);
    StateVector xi= getXi_Minbee(r);
    for(size_t j=0;j<size;j++){
        //Looping over variables
            ubar_L[j]=u_0[j]-0.5*xi[j]*delta[j];
            ubar_R[j]=u_0[j]+0.5*xi[j]*delta[j];
        }  
        return {ubar_L,ubar_R};
    }


template<typename StateVector>
std::tuple< StateVector,StateVector > calc_ubar_plus(const StateVector& ubarL, const StateVector& ubarR,const StateVector& fluxL, const StateVector& fluxR, const double dt, const double dx){
    //returns {ubar_L^(n+1/2),ubar_R^(n+1/2) }
 
    const size_t size= ubarL.size();
    StateVector ubar_L_plus;
    StateVector ubar_R_plus;

    for(size_t j=0;j<size;j++){//Looping over variables
            ubar_L_plus[j]= ubar_L[j]-0.5*dt/dx*(fluxR[j]-fluxL[j]);
            ubar_R_plus[j]= ubar_R[j]-0.5*dt/dx*(fluxR[j]-fluxL[j]);
            
        }  
        return {ubar_L_plus,ubar_R_plus};

    }



