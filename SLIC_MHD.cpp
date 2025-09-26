#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>

typedef std::array< double ,8> StateVector;

//Primitive to Conservative

StateVector PrimitiveToConservative(const StateVector& u , double gamma){
    StateVector v;
    v[0] = u[0]; //density
    v[1] = u[0] * u[1]; //vx
    v[2] = u[0] * u[2]; //vy
    v[3] = u[0] * u[3]; //vz
    v[5] = u[5]; //Bx
    v[6] = u[6]; //By
    v[7] = u[7]; //Bz
    v[4] = u[4] / (gamma-1) + 0.5*u[0]*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]) + 0.5*(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]); //energy
    return v;
}

//Conservative to Primative

StateVector ConservativeToPrimitive(const StateVector& u , double gamma){
    StateVector v;
    v[0] = u[0];
    v[1] = u[1] / u[0];
    v[2] = u[2] / u[0];
    v[3] = u[3] / u[0];
    v[5] = u[5];
    v[6] = u[6];
    v[7] = u[7];
    v[4] = (u[4] - 0.5*v[0]*(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]) - 0.5*(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]))*(gamma -1);//pressure
    return v;
}

void save_to_file(std::vector<StateVector >& u,std::vector<double> x, double num_xcells, double gamma, std::string filename,std::string dir){
    //Convert to primitive variables
    std::vector<StateVector > output(num_xcells+4);
    for(size_t i=0;i<num_xcells+4;i++){

        output[i]=ConservativeToPrimitive(u[i],gamma);
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


StateVector FluxDef(const StateVector& u , double gamma ){//conservative variables go into this function
    StateVector v = ConservativeToPrimitive(u , gamma);
    StateVector flux;
    flux[0] = u[1]; //rho v_x
    flux[1] = u[0]*v[1]*v[1] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ) - u[5]*u[5]; 
    flux[2] = u[0]*v[1]*v[2] - u[5]*u[6];
    flux[3] = u[0]*v[1]*v[3] - u[5]*u[7];
    flux[4] = (u[4] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ))*v[1] - (v[1]*v[5]+v[2]*v[6]+v[3]*v[7])*v[5];
    flux[5] = 0;
    flux[6] = v[6]*v[1]-v[5]*v[2];
    flux[7] = v[7]*v[1]-v[5]*v[3];
    return flux;
}
double ComputeTimeStep(const std::vector<StateVector>& u , double C, double dx, double gamma) {
    double maxSpeed = 0.0;

    for (const auto& state : u) {
        double rho = state[0];
        double mom_x = state[1];
        double mom_y = state[2];
        double mom_z = state[3];
        double E = state[4];
        double Bx = state[5];
        double By = state[6];
        double Bz = state[7];

        double u_x = mom_x / rho;
        double u_y = mom_y / rho;
        double u_z = mom_z / rho;
        double intermediate = 0.5 * rho *( u_x * u_x + u_y*u_y + u_z*u_z) + 0.5*(Bx*Bx + By*By + Bz*Bz);
        double BmagSquared = Bx*Bx + By*By + Bz*Bz;
        double pressure = (gamma - 1.0) * (E - intermediate);

        double sound_speed = std::sqrt(gamma * pressure / rho);
        double alfven_speed = std::abs(Bx) / std::sqrt(rho);
        double slow_ma_speed = std::sqrt( 0.5*(sound_speed*sound_speed + (BmagSquared/rho) - std::sqrt((sound_speed*sound_speed + BmagSquared/rho)*(sound_speed*sound_speed + BmagSquared/rho) - 4.0*sound_speed*sound_speed*Bx*Bx / rho)));
        
        double fast_ma_speed = std::sqrt( 0.5*(sound_speed*sound_speed + (BmagSquared/rho) + std::sqrt((sound_speed*sound_speed + BmagSquared/rho)*(sound_speed*sound_speed + BmagSquared/rho) - 4.0*sound_speed*sound_speed*Bx*Bx / rho)));
        
        double speed = std::abs(u_x) + fast_ma_speed;

        if (speed > maxSpeed) {
            maxSpeed = speed;
        }
    }

    double dt = C * dx / maxSpeed;
     return dt;
}
//define a function to get the flux x is ui and y is ui+1 they are both arrays of length 3
//it will spit out an array of size 3 as well give each variable their flux

StateVector  getFlux(StateVector  x , StateVector  y, double dx , double dt,double gamma){
    //impliment general flux function 
    StateVector  f_1 = FluxDef(x, gamma); //f(ui)
    StateVector  f_2 = FluxDef(y, gamma); //f(i+1)
    StateVector  uPlusHalf; //u_{i+1/2}

    for(int i=0; i<=7; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dx)*(f_2[i] - f_1[i]);
    }

    StateVector RI_flux = FluxDef(uPlusHalf, gamma); //richtmyer flux
    StateVector LF_flux; //set up 3 length array for LF flux
    StateVector FORCE_flux; //set up 3 length array for FORCE

    for (int i=0; i<=7; ++i) {
        LF_flux[i] = (0.5 * (dx/dt) * (x[i] - y[i])) + 0.5 * (f_1[i] + f_2[i]);
        FORCE_flux[i] = 0.5 * (LF_flux[i] + RI_flux[i]);
    }

    return FORCE_flux;
}

//using the minbee limiter
double minbee(double deltaMinus , double deltaPlus , double omega = 0.0){
    if (deltaPlus == 0.0) {
        return 0.0;  
    }
    double r = deltaMinus / deltaPlus;
    double xi;
    if(r <=0){
        double xi =0.0;
    }
    else if(r <=1.0){
        xi = r;
    }
    else{
        xi = std::min(1.0, 2.0/(1.0+r));
    }
    
    double Delta = 0.5 * (1.0 + omega) * deltaMinus + 0.5 * (1.0 - omega) * deltaPlus;
    
    return xi * Delta;
}

StateVector getxiDeltas(const StateVector& u_left, const StateVector& u_center,  const StateVector& u_right) {
    StateVector slopes;
    
    for (int k = 0; k < 8; k++) {
        double backward_diff = u_center[k] - u_left[k];
        double forward_diff = u_right[k] - u_center[k];
        slopes[k] = minbee(backward_diff, forward_diff);
    }
    
    return slopes;
}

// this function takes in u and spits out uBarL and uBarR
void getUbar(const std::vector<StateVector>& u, int i,  StateVector& uBarL,  StateVector& uBarR) {
    
    // Compute slopes using MinBee 
    StateVector xiDeltaL = getxiDeltas(u[i-1], u[i], u[i+1]);
    StateVector xiDeltaR = getxiDeltas(u[i], u[i+1], u[i+2]);
    
    // get the ubars but just replace the others with them
    for (int k = 0; k < 8; k++) {
        uBarL[k] = u[i][k] + 0.5 * xiDeltaL[k];      
        uBarR[k] = u[i+1][k] - 0.5 * xiDeltaR[k];  
    }
}




int main() { 
    int nCells = 800; //the distance between points is 0.01
    double x0 = 0.0;
    double x1 = 800;
    double tStart = 0.0; //set the start and finish time steps the same
    double tStop = 1.0;
    double C = 1.0;
    double gamma = 2.0;
    double omega = 0;

    // Allocate matrices with 2 extra points for transmissive BCs
    std::vector<StateVector> u(nCells+2);
    std::vector<StateVector> uPlus1(nCells+2);
    std::vector<StateVector> flux(nCells+2);
    std::vector<StateVector> uBarL(nCells+2);
    std::vector<StateVector> uBarR(nCells+2);
    std::vector<StateVector> uBarHalfL(nCells+2);
    std::vector<StateVector> uBarHalfR(nCells+2);
    double dx = (x1 - x0) / nCells; //the space steps 

    // Initial conditions!
    std::vector<double> x_vector(nCells+4);
    for(int i = 0; i < u.size(); i++) {
        // x 0 is at point i=1/2
        double x = x0 + (i-0.5) * dx;
        x_vector[i]=x;
        StateVector prim;
        if(x <= 400) {
            prim[0] = 1; // Density
            prim[1] = 0; // Velocity
            prim[2] = 0;
            prim[3] = 0;
            prim[4] = 1; // pressure
            prim[5] = 0.75; // magnetic field
            prim[6] = 1;
            prim[7] = 0; 
            } else {
            prim[0] = 0.125; // Density
            prim[1] = 0; // Velocity
            prim[2] = 0;
            prim[3] = 0;
            prim[4] = 0.1; // pressure
            prim[5] = 0.75; // magnetic field
            prim[6] = -1;
            prim[7] = 0; 
        }

        u[i] = PrimitiveToConservative(prim, gamma);
    }

    

    double dt = ComputeTimeStep(u , C , dx, gamma); //the time steps

    double t = tStart;
    int save_interval=3;
    int count=0;
    do {

        dt = ComputeTimeStep(u , C , dx, gamma); 
        t = t + dt;
        count++;
        std::cout<<"t= "<< t<< " dt= "<< dt<< std::endl;

        // Trasmissive boundary conditions
        u[0] = u[1];
        u[nCells + 1] = u[nCells];

        //find ubar


        for(int i=1; i<=nCells; ++i){
            for(int j=0; j<8; ++j){
                double DeltaPlus = u[i+1][j] - u[i][j];
                double DeltaMinus = u[i][j] - u[i-1][j];
                // std::cout<< DeltaMinus << " " << DeltaPlus << std::endl;
                double r = DeltaMinus / (DeltaPlus  + 1e-8);
                
                double xi_L = 2.0*r/(1+r);
                double xi_R = 2.0/(1+r);
                double xi;
                double Delta = 0.5*(1+omega)*DeltaMinus + 0.5*(1-omega)*DeltaPlus;
                // std::cout << Delta << std::endl;
               
                if(r<=0){ xi=0;}
                else if(r>0 && r<=1){ xi=r;}
                else{ 
                    xi=std::fmin(1, xi_R);
                    
                }

                // xi = 0.2;
                // std::cout << xi << " " << xi_R << " " << r <<  std::endl;
                uBarL[i][j] = u[i][j] - 0.5 * xi * Delta;
                uBarR[i][j] = u[i][j] + 0.5 * xi * Delta;
                
                
            }
        }

        for(int i=1; i<=nCells; ++i){
            for(int j=0; j<8; ++j){
                uBarHalfL[i][j] = uBarL[i][j] - 0.5*(dt/dx)*(FluxDef(uBarR[i] , gamma)[j]-FluxDef(uBarL[i] , gamma)[j]);
                uBarHalfR[i][j] = uBarR[i][j] - 0.5*(dt/dx)*(FluxDef(uBarR[i] , gamma)[j]-FluxDef(uBarL[i] , gamma)[j]);
                
            }
        }

        // for(int i=1; i<=nCells; ++i){
        //     auto fL = FluxDef(uBarL[i], gamma);
        //     auto fR = FluxDef(uBarR[i], gamma);
        //     for (int j = 0; j < 8; ++j) {
        //         uBarHalfL[i][j] = uBarL[i][j] - 0.5 * (dt/dx) * (fR[j] - fL[j]);
        //         uBarHalfR[i][j] = uBarR[i][j] - 0.5 * (dt/dx) * (fR[j] - fL[j]);
        //     }
        // }
        

        uBarHalfL[0] = uBarHalfL[1];
        uBarHalfL[nCells + 1] = uBarHalfL[nCells];

        uBarHalfR[0] = uBarHalfR[1];
        uBarHalfR[nCells + 1] = uBarHalfR[nCells];


        for(int i = 0; i < nCells+1; i++) { //Define the fluxes
            // flux[i] corresponds to cell i+1/2 
            flux[i] = getFlux( uBarHalfR[i], uBarHalfL[i+1] , dx , dt, gamma);
        }

        //the below has a i-1 flux which means we need to define a flux at 0 so make sure the above ^ starts at 0! this is because we have another edge with the number of cells (like the walls)

        for(int i = 1; i <= nCells+1; i++) { //Update the data
            for(int j=0; j<8; ++j){
                uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
            }
        }
    
        // Now replace u with the updated data for the next time step

        
        u = uPlus1;

        if(count%save_interval==0){
            std::string name= "Plot_" + std::to_string(t);
             //Output the data
            std::string dir = "NewSLIC/";
            save_to_file(u,x_vector,nCells,gamma,name,dir);
        }
    } while (t < tStop);


}