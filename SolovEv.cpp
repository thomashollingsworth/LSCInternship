#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <filesystem>  // C++17
#include <initializer_list>

#include <Eigen/Sparse>
using namespace Eigen; //Preferentially uses Eigen library functions in case of conflicts
typedef Triplet<double> T; //Shorthand for the eigen Triplet class which stores data as (row,column,value) where the data is of double type


std::tuple< std::vector<double>, double> get_coordinate_vector(const int& num_points, const double& point_min, const double& point_max){
    std::vector<double> output(num_points);
    double dX= (point_max-point_min)/num_points;
    for(int i=0;i<num_points;i++){
        output[i]=point_min+(i-0.5)*dX;
    };
    return {output,dX};

}
int index(const int& i, const int& j, const int& num_R){
    //Returns the 1D index for a 2D array defined with row i in Z direction and column j for R direction
    return i + j*num_R;
}

std::tuple<int, int> indicesFromRow(const int& row, const int& num_R){
    int j = std::floor(row/num_R);
    int i = row - j*num_R;
    return {i,j};
}


//Functions for calculating matrix elements of the Toroidal Elliptical Operator
double D1(const double& R, const double& dR ){
    return 1./(dR*dR)-1/(2*R*dR);
}
double D2(const double& R, const double& dR ){
    return 1./(dR*dR)+1/(2*R*dR);
}
double D3(const double& dZ ){
    return 1./(dZ*dZ);
}
double D4(const double& dZ,const double& dR ){
    return 1./(dZ*dZ)+1./(dR*dR);
}


SparseMatrix<double> getTEOperator(const std::vector<double>& R_vector,const std::vector<double>& Z_vector, const double& dR, const double& dZ ){
    //Returns a sparse matrix for performing the toroidal elliptic operator in 2D over a cartesian grid (num_R x num_Z)
    //Matrix will have dimensions (N x N) where N is num_R * num_Z
    //Matrix is constructed using Triplets
    const double num_R=R_vector.size();
    const double num_Z=Z_vector.size();

    const double N= num_R*num_Z;
    
    std::vector<T> triplets; //A vector of triplet types
    triplets.reserve(N * 5);  // ~5 non-zeros per row

    for(int j=0;j<num_Z;j++){
        for(int i=0;i<num_R;i++){
            int row = index(i, j,num_R);//gives the row index for a (R,Z) point (i,j)
            if (i == 0 || i == num_R-1 || j == 0 || j == num_Z-1){
                //If on the boundary the row should just be taken from the identity matrix
                triplets.emplace_back(row,row,1.);
                continue;
            }
            //Set standard matrix elements
            triplets.emplace_back(row,index(i+1,j,num_R),D1(R_vector[i],dR));
            triplets.emplace_back(row,index(i-1,j,num_R),D2(R_vector[i],dR));
            triplets.emplace_back(row,index(i,j+1,num_R),D3(dZ));
            triplets.emplace_back(row,index(i,j-1,num_R),D3(dZ));
            triplets.emplace_back(row,row,-2*D4(dZ,dR));

        }
    }
    SparseMatrix<double> A(N, N);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

double getAnalyticVal(const double& R, const double& Z, const double& alpha, const double& beta){
    return 0.5*beta*R*R*Z*Z +1./8.*(alpha-beta)*(R*R-1)*(R*R-1);

}
std::vector<double> getAnalyticVector(const std::vector<double>& R_vector,const std::vector<double>& Z_vector, double alpha,double beta){
    //Returns a 1D array of the results
    int numR=R_vector.size();
    int numZ=Z_vector.size();
    int N=numZ*numR;
    std::vector<double> results(N);

    for(int j=0;j<numZ;j++){
        for(int i=0;i<numR;i++){
            results[index(i,j,numR)]= getAnalyticVal(R_vector[i],Z_vector[j],alpha,beta);
    }
    
}
return results;
}

double getRHSVal(const double& R,const double& alpha){
    //gives the elements for the right hand side of the linear non-dimensional Solov'Ev equation 
    return alpha*R*R;
}
VectorXd getRHSVector(const std::vector<double>& R_vector,const std::vector<double>& Z_vector, double alpha,double beta){
    //gives the vector for the right hand side of the linear non-dimensional Solov'Ev equation
    int numR=R_vector.size();
    int numZ=Z_vector.size();
    int N=numZ*numR;

    VectorXd results(N);

    for(int j=0;j<numZ;j++){
        for(int i=0;i<numR;i++){
            if (i == 0 || i == numR-1 || j == 0 || j == numZ-1){
                //If on the boundary the row should just be taken from the identity matrix
                results[index(i,j,numR)]=getAnalyticVal(R_vector[i],Z_vector[j],alpha,beta);
                continue;
            }
            results[index(i,j,numR)]= getRHSVal(R_vector[i],alpha);
    }}
    return results;
}

void save_to_file(const std::vector<double>& R_vector,const std::vector<double>& Z_vector, const std::vector<double>& x, double numR, double numZ, std::string dir, std::string filename,double kappa){
    //Convert to primitive variables
    double N=numR*numZ;
    std::vector< std::vector<double> > results_array(numR,std::vector<double>(numZ));

    //Output the data
    std::filesystem::create_directories(dir);
    
    filename=dir+filename;

    std::ofstream outfile(filename);
    
    for(size_t j=0; j<numZ;j++){
        for(size_t i=0; i<numR;i++){
        
            outfile << R_vector[i] << " " << Z_vector[j] << " " << x[index(i,j,numR)] << std::endl;
}
outfile<<std::endl;; //Blank line between y rows for correct gnuplot formatting


}

    std::string gnuplot_cmd = "gnuplot -e \"filename='" + filename + "'; kappa=" + std::to_string(kappa) + "\" SolovEv_plot_script.plt";;
    std::system(gnuplot_cmd.c_str());
    std::cout << "Plot saved to " <<filename<< ".png" << std::endl;
    //std::filesystem::remove(filename);
}




int main(void){

    const double kappa = 1.7; //Elongation of the plasma i.e. numerical ratio of ellipse axes vertical/horizontal
    const double q_0 =1; //Safety factor at magnetic axis

    const int numR=200; //Number of points in R direction
    const int numZ=200; //Number of points in Z direction

    int N=numR*numZ;

    //Quantities measured in dimensionless units i.e. Z, R measured in R_0
    double R_min =0.6;
    double R_max=1.4;
 
    double Z_min =-kappa*(R_max-R_min)/2.;
    double Z_max=+kappa*(R_max-R_min)/2.;


    double alpha = (kappa*kappa+1)/(kappa*q_0); //Sole parameter in Dimensionless Solov'ev equation
    double beta= 1./(kappa*q_0); //Additional parameter used for analytic solution

    auto [R_vector, dR] =get_coordinate_vector(numR,R_min,R_max);
    auto [Z_vector, dZ] =get_coordinate_vector(numZ,Z_min,Z_max);

    SparseLU<SparseMatrix<double>> LUsolver;
    SparseMatrix<double> A = getTEOperator(R_vector,Z_vector,dR,dZ);
    VectorXd b = getRHSVector(R_vector,Z_vector,alpha,beta);
    VectorXd x = VectorXd::Zero(N);//Stores the output
    
    LUsolver.compute(A);
    if(LUsolver.info() != Eigen::Success) {
    std::cerr << "Decomposition failed!" << std::endl;
    // Handle error
}
    x = LUsolver.solve(b);

    if(LUsolver.info() != Eigen::Success) {
    std::cerr << "Solving failed!" << std::endl;
    // Handle error
}
    std::vector<double> x_solution(x.data(), x.data() + x.size());
    std::vector<double> x_analytic=getAnalyticVector(R_vector,Z_vector,alpha,beta);
  

    std::string dir = "SolovEvPlots/";

    save_to_file(R_vector,Z_vector,x_solution,numR,numZ,dir,"solution",kappa);
    save_to_file(R_vector,Z_vector,x_analytic,numR,numZ,dir,"analytic",kappa);
    
}


    










