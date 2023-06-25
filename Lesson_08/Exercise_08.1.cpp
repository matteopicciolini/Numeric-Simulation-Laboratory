#include "Library_08.h"


int main(int argc, char *argv[]){
    Random rnd;
    Random_Start(rnd);

    int M = pow(10, 6);    
    int N = 100;            
    int L = int(M / N);

    //Progress bar
    std::string perc = "▪▪▪▪▪▪▪▪▪▪";

    //Acceptance Rate
    double accepted = 0.;
    double attempted = 0.;

    //MC 
    double delta = 2 ;
    double integral = 0;

    // Try some parameters
    double mu = 1.;
    double sigma = 0.5;         

    // Open files
    std::ofstream Integral;
    std::ofstream Coordinates;
    Integral.open(std::string(ROOT_PATH) + "/Data/08.1_integral_" + std::to_string(mu).substr(0, 3) + "_" + std::to_string(sigma).substr(0, 3) + ".dat");
    Coordinates.open(std::string(ROOT_PATH) + "/Data/08.1_coordinates_" + std::to_string(mu).substr(0, 3) + "_" + std::to_string(sigma).substr(0, 3) + ".dat");


    

    // set starting point
    double x = 0.;
    
    Equilibrate(100, pow(10, 3), x, rnd, delta, mu, sigma);

    //Blocking average variables
    double sum_prog = 0., sum_prog_2 = 0., acc = 0., acc_2 = 0., error = 0.;
    for(int j = 0; j < N; ++j){
        integral = 0;
        accepted = 0;
        attempted = 0;
        for (int i = 0; i < L; ++i){
            Metropolis(x, rnd, delta, accepted, attempted, mu, sigma);
            integral += (( -0.5 * EvalWaveFunctionSecondDerivative(x, mu, sigma) ) / EvalWaveFunction(x, mu, sigma))  + EvalPotential(x);
            Coordinates << x << std::endl;
        }
        		
		acc += integral / static_cast<double>(L);
		acc_2 += pow(integral / static_cast<double>(L), 2);
		
		sum_prog = acc / static_cast<double>(j + 1);
		sum_prog_2 = acc_2 / static_cast<double>(j + 1);
		error = Error(sum_prog, sum_prog_2, j);
        
        Integral << j << " " << sum_prog << " " << error << std::endl;
        
        //----------------------------Progress Bar------------------------------
		Progress_bar(j, N, perc);
		//----------------------------------------------------------------------
    }

    // Close files
    Integral.close();
    Coordinates.close();
    rnd.SaveSeed();
    return 0;
}