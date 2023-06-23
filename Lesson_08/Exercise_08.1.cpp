#include "Library_08.h"


int main(int argc, char *argv[]){
    Random rnd;
    Random_Start(rnd);

    int M = pow(10, 6);    
    int N = 100;            
    int L = int(M / N);

    std::cout << "Working with " << L << " throws in each block" << std::endl;

    //variables for acceptance rate of the algorithm
    double accepted = 0.;
    double attempted = 0.;

    //Define global variables
    double a0 = 0.0529 * pow(10, -9); // reduced units
    double c = 2 ;       //to tune the acceptance rate
    double integral = 0;             

    //Blocking average variables
    double runningSum = 0.;
    double runningSquared = 0.;
    double error = 0.;

    std::ofstream Averages;
    std::ofstream Coordinates;

    double mu = 1;
    double sigma = 0.5;

    double x = 0.;
    Averages.open(std::string(ROOT_PATH) + "/Data/08.1_averages_" + std::to_string(mu).substr(0, 3) + "_" + std::to_string(sigma).substr(0, 3) + ".dat");
    Coordinates.open(std::string(ROOT_PATH) + "/Data/08.1_coordinates_" + std::to_string(mu).substr(0, 3) + "_" + std::to_string(sigma).substr(0, 3) + ".dat");

    EquilibrateUN(100, pow(10,4), x, rnd, c, mu, sigma);
    double sum_prog = 0., sum_prog_2 = 0., acc = 0., acc_2 = 0.;
    for(int j = 0; j < N; ++j){
        integral = 0;
        accepted = 0;
        attempted = 0;
        for (int i = 0; i < L; ++i){
            MetropolisUniform(x, rnd, c, accepted, attempted, mu, sigma);
            integral += (( -0.5 * EvalWaveFunctionSecondDerivative(x, mu, sigma) ) / EvalWaveFunctionNoAbs(x, mu, sigma))  + EvalPotential(x);
            Coordinates << x << std::endl;
        }
        		
		acc += integral / static_cast<double>(L);
		acc_2 += pow(integral / static_cast<double>(L), 2);
		
		sum_prog = acc / static_cast<double>(j + 1);
		sum_prog_2 = acc_2 / static_cast<double>(j + 1);
		error = Error(sum_prog, sum_prog_2, j);
        
        Averages << j << " "<< sum_prog << " " << error << std::endl;
        
        if (j%10 == 0){
           std::cout <<"Block: "<< j << " Uniform Acceptance rate " << accepted / attempted << std::endl;
        }
    }
    Averages.close();
    Coordinates.close();
    rnd.SaveSeed();
    return 0;
}