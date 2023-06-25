#include "Library_08.h"

int main(int argc, char *argv[]){
    Random rnd;
    Random_Start(rnd, 2);

    // M thows, N blocks, L throws in each block
    int M = pow(10,6);    
    int N = 1000;            
    int L = int(M / N);

    //temperature
    double T = 2.;

    //ANNEALING
    double beta = 1./T;
    double energy = 1.;
    double currentEnergy = energy;
    double currentError = 0.;

    double mu = 1.;
    double sigma = 1.;
    double position = 0;     //initial position x=0
    double delta = 2.;   //step of the proposed Metropolis move
    
    std::ofstream Integral;
    std::ofstream Coordinates;
    std::ofstream Parameters;
    Parameters.open(std::string(ROOT_PATH) + "/Data/08.2_parameters.dat");
    Integral.open(std::string(ROOT_PATH) + "/Data/08.2_integral.dat");


    Integral << T << " " << energy << " " << currentError << " " << energy<< std::endl;
    Parameters << T << " " << mu << " " << sigma << std::endl;

    double bestEnergy = energy;
    double bestParameters[2] = {mu, sigma};

    while(T >= 0.01){
        std::cout << "Temperature: " << T << std::endl;

        //-----------------------------------------------------
        //Annealing
        double oldMu = mu;
        double oldSigma = sigma;

        mu = std::abs(oldMu + rnd.Rannyu(-1, 1) * .5 * (1./beta));
        sigma = std::abs(oldSigma + rnd.Rannyu(-1, 1) * .25 * (1./beta));

        Equilibrate(100, pow(10, 3), position, rnd, delta, mu, sigma);

        double sum_prog = 0., sum_prog_2 = 0., acc = 0., acc_2 = 0., error = 0., integral = 0.;
        for(int j = 0; j < N; ++j){
            integral = 0;
            double accepted = 0;
            double attempted = 0;
            for (int i = 0; i < L; ++i){
                Metropolis(position, rnd, delta, accepted, attempted, mu, sigma);
                integral += (( -0.5 * EvalWaveFunctionSecondDerivative(position, mu, sigma) ) / EvalWaveFunction(position, mu, sigma))  + EvalPotential(position);
            }

            acc += integral / static_cast<double>(L);
            acc_2 += pow(integral / static_cast<double>(L), 2);
            
            sum_prog = acc / static_cast<double>(j + 1);
            sum_prog_2 = acc_2 / static_cast<double>(j + 1);
            error = Error(sum_prog, sum_prog_2, j);
        }

        currentEnergy = sum_prog;
        currentError = error;

        double alpha = std::min(1., (exp(- beta * (currentEnergy - energy))));
        double p = rnd.Rannyu();

        if (alpha > p)
            energy = sum_prog;
        else{
            mu = oldMu;
            sigma = oldSigma;
        }

        //------------------------------------------------------------------
        
        Integral << T << " " << integral << " " << energy << " " << currentError << std::endl;
        Parameters << T << " " << mu << " " << sigma << std::endl;

        T = T * 0.997;
        beta = 1/T;
        if (bestEnergy > energy) {
            bestEnergy = energy;
            bestParameters[0] = mu;
            bestParameters[1] = sigma;
        }
    }

    std::cout << "Lowest energy found: " << bestEnergy << std::endl;
    std::cout << "Optimized parameters: mu = " << bestParameters[0] << " sigma = " << bestParameters[1] << std::endl;

    Integral.close();
    Parameters.close();

    //Exercise_08.1 with best parameters

    Integral.open(std::string(ROOT_PATH) + "/Data/08.2_optimizedEnergy.dat");
    Coordinates.open(std::string(ROOT_PATH) + "/Data/08.2_optimizedCoordinates.dat");

    std::cout << "Finding the optimized energy" << std::endl;

    //Progress bar
    std::string perc = "▪▪▪▪▪▪▪▪▪▪";

    //Acceptance Rate
    double accepted = 0.;
    double attempted = 0.;

    mu = bestParameters[0];
    sigma = bestParameters[1];

    Equilibrate(100, pow(10, 3), position, rnd, delta, mu, sigma);
    double sum_prog = 0., sum_prog_2 = 0., acc = 0., acc_2 = 0., error = 0., integral = 0;
    for (int j = 0; j < N; ++j){
        integral = 0;
        accepted = 0;
        attempted = 0;
        for(int i = 0; i < L; ++i){
            Metropolis(position, rnd, delta, accepted, attempted, mu, sigma);
            integral += (( -0.5 * EvalWaveFunctionSecondDerivative(position, mu, sigma) ) / EvalWaveFunction(position, mu, sigma))  + EvalPotential(position);
            Coordinates << position << std::endl;
        }
        acc += integral / static_cast<double>(L);
		acc_2 += pow(integral / static_cast<double>(L), 2);
		
		sum_prog = acc / static_cast<double>(j + 1);
		sum_prog_2 = acc_2 / static_cast<double>(j + 1);
		error = Error(sum_prog, sum_prog_2, j);
        
        Integral << j << " "<< sum_prog << " " << error << std::endl;

        //----------------------------Progress Bar------------------------------
		Progress_bar(j, N, perc);
		//----------------------------------------------------------------------
    }

    Integral.close();
    Coordinates.close();
    rnd.SaveSeed();
    return 0;
}
