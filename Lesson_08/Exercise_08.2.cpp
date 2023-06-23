#include "Library_08.h"

int main(int argc, char *argv[]){
    

    Random rnd;
    Random_Start(rnd, 2);

    // M thows, N blocks, L throws in each block
    int M = pow(10,6);    
    int N= 1000;            
    int L = int(M / N);

    std::cout << "Working with " << L << " throws in each block" << std::endl;

    //temperature
    double T = 2.;




    //ANNEALING
    double m_beta = 1./T;
    double m_energy = 1.;
    double m_currentEnergy = m_energy;
    double m_currentError = 0.;

    double m_mu = 1.;
    double m_sigma = 1.;
    double m_position = 0;     //initial position x=0
    double m_stepLength = 2.;   //step of the proposed Metropolis move
    
    std::ofstream Averages;
    std::ofstream Coordinates;
    std::ofstream Parameters;

    Parameters.open(std::string(ROOT_PATH) + "/Data/08.2_parameters.dat");
    Averages.open(std::string(ROOT_PATH) + "/Data/08.2_averages.dat");


    Averages << T << " " << m_energy << " " << m_currentError << " " << m_energy<< std::endl;
    Parameters << T << " " << m_mu << " " << m_sigma << std::endl;

    double bestEnergy = m_energy;
    double bestParameters[2] = {m_mu, m_sigma};

    while(T >= 0.01){
        std::cout << "T: " << T << std::endl;

        //-----------------------------------------------------
        //Annealing
        double oldMu = m_mu;
        double oldSigma = m_sigma;

        //std::cout << "mu " << oldMu << std::endl;
        //std::cout << "sigma " << oldSigma << std::endl;

        m_mu = std::abs(oldMu + rnd.Rannyu(-1, 1) * .5 * (1./m_beta));
        m_sigma = std::abs(oldSigma + rnd.Rannyu(-1, 1) * .25 * (1./m_beta));


        //std::cout << "m_mu " << m_mu << std::endl;
        //std::cout << "m_sigma " << m_sigma << std::endl;

        double runningSum = 0.;
        double runningSquared = 0.;
        double error = 0.;
        double integral = 0.;

        EquilibrateUN(100, pow(10, 3), m_position, rnd, m_stepLength, m_mu, m_sigma);

        double sum_prog = 0., sum_prog_2 = 0., acc = 0., acc_2 = 0.;
        for(int j = 0; j < N; ++j){
            integral = 0;
            double accepted = 0;
            double attempted = 0;
            for (int i = 0; i < L; ++i){
                MetropolisUniform(m_position, rnd, m_stepLength, accepted, attempted, m_mu, m_sigma);
                integral += (( -0.5 * EvalWaveFunctionSecondDerivative(m_position, m_mu, m_sigma) ) / EvalWaveFunctionNoAbs(m_position, m_mu, m_sigma))  + EvalPotential(m_position);
            }

            acc += integral / static_cast<double>(L);
            acc_2 += pow(integral / static_cast<double>(L), 2);
            
            sum_prog = acc / static_cast<double>(j + 1);
            sum_prog_2 = acc_2 / static_cast<double>(j + 1);
            error = Error(sum_prog, sum_prog_2, j);
            
        }

        m_currentEnergy = sum_prog;
        m_currentError = error;

        double alpha = std::min(1., (exp(- m_beta * (m_currentEnergy - m_energy))));
        //std::cout << alpha <<std::endl;
        double p = rnd.Rannyu();

        if (alpha > p){
            m_energy = sum_prog;
        }

        else{
            m_mu = oldMu;
            m_sigma = oldSigma;
        }

        //------------------------------------------------------------------

        Averages << T << " " << integral << " " << m_energy << " " << m_currentError << std::endl;
        Parameters << T << " " << m_mu << " " << m_sigma << std::endl;

        T = T * 0.997;
        m_beta = 1/T;
        if (bestEnergy > m_energy) {
            bestEnergy = m_energy;
            bestParameters[0] = m_mu;
            bestParameters[1] = m_sigma;
        }
    }

    std::cout << "Lowest energy found: " << bestEnergy << std::endl;
    std::cout << "Optimized parameters: mu = " << bestParameters[0] << " sigma = " << bestParameters[1] << std::endl;

    Averages.close();
    Parameters.close();

    Averages.open(std::string(ROOT_PATH) + "/Data/08.2_optimizedEnergy.dat");
    Coordinates.open(std::string(ROOT_PATH) + "/Data/08.2_optimizedCoordinates.dat");

    std::cout << "Finding the optimized energy" << std::endl;


    m_mu = bestParameters[0];
    m_sigma = bestParameters[1];

    double runningSum = 0.;
    double runningSquared = 0.;
    double error = 0.;

    EquilibrateUN(100, pow(10, 3), m_position, rnd, m_stepLength, m_mu, m_sigma);

    double sum_prog = 0., sum_prog_2 = 0., acc = 0., acc_2 = 0.;
    for (int j = 0; j < N; ++j){
        double integral = 0;
        double accepted = 0;
        double attempted = 0;
        for(int i = 0; i < L; ++i){
            MetropolisUniform(m_position, rnd, m_stepLength, accepted, attempted, m_mu, m_sigma);
            integral += (( -0.5 * EvalWaveFunctionSecondDerivative(m_position, m_mu, m_sigma) ) / EvalWaveFunctionNoAbs(m_position, m_mu, m_sigma))  + EvalPotential(m_position);
            Coordinates << m_position << std::endl;
        }
        acc += integral / static_cast<double>(L);
		acc_2 += pow(integral / static_cast<double>(L), 2);
		
		sum_prog = acc / static_cast<double>(j + 1);
		sum_prog_2 = acc_2 / static_cast<double>(j + 1);
		error = Error(sum_prog, sum_prog_2, j);
        
        Averages << j << " "<< sum_prog << " " << error << std::endl;
    }

    Averages.close();
    Coordinates.close();

    rnd.SaveSeed();
    return 0;
}
