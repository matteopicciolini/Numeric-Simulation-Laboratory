#include "Library_08.h"

void Random_Start(Random &random_generator, int primes_selector){
	int seed[4];
	int p1, p2;
	std::ifstream Primes(std::string(ROOT_PATH) + "/random-library/Primes");
	if (Primes.is_open()){
		for (int i = 0; i < 1 + primes_selector; ++i) Primes >> p1 >> p2 ;
	} else {std::cerr << "PROBLEM: Unable to open Primes" << std::endl;}
	
	Primes.close();
	
	std::ifstream input(std::string(ROOT_PATH) + "/random-library/seed.in");
	std::string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				random_generator.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else {std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;}
}

double EvalPotential(double x, double a, double b){
    return a * std::pow(x, 2) + b * std::pow(x, 4);
}

double EvalWaveFunction(double x, double mu, double sigma){
    return std::pow(std::abs(exp(-std::pow(x - mu, 2) / (2 * std::pow(sigma, 2))) + exp(-std::pow(x + mu, 2) / (2 * std::pow(sigma, 2)))), 2);
}

void EquilibrateUN(int nblocks, int L, double &initialPosition, Random &rnd, double c, double mu, double sigma){
    double accepted = 0.;
    double attempted = 0.;
    for (int j = 0; j < nblocks; j++){
        for (int i = 0; i < L; i++){
            MetropolisUniform(initialPosition, rnd, c, accepted, attempted, mu, sigma);
        }
    }
}

void MetropolisUniform(double &initialPosition, Random &rnd, double c, double &accepted, double &attempted, double mu, double sigma){

    double tempPosition = initialPosition + rnd.Rannyu(-1, 1) * c;
    double alpha = std::min(1., (EvalWaveFunction(tempPosition, mu, sigma) / EvalWaveFunction(initialPosition, mu, sigma)));

    // accepting the configuration with probability \alpha:
    double p = rnd.Rannyu();
    if (p < alpha){
        initialPosition = tempPosition;
        accepted++;
    }
    attempted++;
}

double EvalWaveFunctionSecondDerivative(double x, double mu, double sigma){

    double minusExp = exp(-0.5 * (std::pow(x - mu, 2) / (std::pow(sigma, 2))));
    double plusExp = exp(-0.5 * (std::pow(x + mu, 2) / (std::pow(sigma, 2))));

    return ((-1 / std::pow(sigma, 2)) * minusExp) + ((-1 / std::pow(sigma, 2)) * plusExp) + ((std::pow(x - mu, 2) / std::pow(sigma, 4)) * minusExp) + ((std::pow(x + mu, 2) / std::pow(sigma, 4)) * plusExp);
}

double EvalWaveFunctionNoAbs(double x, double mu, double sigma){
    return std::abs(exp(-std::pow(x - mu, 2) / (2 * std::pow(sigma, 2))) + exp(-std::pow(x + mu, 2) / (2 * std::pow(sigma, 2))));
}

//double Error(double sum, double sum2, int iblk){
//    if(iblk == 1) return 0.0;
//    else return sqrt((sum2/(double)iblk - std::pow(sum/(double)iblk,2))/(double)(iblk-1));
//}

double Error(double AV, double AV2, int n){
	if (n == 0){
		return 0;
	}
	else {
		return sqrt((AV2 - pow(AV,2)) / static_cast<double>(n));
	}
}