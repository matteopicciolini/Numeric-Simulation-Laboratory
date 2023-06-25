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

double EvalWaveFunctionSquared(double x, double mu, double sigma){
    return std::pow(std::abs(exp(-std::pow(x - mu, 2) / (2 * std::pow(sigma, 2))) + exp(-std::pow(x + mu, 2) / (2 * std::pow(sigma, 2)))), 2);
}

void Equilibrate(int nblocks, int L, double &position, Random &rnd, double delta, double mu, double sigma){
    double accepted = 0.;
    double attempted = 0.;
    for (int j = 0; j < nblocks; j++){
        for (int i = 0; i < L; i++){
            Metropolis(position, rnd, delta, accepted, attempted, mu, sigma);
        }
    }
}

void Metropolis(double &position, Random &rnd, double delta, double &accepted, double &attempted, double mu, double sigma){

    double future_position = position + rnd.Rannyu(-1, 1) * delta;
    double alpha = std::min(1., (EvalWaveFunctionSquared(future_position, mu, sigma) / EvalWaveFunctionSquared(position, mu, sigma)));

    double p = rnd.Rannyu();
    if (p < alpha){
        position = future_position;
        accepted++;
    }
    attempted++;
}

double EvalWaveFunctionSecondDerivative(double x, double mu, double sigma){

    double minusExp = exp(-0.5 * (std::pow(x - mu, 2) / (std::pow(sigma, 2))));
    double plusExp = exp(-0.5 * (std::pow(x + mu, 2) / (std::pow(sigma, 2))));

    return ((-1 / std::pow(sigma, 2)) * minusExp) + ((-1 / std::pow(sigma, 2)) * plusExp) + ((std::pow(x - mu, 2) / std::pow(sigma, 4)) * minusExp) + ((std::pow(x + mu, 2) / std::pow(sigma, 4)) * plusExp);
}

double EvalWaveFunction(double x, double mu, double sigma){
    return exp(-std::pow(x - mu, 2) / (2 * std::pow(sigma, 2))) + exp(-std::pow(x + mu, 2) / (2 * std::pow(sigma, 2)));
}

double Error(double AV, double AV2, int n){
	if (n == 0){
		return 0;
	}
	else {
		return sqrt((AV2 - pow(AV,2)) / static_cast<double>(n));
	}
}

//-----------------------Progress bar ----------------------------------

std::string Green(std::string green){
	return std::string("\033[1;32m") + green + "\033[0m";
}

std::string Gray(std::string gray){
	return std::string("\033[1;90m") + gray + "\033[0m";
}

std::string Red(std::string red){
	return std::string("\033[1;31m") + red + "\033[0m";
}

void Progress_bar(int& prog, int N, std::string perc){
	std::cout << Red("Block number " + std::to_string(prog) + ". Progress: ");
	std::cout << Green(perc.substr(0, 3 * std::floor(prog * 10.0 / N)));
	std::cout << Gray(perc.substr(3 * std::floor(prog * 10.0 / N), 3 * 10));
	std::cout << " " << int(prog * 100.0 / N) << " %\r";
	std::cout.flush();
}