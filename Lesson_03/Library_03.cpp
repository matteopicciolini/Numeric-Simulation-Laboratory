#include "Library_03.h"

void Random_Start(Random &random_generator){

	int seed[4];
	int p1, p2;
	std::ifstream Primes(std::string(ROOT_PATH) + "/random-library/Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
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

double error(double AV, double AV2, int n){
	if (n == 0){
		return 0;
	}
	else {
		return sqrt((AV2 - pow(AV,2)) / static_cast<double>(n));
	}
}
