#include "Library_02.h"

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

void blocks_averege(int n_iteration, int n_block, Random &random_generator, std::ofstream &file_output, void (*p_to_function)(double&, Random &), void (*p_to_function_2)(double&, double&, int&)){
	int L = n_iteration / n_block;
	double sum = 0., ave_i = 0., ave2_i = 0.;
	double acc = 0., acc_2 = 0.;
	double sum_prog, su2_prog, err_prog;

	for(int i = 0; i < n_block; ++i){
		sum = 0.;
		for (int j = 0; j < L; ++j){
			p_to_function(sum, random_generator);
		}
        p_to_function_2(ave_i, sum, L);
		ave2_i = pow(ave_i, 2);
		
		acc += ave_i;
		acc_2 += ave2_i;
		
		sum_prog = acc / static_cast<double>(i + 1);
		su2_prog = acc_2 / static_cast<double>(i + 1);
		err_prog = error(sum_prog, su2_prog, i);
		file_output << (i+1)*L << std::setw(20) << sum_prog << std::setw(20) << err_prog << std::endl;
	}
}

void characteristic_function_uniform_1(double& sum, Random & random_generator){
    double x = 0.;
    x = random_generator.Rannyu();
    sum += 0.5 * M_PI * cos(M_PI * 0.5 * x);
}

void characteristic_function_uniform_2(double& ave_i, double& sum, int& L){
    double xmax = 1.;
    double xmin = 0.;
    ave_i = (xmax - xmin) * sum / static_cast<double>(L);
}

void characteristic_function_sampling_1(double& sum, Random & random_generator){
    double x = 0.;
    x = 1 - pow(1 - random_generator.Rannyu(), 0.5);
    sum += 0.5 * M_PI * cos(M_PI * 0.5 * x) / (2 * (1 - x));
}

void random_walk_discrete_step(Random &random_generator, std::vector<double> &cartesian_position, double step_lenght){
	int x = random_generator.Rannyu(0,3);
	double var_sign = random_generator.Rannyu(-1,1);
	var_sign = (var_sign >= 0) - (var_sign < 0); //boolean value is 0 when false, 1 when true
	cartesian_position[x] += var_sign * step_lenght;
}
void random_walk_continue_step(Random &random_generator, std::vector<double> &cartesian_position, double step_lenght){
	// Generate three random numbers between -1 and 1

	double theta = random_generator.Rannyu(0, 2.* M_PI);
	double u = random_generator.Rannyu(-1, 1);
	cartesian_position[0] += step_lenght * sqrt(1 - pow(u, 2)) * cos(theta);
	cartesian_position[1] += step_lenght * sqrt(1 - pow(u, 2)) * sin(theta);
	cartesian_position[2] += step_lenght * u;
	
	/*
	double theta = random_generator.Rannyu(0, 2.* M_PI);
	double phi = random_generator.Rannyu(0, M_PI);
	cartesian_position[0] += step_lenght * sin(phi) * cos(theta);
	cartesian_position[1] += step_lenght * sin(phi) * sin(theta);
	cartesian_position[2] += step_lenght * cos(phi);*/
}
