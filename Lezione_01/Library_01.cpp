#include "Library_01.h"

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

void average_calc(int n_iteration, int n_block, std::vector<double> &ave, std::vector<double> &ave_2, Random &random_generator){
	int L = n_iteration / n_block;
	double sum = 0., rand;

	for (int i = 0; i < n_block; ++i){
		sum = 0.;
		for (int j = 0; j < L; ++j){
			rand = random_generator.Rannyu();
			sum += rand;
		}
		ave.push_back(sum/static_cast<double>(L));
		ave_2.push_back(pow(ave[i], 2));
	}
}

void blocks_averege(int n_iteration, int n_block, std::vector<double> &ave, std::vector<double> &ave_2, std::ofstream &file_output){
	double sum_prog = 0., sum2_prog = 0.; //progressive sum
	for (int i = 0; i < n_block; ++i){
        sum_prog = 0.;
        sum2_prog = 0.;
		for (int j = 0; j < i + 1; ++j){
			sum_prog += ave[j];
			sum2_prog += ave_2[j];
		}
		sum_prog /= i + 1;
		sum2_prog /= i + 1;
		file_output << (i + 1) << std::setw(20) << sum_prog << std::setw(20) << error(sum_prog, sum2_prog, i) << std::endl;
	}
}

void characteristic_function_simple(double &sum, Random &random_generator){
	sum += random_generator.Rannyu();
}

void blocks_averege(int n_iteration, int n_block, Random &random_generator, std::ofstream &file_output, void (*p_to_function)(double&, Random &)){
	int L = n_iteration / n_block;
	double sum = 0., ave_i = 0., ave2_i = 0.;
	double acc = 0., acc_2 = 0.;
	double sum_prog, su2_prog, err_prog;

	for(int i = 0; i < n_block; ++i){
		sum = 0.;
		for (int j = 0; j < L; ++j){
			p_to_function(sum, random_generator);
		}
		ave_i = sum / static_cast<double>(L);
		ave2_i = pow(ave_i, 2);
		
		acc += ave_i;
		acc_2 += ave2_i;
		
		sum_prog = acc / static_cast<double>(i + 1);
		su2_prog = acc_2 / static_cast<double>(i + 1);
		err_prog = error(sum_prog, su2_prog, i);
		file_output << (i+1) << std::setw(20) << sum_prog << std::setw(20) << err_prog << std::endl;
	}
}

void characteristic_function_sigma(double &sum, Random &random_generator){
    sum += pow(random_generator.Rannyu() - 0.5, 2);
}

double exponential_distribution(double lambda, double random){
    return - 1./lambda * std::log(1 - random);
}

double Lorentz_distribution(double mean, double gamma, double random){
    return mean + gamma * tan(M_PI * (random - 0.5));
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
		file_output << (i+1) << std::setw(20) << sum_prog << std::setw(20) << err_prog << std::endl;
	}
}

void characteristic_function_buffon(double& nhint, Random & random_generator){
    double L = 0.8;
    int N_HINT = static_cast<int>(nhint);
    double x, y;
    double y_0 = random_generator.Rannyu(0, 1);
    do{
        x = random_generator.Rannyu();
        y = random_generator.Rannyu(-1.,1.);
        
    }while((pow(x,2) + pow(y,2)) > 1.);
    
    double cos_theta = y / sqrt((pow(x,2) + pow(y,2)));
    
    if(std::floor(y_0 + (L/2. * cos_theta)) != std::floor(y_0 - (L/2. * cos_theta))){
        N_HINT++;
    }
    nhint = static_cast<double>(N_HINT);
}

void characteristic_function_buffon_2(double& ave_i, double& n_hint, int& L){
    ave_i =  2 * 0.8 * static_cast<double>(L)/ n_hint;
    n_hint = 0;
}