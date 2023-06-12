#include "Library_08.h"

int main(int argc, char* argv[]){

	Random_Start(rnd);
    
    file_out_smt.open(std::string(ROOT_PATH) + "/Data/08.1/smt.dat");

    //Equilibration
    for (int j = 0; j < 10000; ++j) {
        Move();
        Energy();
    }

    std::cout << "Acceptance rate: " << (double) i / accepted << std::endl;

    Reset();
    for(int j = 0; j < N; ++j){
        blk_ave = 0;
        for(int h = 0; h < L; ++h){
            Move();
            Energy();
        }
        sum_prog += blk_ave / L;
        sum_prog_2 += std::pow(blk_ave / L, 2);
        if(j == 0){
            file_out_smt << j + 1 << " " << sum_prog / (j + 1) << " " << 0. << std::endl;
        }
        else{
            file_out_smt << j + 1 << " " << sum_prog/(j + 1) << " " << std::sqrt((std::abs(sum_prog_2 / (j + 1) - std::pow(sum_prog / (j + 1), 2))) / j) << std::endl;
        }
        mu_old = mu;
        sigma_old = sigma;
        energy_old = sum_prog/N;
        sigma_energy_old = std::sqrt((std::abs(sum_prog_2 / N - std::pow(sum_prog / N, 2))) / (N - 1));  
    }


    rnd.SaveSeed();
    return 0;
}

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

void Move(){
    x = x_old + rnd.Rannyu(-step, step);

    double p_old = std::pow(Psi(x_old),2);
    double p_new = std::pow(Psi(x), 2);

    A = std::min(1, p_new / p_old);
    if (rnd.Rannyu() < A) {
        x_old = x;
        i++;
        
    }
    accepted++;
}

double Psi(double x){
    double alpha_1, alpha_2;
    alpha_1 = 0.5 * std::pow(x - mu, 2) / std::pow(sigma, 2);
    alpha_2 = 0.5 * std::pow(x + mu, 2) / std::pow(sigma, 2);
    return std::exp(-alpha_1) + std::exp(-alpha_2);
}

void Energy(){
    double alpha1 = pow(x_old - mu, 2)/pow(sigma, 2) * 0.5;
    double alpha2 = pow(x_old + mu, 2)/pow(sigma, 2) * 0.5;
    double a_1 = pow(mu, 2) - pow(sigma, 2) + pow(x_old, 2) - 2 * mu * x_old;
    double a_2 = pow(mu, 2) - pow(sigma, 2) + pow(x_old, 2) + 2 * mu * x_old;
    double v = pow(x_old, 4) - 5./2. * pow(x_old, 2);
    
    E = -0.5/(psi(x_old) * pow(sigma,4)) * (exp(-alpha1) * a_1 + exp(-alpha2) * a_2) + v;
    
    blk_ave += E;
}