#include "Library_01.h"


int main (int argc, char *argv[]){

	//Random Start
	Random rnd;
	Random_Start(rnd);

	double ave_exp = 0., sum_exp = 0.;
	double ave_std = 0., sum_std = 0.;
	double ave_lrt = 0., sum_lrt = 0.;

	//parameters
	double lambda = 1.;
	double mean = 0., gamma = 1.;

	int N = 10000, sum_dim = 0;
	int dim[4]= {1, 2, 10, 100};
	double rand;

	std::ofstream file_out_exp("../../Data/01.2_exponential.dat");
	std::ofstream file_out_std("../../Data/01.2_standard.dat");
	std::ofstream file_out_lrt("../../Data/01.2_lorentzian.dat");

	for(int i = 0; i < N; ++i){
		sum_exp = 0.;
		sum_std = 0.;
		sum_lrt = 0.;
		sum_dim = 0;

		// in this way for each loop I avoid redoing the sum, 
		// but I use the sum I already had before
		for(int j = 0; j < 4; ++j){
			int D = dim[j] - sum_dim;
			for(int k = 0; k < D; ++k){ 
				rand = rnd.Rannyu();
				sum_exp += exponential_distribution(lambda, rand);
				sum_lrt += Lorentz_distribution(mean, gamma, rand);
				sum_std += rand;
				
				sum_dim++;
			}
			
			ave_exp = sum_exp / dim[j];
			ave_std = sum_std / dim[j];
			ave_lrt = sum_lrt / dim[j];

			file_out_exp << ave_exp;
			file_out_std << ave_std;
			file_out_lrt << ave_lrt;

			if(j!=3){
				file_out_exp << std::setw(20);
				file_out_std << std::setw(20);
				file_out_lrt << std::setw(20);
			}
		}
		file_out_exp << std::endl;
		file_out_std << std::endl;
		file_out_lrt << std::endl;
	}
	file_out_exp.close();
	file_out_std.close();
	file_out_lrt.close();
	
	
	//Random Save Seed
    rnd.SaveSeed();

	return 0;
}