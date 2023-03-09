#include "Library_01.h"


int main (int argc, char *argv[]){

    //Random Start
	Random rnd;
	Random_Start(rnd);
    
    //-------------------- PART 1 --------------------

    //Method 1 - Python mode (2 iteration needed to obtain 1 vector of blocks average)
	std::vector<double> ave, ave_2;
	average_calc(1e5, 100, ave, ave_2, rnd);

	std::ofstream file_out(std::string(ROOT_PATH) + "/Data/01.1a_mean_value_random_generator_M1.dat");
	blocks_averege(1e5, 100, ave, ave_2, file_out);
	file_out.close();
    
    //Method 2 - Pointer function mode (1 iteration needed to obtain 1 vector of blocks average)
	file_out.open(std::string(ROOT_PATH) +  "/Data/01.1a_mean_value_random_generator_M2.dat");
	blocks_averege(1e5, 1e2, rnd, file_out, characteristic_function_simple);
	file_out.close();


    //-------------------- PART 2 --------------------

    //Method 2 - Pointer function mode
	file_out.open(std::string(ROOT_PATH) +  "/Data/01.1b_mean_value_standard_deviation_M2.dat");
	blocks_averege(1e5, 1e2, rnd, file_out, characteristic_function_sigma);
	file_out.close();


    //-------------------- PART 3 --------------------

    std::ofstream file_out_chi(std::string(ROOT_PATH) +  "/Data/01.1c_chi_squared.dat"); 
	double sum_chi = 0., rand_chi;
	int M_chi = 100;
	int n = 10e4;
	double step = static_cast<double>(M_chi) / static_cast<double>(n);
	double expected = 1./step;
	std::vector<int> n_i(M_chi);

	for(int i = 0; i < M_chi; ++i){
		sum_chi = 0.;
        std::fill(n_i.begin(), n_i.end(), 0);
		for(int j = 0; j < n; ++j){
			rand_chi = rnd.Rannyu();

            //using the floor function avoids nesting of for loops
			n_i[floor(static_cast<double>(M_chi) * rand_chi)]++;

            //nesting loops - explicit control (deprecated)
			/*for(int k = 0;  k < M_chi; ++k){
				if (rand_chi < step * 10 * static_cast<double>(k+1) && 
					rand_chi > step * 10 * static_cast<double>(k)){
                    std::cout << k <<std::endl;
					n_i[k]++;
					break;
				}
			}*/
		}
		for(int j = 0; j < M_chi; ++j){
			sum_chi += pow(n_i[j] - expected, 2) * step;
            
		}
		file_out_chi << i << std::setw(20) << sum_chi << std::endl;
	}
	file_out_chi.close();

    //Random Save Seed
    rnd.SaveSeed();
}