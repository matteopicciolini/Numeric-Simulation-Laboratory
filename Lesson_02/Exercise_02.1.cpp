#include "Library_02.h"

int main (int argc, char *argv[]){

    //Random Start
    Random rnd;
	Random_Start(rnd);

    //std::ofstream file_out("../../Data/02_integral_uniform_sampling.dat");
	
	
    //double xmax = 1.;
	//double xmin = 0.;
    /*
    int point = 1e7;
    double x;
    double sum = 0;
    for(int i = 0; i < point; ++i){
        x = rnd.Rannyu();
        sum += 0.5 * M_PI * cos(M_PI * 0.5 * x);
    }
    std::cout << (xmax - xmin) * sum / static_cast<double>(point) << std::endl;*/

    double M = 1e7;
	double N = 100;
    std::ofstream file_out(std::string(ROOT_PATH) + "/Data/02.1_integral_uniform_sampling.dat");
    std::ofstream file_out_imp_smp(std::string(ROOT_PATH) + "/Data/02.1_integral_importance_sampling.dat"); 
    blocks_averege(M, N, rnd, file_out, characteristic_function_uniform_1, characteristic_function_uniform_2);
    blocks_averege(M, N, rnd, file_out_imp_smp, characteristic_function_sampling_1, characteristic_function_uniform_2);
    file_out.close();
	file_out_imp_smp.close();

    //Random Save Seed
    rnd.SaveSeed();
}