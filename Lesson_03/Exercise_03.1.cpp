#include "Library_03.h"


int main (int argc, char *argv[]){

    Random rnd;
	Random_Start(rnd);

    int M = 1E7; //n_steps
	int N = 100; //n block
	int L = M / N; //steps per block

    //open files
    std::ofstream file_output_directly_call(std::string(ROOT_PATH) + "/Data/03.1a_directly_call_option.dat");
    std::ofstream file_output_directly_put(std::string(ROOT_PATH) + "/Data/03.1a_directly_put_option.dat");
    std::ofstream file_output_discretized_call(std::string(ROOT_PATH) + "/Data/03.1b_discretized_call_option.dat");
    std::ofstream file_output_discretized_put(std::string(ROOT_PATH) + "/Data/03.1b_discretized_put_option.dat");

    //-------------------------------------------- direct sampling ------------------------------------------
    //initialisation
	double sum_call = 0., sum_put = 0.;

	double S_0 = 100.;
	double T = 1.;
    double K = 100.;
    double r = 0.1;
	double sigma = 0.25;
	double S;

    double ave_i_call = 0., ave2_i_call = 0.;
    double ave_i_put = 0., ave2_i_put = 0.;

    double acc_call = 0., acc_2_call = 0.;
    double acc_put = 0., acc_2_put = 0.;

	double sum_prog_call, su2_prog_call, err_prog_call;
    double sum_prog_put, su2_prog_put, err_prog_put;



	for (int i = 0; i < N; ++i){ //iteration on blocks
		sum_call = 0.;
		sum_put = 0.;
		for (int j = 0; j < L; ++j){ //step in block
			double rand = rnd.Gauss(0, T);
            S = S_0 * exp((r - 0.5 * pow(sigma, 2)) * T + sigma * rand);
            sum_call += exp(-r * T) * static_cast<double>(std::max(0., S - K));
            sum_put += exp(-r * T) * static_cast<double>(std::max(0., K - S));
		}

        //data blocking call
        ave_i_call = sum_call/static_cast<double>(L);
        ave2_i_call = pow(ave_i_call, 2);

        acc_call += ave_i_call;
		acc_2_call += ave2_i_call;

        sum_prog_call = acc_call / static_cast<double>(i + 1);
		su2_prog_call = acc_2_call / static_cast<double>(i + 1);
		err_prog_call = error(sum_prog_call, su2_prog_call, i);
		file_output_directly_call << i + 1<< std::setw(20) << sum_prog_call << std::setw(20) << err_prog_call << std::endl;

        //data blocking put
        ave_i_put = sum_put/static_cast<double>(L);
        ave2_i_put = pow(ave_i_put, 2);
	
        acc_put += ave_i_put;
		acc_2_put += ave2_i_put;

        sum_prog_put = acc_put / static_cast<double>(i + 1);
		su2_prog_put = acc_2_put / static_cast<double>(i + 1);
		err_prog_put = error(sum_prog_put, su2_prog_put, i);
		file_output_directly_put << i + 1 << std::setw(20) << sum_prog_put << std::setw(20) << err_prog_put << std::endl;
	}
    file_output_directly_call.close();
    file_output_directly_put.close();


    //-------------------------------------------- indirect sampling ------------------------------------------
    //initialisation
    ave_i_call = 0.;
    ave2_i_call = 0.;
    ave_i_put = 0.;
    ave2_i_put = 0.;
    acc_call = 0.; 
    acc_2_call = 0.;
    acc_put = 0.;
    acc_2_put = 0.;
    
    double n = 0.;
    double t = T / n;
    for (int i = 0; i < N; ++i){ //iteration on blocks
		sum_call = 0.;
		sum_put = 0.;
		for (int j = 0; j < L; ++j){ //step in block
			double rand = rnd.Gauss(0, T);
            S = S_0 * exp((r - 0.5 * pow(sigma, 2)) * T + sigma * rand);
            for(int z = 0; z < n; ++z){ // indirect sampling
       			S = S * exp((r - 0.5 * pow(sigma, 2.)) * (t * (z + 1) - t * z) + sigma * rnd.Gauss(0., 1.) * sqrt(t * (z + 1) - t * z));
      		}
			sum_call += exp(-r * T) * static_cast<double>(std::max(0., S - K));
			sum_put += exp(-r * T) * static_cast<double>(std::max(0., K - S));
		}

        //data blocking call
        ave_i_call = sum_call/static_cast<double>(L);
        ave2_i_call = pow(ave_i_call, 2);

        acc_call += ave_i_call;
		acc_2_call += ave2_i_call;

        sum_prog_call = acc_call / static_cast<double>(i + 1);
		su2_prog_call = acc_2_call / static_cast<double>(i + 1);
		err_prog_call = error(sum_prog_call, su2_prog_call, i);
		file_output_discretized_call << i + 1<< std::setw(20) << sum_prog_call << std::setw(20) << err_prog_call << std::endl;


        //data blocking put
        ave_i_put = sum_put/static_cast<double>(L);
        ave2_i_put = pow(ave_i_put, 2);
	
        acc_put += ave_i_put;
		acc_2_put += ave2_i_put;

        sum_prog_put = acc_put / static_cast<double>(i + 1);
		su2_prog_put = acc_2_put / static_cast<double>(i + 1);
		err_prog_put = error(sum_prog_put, su2_prog_put, i);
		file_output_discretized_put << i + 1 << std::setw(20) << sum_prog_put << std::setw(20) << err_prog_put << std::endl;
	}

    file_output_discretized_call.close();
    file_output_discretized_put.close();


    //save seed
    rnd.SaveSeed();
    
    return 0;
}