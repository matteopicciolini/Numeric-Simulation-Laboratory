#include "Library_01.h"

int main (int argc, char *argv[]){

    //Random Start
	Random rnd;
	Random_Start(rnd);

    //Pointer function method
    std::ofstream file_out("../../Data/01.3_buffon.dat");
    blocks_averege(1e5, 100, rnd, file_out, characteristic_function_buffon, characteristic_function_buffon_2);

    //explicit method
    /*
	std::vector<double> ave;
	int nhit = 0;
	int M = 1e5;
	double x, y;
	double ave_i = 0., ave2_i = 0.;
	double acc = 0., acc_2 = 0.;
	double sum_prog, su2_prog, err_prog;
    int N = M/100;
	for(int i = 0; i < 100; ++i){
        nhit = 0.;
		for (int j = 0; j < N; ++j){
			double y_0 = rnd.Rannyu(0, 1);
            do{
                x = rnd.Rannyu();
                y = rnd.Rannyu(-1.,1.);
                
            }while((pow(x,2) + pow(y,2)) > 1.);
            
            double cos_theta = y / sqrt((pow(x,2) + pow(y,2)));
            
            if(std::floor(y_0 + (L/2. * cos_theta)) != std::floor(y_0 - (L/2. * cos_theta))){
                nhit++;
            }
		}
        ave_i = 2. * L * static_cast<double>(N) / (static_cast<double>(nhit) * 1.);
		ave2_i = pow(ave_i, 2);
		
		acc += ave_i;
		acc_2 += ave2_i;
		
		sum_prog = acc / static_cast<double>(i + 1);
		su2_prog = acc_2 / static_cast<double>(i + 1);
		err_prog = error(sum_prog, su2_prog, i);
		file_out << (i+1)*N << std::setw(20) << sum_prog << std::setw(20) << err_prog << std::endl;
	}
*/

    file_out.close();


    //Random Save Seed
    rnd.SaveSeed();

    return 0;
}