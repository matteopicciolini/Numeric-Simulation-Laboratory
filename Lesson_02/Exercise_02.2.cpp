#include "Library_02.h"

int main (int argc, char *argv[]){

    //Random Start
    Random rnd;
	Random_Start(rnd);

    int M = 1e4;
	double step_lenght = 1.;
	int n_steps = 100;
	int n_blocks = 100;
    int walks_per_block = M / n_blocks;

    std::ofstream file_out_r_module_square_discrete(std::string(ROOT_PATH) + "/Data/02.2_discrete_random_walk_r_module_square.dat");
	std::vector<std::vector<double>> matrix_path_discrete(n_blocks, std::vector<double>(n_steps));

    // Blocks cycle
    for(int i = 0; i < n_blocks; ++i){
        // Walk in a block cycle
        for(int j = 0; j < walks_per_block; ++j){
            // Restart from (0, 0, 0)
            std::vector<double> pos_discr(3);
	        std::fill(pos_discr.begin(), pos_discr.end(), 0);

            // Step in a walk cycle
            for(int k = 0; k < n_steps; ++k){ 
                // 1 step
                random_walk_discrete_step(rnd, pos_discr, step_lenght);
                // for each walk, save (|r_k|^2)/walks_per_block
                matrix_path_discrete[k][i] += ( pow(pos_discr[0],2) + pow(pos_discr[1],2) + pow(pos_discr[2],2) ) / walks_per_block;
            }    
        }
    }
    //blocking average
    double sum, sum2, err_prog;
    for (int k = 0; k < n_steps; ++k){
        // fix k, calc average in a block
        sum = 0.;
        sum2 = 0.;
        for (int i = 0; i < n_blocks; ++i){
            sum += matrix_path_discrete[k][i]/n_blocks;
            sum2 += pow(matrix_path_discrete[k][i], 2)/n_blocks;
        }
        err_prog = error(sum, sum2, n_blocks)/(2 * sqrt(sum));
        file_out_r_module_square_discrete << (k + 1)<< std::setw(20) << sqrt(sum) << std::setw(20) << err_prog << std::endl;
    }




    std::ofstream file_out_r_module_square_continue(std::string(ROOT_PATH) + "/Data/02.2_continue_random_walk_r_module_square.dat");
	std::vector<std::vector<double>> matrix_path_continue(n_blocks, std::vector<double>(n_steps));

    // Blocks cycle
    for(int i = 0; i < n_blocks; ++i){
        // Walk in a block cycle
        for(int j = 0; j < walks_per_block; ++j){
            // Restart from (0, 0, 0)
            std::vector<double> pos_cont(3);
	        std::fill(pos_cont.begin(), pos_cont.end(), 0);

            // Step in a walk cycle
            for(int k = 0; k < n_steps; ++k){ 
                // 1 step
                random_walk_continue_step(rnd, pos_cont, step_lenght);
                // for each walk, save (|r_k|^2)/walks_per_block
                matrix_path_continue[k][i] += ( pow(pos_cont[0],2) + pow(pos_cont[1],2) + pow(pos_cont[2],2) ) / walks_per_block;
            }    
        }
    }

    //blocking average
    for (int k = 0; k < n_steps; ++k){
        // fix k, calc average in a block
        sum = 0.;
        sum2 = 0.;
        for (int i = 0; i < n_blocks; ++i){
            sum += matrix_path_continue[k][i]/n_blocks;
            sum2 += pow(matrix_path_continue[k][i], 2)/n_blocks;
        }
        err_prog = error(sum, sum2, n_blocks)/(2 * sqrt(sum));
        file_out_r_module_square_continue << (k + 1)<< std::setw(20) << sqrt(sum) << std::setw(20) << err_prog << std::endl;
    }




    //Let's simulate a 3D Random Walk

	std::vector<double> cartesian_coordinates_continue(3);
	std::vector<double> cartesian_coordinates_discrete(3);
    std::ofstream file_out_discr(std::string(ROOT_PATH) + "/Data/02.2_discrete_random_walk_path.dat");
	std::ofstream file_out_cont(std::string(ROOT_PATH) + "/Data/02.2_continue_random_walk_path.dat");

	for (int i = 0; i < M; ++i){
        
		random_walk_discrete_step(rnd, cartesian_coordinates_discrete, step_lenght);
		file_out_discr << cartesian_coordinates_discrete[0] << std::setw(20) 
                       << cartesian_coordinates_discrete[1] << std::setw(20) 
                       << cartesian_coordinates_discrete[2] << std::endl;

		random_walk_continue_step(rnd, cartesian_coordinates_continue, step_lenght);
		file_out_cont << cartesian_coordinates_continue[0] << std::setw(20) 
                      << cartesian_coordinates_continue[1] << std::setw(20) 
                      << cartesian_coordinates_continue[2] << std::endl;
	}
    file_out_discr.close();
	file_out_cont.close();


    //Extractions on the unit sphere

    std::ofstream file_out_distr_unif_sphere(std::string(ROOT_PATH) + "/Data/02.2_distr_unif_sphere.dat");
    std::ofstream file_out_distr_non_unif_sphere(std::string(ROOT_PATH) + "/Data/02.2_distr_non_unif_sphere.dat");
    std::vector<double> pos_cart_unif(3);
    std::vector<double> pos_cart_non_unif(3);

    for(int i = 0; i < 5e3; ++i){
        std::fill(pos_cart_unif.begin(), pos_cart_unif.end(), 0);
        random_walk_continue_step(rnd, pos_cart_unif, step_lenght);
        file_out_distr_unif_sphere << (i + 1)<< std::setw(20) 
                                   << pos_cart_unif[0] << std::setw(20) 
                                   << pos_cart_unif[1] << std::setw(20) 
                                   << pos_cart_unif[2] << std::endl;
        
        double theta = rnd.Rannyu(0, 2.* M_PI);
        double phi = rnd.Rannyu(0, M_PI);
        pos_cart_non_unif[0] = step_lenght * sin(phi) * cos(theta);
        pos_cart_non_unif[1] = step_lenght * sin(phi) * sin(theta);
        pos_cart_non_unif[2] = step_lenght * cos(phi);
        file_out_distr_non_unif_sphere << (i + 1)<< std::setw(20) 
                                       << pos_cart_non_unif[0] << std::setw(20) 
                                       << pos_cart_non_unif[1] << std::setw(20) 
                                       << pos_cart_non_unif[2] << std::endl;  
    }
        
    file_out_distr_unif_sphere.close();
    file_out_distr_non_unif_sphere.close();

    rnd.SaveSeed();
}






