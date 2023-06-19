#include "Library_09.h"
#include "Chromosome_test.h"



int main (int argc, char* argv[]){
    std::string circ;
    std::string print_city;
    Usage(argc, argv, circ, print_city);

    Random rnd;
    Random_Start(rnd);
    if (circ == "true"){
        Delete_old_files("09.1_circ_best");
    }
    else if(circ == "false"){
        Delete_old_files("09.1_square_best");
    }
    
    //------------------------Progress bar declaration--------------------------
	int stp2 = NGeneration / 10;
	int s = 0;
	//--------------------------------------------------------------------------

    Task task(circ);
    Population population;
    
    task.generate_cities(rnd);

    population.set_configuration(rnd);
    task.eval(population);    
    task.sort_population(&population);

    for(int i = 0; i < NGeneration; ++i){
        
        population.mutate(rnd);
        task.eval(population);    
        task.sort_population(&population);

        // salvo i risultati
        if(print_city == "true"){
            task.print_cities(i, population.chromosomes[population.get_n_individuals() - 1]);
        }
        //prendo il best
        task.print_bests_len_ave(i, population.get_n_individuals() / 2, population);
        task.print_best_len(i, population);

        //----------------------------Progress Bar------------------------------
        if (i % stp2 == 0) Progress_bar(s, i);
		//----------------------------------------------------------------------
    }
    std::cout << Green("The process has terminated successfully. Progress: ▪▪▪▪▪▪▪▪▪▪ 100%") << std::endl;
    return 0;
}