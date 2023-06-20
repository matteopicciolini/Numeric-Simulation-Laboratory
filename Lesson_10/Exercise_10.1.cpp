#include "Library_10.h"

int main (int argc, char* argv[]){
    std::string print_city;
    std::string migr;
    bool migrate;
    int migration_freq;
    Usage(argc, argv, print_city, migr, migration_freq);
    std::cout << migr << std::endl;
    if (migr == "true"){
        migrate = true;
    }else if (migr == "false"){
        migrate = false;
    }
    else{exit(-1);}

    Random rnd;
    Random_Start(rnd);
    Delete_old_files("10.1_");
    
    //------------------------Progress bar declaration--------------------------
	int stp2 = NGeneration / 10;
	int s = 0;
	//--------------------------------------------------------------------------

    Task task;
    Population population;
    
    task.load_cities("American_capitals.dat");

    
    population.set_configuration(rnd);
    task.eval(population);    
    task.sort_population(&population);

    //MPI
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;

    
    Random_Start(rnd, rank); //inizializzo ciascuno nodo su coppie di numeri primi differenti    

    for(int i = 0; i < NGeneration; ++i){
        if(!migrate || (migrate && i % migration_freq != 0) || i == 0){ // continenti indipendenti -> no migrazioni
            population.mutate(rnd);
        
        }
        else if(migrate && i % migration_freq == 0){
            int h = 0, giver = 0, receiver = 1;
            int migrator[Ngenes];
            int itag = 1;

            if(rank == 0){
                h = (int)(population.get_n_individuals() * (1 - pow(rnd.Rannyu(), 20))) - 1; // scelgo un individuo come al solito
                do{
                    giver = (int)rnd.Rannyu(0, size);             // scelgo il core di partenza
                    receiver = (int)rnd.Rannyu(0, size);          // scelgo il core di arrivo... 
                }while(receiver == giver); // finchè giver e receiver sono uguali continua il ciclo
                std::cout << std::endl << "Migration (h, giver, receiver): " << h << ", " << giver << ", " << receiver << std::endl;
            }

            MPI_Bcast(&h, 1, MPI_INTEGER, 0, MPI_COMM_WORLD); // il core 0 manda a tutti gli altri
            MPI_Bcast(&giver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD); 
            MPI_Bcast(&receiver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

            if(rank == giver){
                for(int i = 0; i< Ngenes; i++){
                    migrator[i] = population.chromosomes[h].get_gene(i);
                }
                MPI_Send(migrator, Ngenes, MPI_INTEGER, receiver, itag, MPI_COMM_WORLD);
            }

            // se sono il ricevitore lo salvo
            if(rank == receiver){
                MPI_Recv(migrator, Ngenes, MPI_INTEGER, giver, itag, MPI_COMM_WORLD, &stat);
                population.chromosomes[h].set_gen(migrator); 
                population.chromosomes[h].check();
            }
        }

        task.eval(population);    
        task.sort_population(&population);

        // salvo i risultati
        if(print_city == "true"){
            task.print_cities(i, population.chromosomes[population.get_n_individuals() - 1], rank);
        }
        //prendo il best
        task.print_bests_len_ave(i, population.get_n_individuals() / 2, population, rank);
        task.print_best_len(i, population, rank);

        //----------------------------Progress Bar------------------------------
        //if (i % stp2 == 0) Progress_bar(s, i);
        if(i%10 == 0) std::cout << "Generation n: " << i  << " of " << NGeneration << std::endl;
		//----------------------------------------------------------------------
    }
    std::cout << Green("The process has terminated successfully. Progress: ▪▪▪▪▪▪▪▪▪▪ 100%") << std::endl;

    MPI_Finalize();
    return 0;
}