#include "Library_10.h"

int main (int argc, char* argv[]){
    std::string print_city;
    std::string migr;
    bool migrate;
    int migration_freq;
    Usage(argc, argv, print_city, migr, migration_freq);
    if (migr == "true"){
        migrate = true;
    }else if (migr == "false"){
        migrate = false;
    }
    else{exit(-1);}

    //MPI
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;


    Random rnd;
    Random_Start(rnd);
    Delete_old_files("10.1_" + migr + "_" + std::to_string(rank) + "_best_len");
    
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

    
    Random_Start(rnd, rank); //inizializzo ciascuno nodo su coppie di numeri primi differenti    

    for(int i = 0; i < NGeneration; ++i){
        if(!migrate || (migrate && i % migration_freq != 0) || i == 0){ // continenti indipendenti -> no migrazioni
            population.reproduce(rnd);
        
        }
        else if(migrate && i % migration_freq == 0){
            int h = 0, giver = 0, receiver = 1;
            int migrator[Ngenes];
            int itag = 1;

            if(rank == 0){
                h = (int)(population.get_n_individuals() * (1 - pow(rnd.Rannyu(), 6))) - 1; // scelgo un individuo come al solito
                do{
                    giver = (int)rnd.Rannyu(0, size);             // scelgo il core di partenza
                    receiver = (int)rnd.Rannyu(0, size);          // scelgo il core di arrivo... 
                }while(receiver == giver); // finchè giver e receiver sono uguali continua il ciclo
            }

            MPI_Bcast(&h, 1, MPI_INTEGER, 0, MPI_COMM_WORLD); // il core 0 manda a tutti gli altri
            MPI_Bcast(&giver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD); 
            MPI_Bcast(&receiver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

            if(rank == giver){
                for(int i = 0; i< Ngenes; i++){
                    migrator[i] = population.individuals[h].get_gene(i);
                }
                MPI_Send(migrator, Ngenes, MPI_INTEGER, receiver, itag, MPI_COMM_WORLD);
            }

            // se sono il ricevitore lo salvo
            if(rank == receiver){
                MPI_Recv(migrator, Ngenes, MPI_INTEGER, giver, itag, MPI_COMM_WORLD, &stat);
                population.individuals[h].set_gene(migrator); 
                population.individuals[h].check();
            }
        }

        task.eval(population);    
        task.sort_population(&population);

        // salvo i risultati
        if(print_city == "true"){
            task.print_cities(i, population.individuals[population.get_n_individuals() - 1], rank, migr);
        }
        //prendo il best
        task.print_bests_len_ave(i, population.get_n_individuals() / 2, population, rank, migr);
        task.print_best_len(i, population, rank, migr);

        //----------------------------Progress Bar------------------------------
        if (i % stp2 == 0) Progress_bar(s, i, rank);
		//----------------------------------------------------------------------
    }
    std::cout << Green("The process has terminated successfully for rank " + std::to_string(rank) + ". Progress: ▪▪▪▪▪▪▪▪▪▪ 100%") << std::endl;

    MPI_Finalize();
    return 0;
}