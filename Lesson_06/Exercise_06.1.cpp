#include "Library_06.h"

std::string random_lib_path = std::string(ROOT_PATH) + "/random-library/";
std::string input_path = "input-output/input.dat", metro_str, directory_path, pattern, pattern_2, eq, eq_str;

int main(int argc, char* argv[]){
    //Usage and phase choosing
    Usage(argc, argv);
    std::cout << std::endl;

    if(eq == "false"){
        for (int istep = 1; istep <= 1e5; ++istep){
            Move(metro);
        }
    }
    Input(); //Inizialization
    for(int iblk = 1; iblk <= nblk; ++iblk){
        Reset(iblk);   //Reset block averages
        for(int istep = 1; istep <= nstep; ++istep){
            Move(metro);
            Measure();
            Accumulate(); //Update block averages
        }
        Averages(iblk);   //Print results for current block
    }
    Results();
    ConfFinal(); //Write final configuration

    return 0;
}


void Usage(int argc, char* argv[]){
    if(argc == 2){
		eq = static_cast<std::string>(argv[1]);
		if(eq == "true"){
		    std::cout << "You choose EQUILIBRATION mode" << std::endl;
        }
        else if (eq == "false"){
            std::cout << "You choose SIMULATION mode" << std::endl;
        }
        else{
        std::cout << "You choose SIMULATION mode" << std::endl;
        }
	}
	else{
		std::cerr << "Usage: ./Exercise_06.1 <equilibration>' with <equilibration> = {true, false, null}" << std::endl;
        exit(-1);
	}
}

void Input(void){
    std::ifstream ReadInput;

    std::cout << "Classic 1D Ising model             " << std::endl;
    std::cout << "Monte Carlo simulation             " << std::endl << std::endl;
    std::cout << "Nearest neighbour interaction      " << std::endl << std::endl;
    std::cout << "Boltzmann weight exp(- Beta * H ), Beta = 1/T " << std::endl << std::endl;
    std::cout << "The program uses k_B=1 and mu_B=1 units " << std::endl;

    //Read seed for random numbers
    int p1, p2;
    std::ifstream Primes(random_lib_path + "Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    std::ifstream input(random_lib_path + "seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed, p1, p2);
    input.close();
        
    //Read input informations
    ReadInput.open(input_path);

    ReadInput >> temp;
    Beta = 1.0 / temp;
    std::cout << "Temperature = " << temp << std::endl;

    ReadInput >> nspin;
    std::cout << "Number of spins = " << nspin << std::endl;

    ReadInput >> J;
    std::cout << "Exchange interaction = " << J << std::endl;

    ReadInput >> h;
    std::cout << "External field = " << h << std::endl << std::endl;
        
    ReadInput >> metro; // if = 1 Metropolis else Gibbs

    ReadInput >> nblk;

    ReadInput >> nstep;

    if(metro == 1){
        metro_str = "metro";
    }
    else{
        metro_str = "gibbs";
    }
    if(eq == "true"){
        eq_str = "eq";
    }
    else{
        eq_str = "sim";
    }

    //delete old files
    Delete_old_files();


    if(metro == 1) std::cout << "The program perform Metropolis moves" << std::endl;
    else std::cout << "The program perform Gibbs moves" << std::endl;
    std::cout << "Number of blocks = " << nblk << std::endl;
    std::cout << "Number of steps in one block = " << nstep << std::endl << std::endl;
    ReadInput.close();


    //Prepare arrays for measurements
    iu = 0; //Energy
    ic = 1; //Heat capacity
    im = 2; //Magnetization
    ix = 3; //Magnetic susceptibility
    
    n_props = 4; //Number of observables

    //initial configuration
    for (int i = 0; i < nspin; ++i){
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
    }
        
    //Evaluate energy etc. of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    std::cout << "Initial energy = " << walker[iu]/(double)nspin << std::endl;
}

void Delete_old_files(){
    pattern_2 = "06.1_temps_" + eq_str + "_" + metro_str + "_" + std::to_string(h).substr(0, 4);
    pattern = "06.1_" + eq_str + "_" + metro_str + "_" + std::to_string(temp).substr(0, 3) + "_" + std::to_string(h).substr(0, 4);
    std::filesystem::path directory_path = std::string(ROOT_PATH) + "/Data";
    std::string confirm = "";
    std::cout << "This program will delete files in Data with pattern '" + pattern + "'. Press <enter> to confirm." << std::endl;
    std::cin.ignore();
    for (auto& file : std::filesystem::directory_iterator(directory_path)) {
        if (std::filesystem::is_regular_file(file) && file.path().filename().string().find(pattern) != std::string::npos){
            std::filesystem::remove(file.path());
            std::cout << "Deleted file: " << file.path() << std::endl;
        }	
        }
    std::cout << std::endl;   

    std::string input;
    std::cout << "This program will NOT delete files in Data with pattern '" + pattern_2 + "'. If you wish for these files to be deleted, press 1, otherwise, press Enter." << std::endl;
    getline(std::cin, input);

    if (input == "1") {
        for (auto& file : std::filesystem::directory_iterator(directory_path)) {
        if (std::filesystem::is_regular_file(file) && file.path().filename().string().find(pattern_2) != std::string::npos){
            std::filesystem::remove(file.path());
            std::cout << "Deleted file: " << file.path() << std::endl;
        }	
            }
    } else {
        std::cout << "Files not deleted." << std::endl;
    }
    std::cout << std::endl;   
}

void Move(int metro){
    int o;
    double p, energy_old, energy_new;
    double energy_up, energy_down;

    for(int i = 0; i < nspin; ++i){
        o = static_cast<int>(rnd.Rannyu() * nspin);
        if(metro == 1){
            
            int s_flip = -s[o];

            energy_old = Boltzmann(s[o], o);
            energy_new = Boltzmann(s_flip, o);
            
            double A = std::min(1., exp(-Beta * (energy_new - energy_old)));
            
            if(rnd.Rannyu() < A){
                s[o] = s_flip; 
                accepted++;
            }
            attempted++;
        }
        else{
            energy_up = Boltzmann(1, o);
            energy_down = Boltzmann(-1, o);
            p = exp(-Beta * (energy_down - energy_up));
            if (rnd.Rannyu() < 1 / (1 + p)) {            
                s[o] = 1;
                accepted++;
            }
            else{
                s[o] = -1;
                accepted++;
            }
            attempted++;
        }
    }
}

double Boltzmann(int sm, int ip){
    double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
    return ene;
}

void Measure(){
    //int bin;
    double u = 0.0, m = 0.0;

    //cycle over spins
    for (int i = 0; i < nspin; ++i){
        u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
        m += s[i];
    }
    walker[iu] = u;
    walker[ic] = u * u;
    walker[im] = m;
    walker[ix] = m * m;
}


void Reset(int iblk){ //Reset block averages

    if(iblk == 1){
        for(int i = 0; i < n_props; ++i){
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for(int i = 0; i < n_props; ++i){
        blk_av[i] = 0;
    }
    
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}


void Accumulate(void){ //Update block averages

    for(int i = 0; i < n_props; ++i){
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){ //Print results for current block
    const int wd = 20;
    
    std::cout << "Block number " << iblk << std::endl;
    std::cout << "Acceptance rate " << accepted / attempted << std::endl << std::endl;
    
    if(h == 0.){
        //INTERNAL ENERGY
        stima_u = blk_av[iu] / blk_norm / (double)nspin;
        glob_av[iu]  += stima_u;
        glob_av2[iu] += stima_u * stima_u;
        err_u = Error(glob_av[iu], glob_av2[iu], iblk);

        // CAPACITÀ TERMICA
        stima_c = Beta * Beta * (blk_av[ic] / blk_norm - pow(blk_av[iu] / blk_norm, 2)) / (double)nspin; 
        glob_av[ic]  += stima_c;
        glob_av2[ic] += stima_c * stima_c;
        err_c = Error(glob_av[ic], glob_av2[ic], iblk);
        
        // SUSCETTIVITÀ
        stima_x = Beta * blk_av[ix] / blk_norm / (double)nspin;
        glob_av[ix]  += stima_x;
        glob_av2[ix] += stima_x * stima_x;
        err_x = Error(glob_av[ix], glob_av2[ix], iblk);
        

        if(eq == "true"){
            std::ofstream Ene, Heat, Chi;
            Ene.open(std::string(ROOT_PATH) + "/Data/" + pattern + "_ene.dat", std::ios::app);
            Ene << iblk <<  std::setw(wd) << stima_u << std::setw(wd) << glob_av[iu] / (double)iblk << std::setw(wd) << err_u << std::endl;
            Ene.close();

            Heat.open(std::string(ROOT_PATH) + "/Data/" + pattern + "_heat.dat", std::ios::app);
            Heat << iblk <<  std::setw(wd) << stima_c << std::setw(wd) << glob_av[ic] / (double)iblk << std::setw(wd) << err_c << std::endl;
            Heat.close();

            Chi.open(std::string(ROOT_PATH) + "/Data/" + pattern + "_chi.dat", std::ios::app);
            Chi << iblk <<  std::setw(wd) << stima_x << std::setw(wd) << glob_av[ix] / (double)iblk << std::setw(wd) << err_x << std::endl;
            Chi.close();
        }
    }
    else{
        // MAGNETIZZAZIONE
        stima_m = blk_av[im] / blk_norm / (double)nspin; 
        glob_av[im]  += stima_m;
        glob_av2[im] += stima_m * stima_m;
        err_m = Error(glob_av[im], glob_av2[im], iblk);
        
        if(eq == "true"){
            std::ofstream Mag;
            Mag.open(std::string(ROOT_PATH) + "/Data/" + pattern + "_meg.dat", std::ios::app);
            Mag << iblk <<  std::setw(wd) << stima_m << std::setw(wd) << glob_av[im] / (double)iblk << std::setw(wd) << err_m << std::endl;
            Mag.close();
        }
    }
    std::cout << "----------------------------" << std::endl << std::endl;
}


void ConfFinal(void){

    std::ofstream WriteConf;
    std::cout << "Print final configuration to file config.final " << std::endl << std::endl;
    WriteConf.open("input-output/config.final");
    for (int i = 0; i < nspin; ++i){
        WriteConf << s[i] << std::endl;
    }
    WriteConf.close();
    rnd.SaveSeed("input-output/seed.out");
}

int Pbc(int i){  //Algorithm for periodic boundary conditions
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk){
    if(iblk == 1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void Results(){ // Fill files for plotting

    std::ofstream Ene, Heat, Mag, Chi;
    if(h == 0){
        // ENERGY
        Ene.open(std::string(ROOT_PATH) + "/Data/" + pattern_2 + "_ene.dat", std::ios::app);
        Ene << temp << " " << glob_av[iu] / (double)nblk << " " << err_u << std::endl;
        Ene.close();

        // HEAT CAPACITY
        Heat.open(std::string(ROOT_PATH) + "/Data/" + pattern_2 +"_heat.dat", std::ios::app);
        Heat << temp << " " << glob_av[ic] / (double)nblk << " " << err_c << std::endl;
        Heat.close();

        // MAGN SUSCEPTIBILITY 
        Chi.open(std::string(ROOT_PATH) + "/Data/" + pattern_2 + "_chi.dat", std::ios::app);
        Chi << temp << " " << glob_av[ix] / (double)nblk << " " << err_x << std::endl;
        Chi.close();
    }
    else{
        // MAGNETIZATION
        std::stringstream stream;
        Mag.open(std::string(ROOT_PATH) + "/Data/" + pattern_2 + "_mag.dat", std::ios::app);
        Mag << temp << " " << glob_av[im]/(double)nblk << " " << err_m << " " << h << std::endl;
        Mag.close();
    }
}