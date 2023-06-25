#include "Library_09.h"

void Random_Start(Random &random_generator){
	int seed[4];
	int p1, p2;
	std::ifstream Primes(std::string(ROOT_PATH) + "/random-library/Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else {std::cerr << "PROBLEM: Unable to open Primes" << std::endl;}
	
	Primes.close();
	
	std::ifstream input(std::string(ROOT_PATH) + "/random-library/seed.in");
	std::string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				random_generator.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else {std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;}
}

void swap(Individual& a, Individual& b) {
    Individual temp = a;
    a = b;
    b = temp;
}

int partition(Population &population, int low, int high) {
    double pivot = population.individuals[high].fitness;
    int i = (low - 1);

    for (int j = low; j < high; j++) {
        if (population.individuals[j].fitness <= pivot) {
            i++;
            swap(population.individuals[i], population.individuals[j]);
        }
    }
    swap(population.individuals[i + 1], population.individuals[high]);

    return (i + 1);
}

void quickSort(Population &population, int low, int high) {
  if (low < high) {
    int pi = partition(population, low, high);

    quickSort(population, low, pi - 1); 
    quickSort(population, pi + 1, high); 
   }
}

std::string Green(std::string green){
	return std::string("\033[1;32m") + green + "\033[0m";
}

std::string Gray(std::string gray){
	return std::string("\033[1;90m") + gray + "\033[0m";
}

std::string Red(std::string red){
	return std::string("\033[1;31m") + red + "\033[0m";
}

void Progress_bar(int& s, int& prog){
    std::string perc = "▪▪▪▪▪▪▪▪▪▪";
    std::cout << Red("Genetic Algoritm is running. Progress: ");
	std::cout << Green(perc.substr(0, 3 * s++));
	std::cout << Gray(perc.substr(3 * s, 3 * 10));
	std::cout << " " << int(prog * 100.0 / NGeneration) << " %\r";
	std::cout.flush();
}

void Delete_old_files(std::string pattern){
    std::filesystem::path directory_path = std::string(ROOT_PATH) + "/Data";
	std::string confirm = "";
	std::cout << "This program will delete files in Data with pattern '" + pattern + "'. Press <enter> to confirm." << std::endl;
	std::cin.ignore();
	for (auto& file : std::filesystem::directory_iterator(directory_path)) {
			if (std::filesystem::is_regular_file(file) && file.path().filename().string().find(pattern) != std::string::npos) {
					std::filesystem::remove(file.path());
					std::cout << "Deleted file: " << file.path() << std::endl;
			}
	}
  std::cout << std::endl;
}

void Usage(int argc, char* argv[], std::string &circ, std::string &print_city){
    if(argc == 3){
        circ = static_cast<std::string>(argv[1]);
        print_city = static_cast<std::string>(argv[2]);
        if (circ != "true" && circ != "false" && print_city != "true" && print_city != "false"){ 
            std::cerr << "Usage: ./Exercise_09.1 <circ> <print_city> with <circ> and <print_city> = {true, false}" << std::endl;
            exit(-1);
        }
    }
    else{
        std::cerr << "Usage: ./Exercise_09.1 <circ> <print_city> with <circ> and <print_city> = {true, false}" << std::endl;
        exit(-1);
    }
}

std::string intToStringWithLeadingZeros(int value, int width) {
    std::stringstream ss;
    ss << std::setw(width) << std::setfill('0') << value;
    return ss.str();
}


// ----------------------------- CLASSES ---------------------------------

Individual::Individual(){}
Individual::~Individual(){}

Individual& Individual::operator= (const Individual& individual){
    for(int i = 0; i < n_genes; ++i) genes[i] = individual.genes[i];
    // fitness
    fitness = individual.fitness;
    //lughezza
    len = individual.len;
    return *this;
}

void Individual::set_gene(int vec[n_genes]){for(int i = 0; i < n_genes; ++i) genes[i] = vec[i];}
void Individual::set_fitness(double fitness){this->fitness = fitness;}
void Individual::set_len(double len){this->len = len;}
double Individual::get_fitness(){return fitness;}
double Individual::get_len(){return len;}
int Individual::get_n_genes(){return n_genes;}
int Individual::get_gene(int i){return genes[i];}

void Individual::print(){
   for(int i = 0; i < n_genes; ++i) std::cout << genes[i] << " ";
   std::cout << std::endl;
}


void Individual::fill(Random &rnd){
    if(genes[0] != 0) {
            std::cerr << Red("Alert:") << " individual already full." << std::endl; 
            this->empty();
        }
    
    int r = (int)(rnd.Rannyu() * n_genes);
    genes[0] = 1;
    // Take integers in order
    for(int i = 2; i < n_genes + 1; ++i){
        
        // finchè non trovo una posizione vuota in gen, genero un nuovo numero random
        while(genes[r]!=0){
            r = (int)(rnd.Rannyu() * n_genes);
        }
        genes[r] = i; //riempio la posizione
    }
}

void Individual::empty(){
    for(int i = 0; i < n_genes; ++i){
        genes[i] = 0; 
    }
}


void Individual::swap(int position_1, int position_2){
   if(position_1 == 0) {
      position_1 = 1;
   }
   if(position_2 == 0) {
      position_2 = 1;
   }

   int g1 = genes[pbc(position_1)];
   int g2 = genes[pbc(position_2)];
   
   genes[pbc(position_1)] = g2;
   genes[pbc(position_2)] = g1;
}



void Individual::shift(int position, int m, int n){

    if(position == 0) {
        position = 1;
    }
    if(position > (n_genes - 1)) {
        position = pbc(position);
    }

    if(m >= n_genes - 1) {
        m = 1;
    }

    int block[m];

    // ESEMPIO: [1 4 3 2 5]
    // salvo il blocco da shiftare
    // [4, 3]
    for(int i = 0; i < m; ++i){
        block[i] = genes[pbc(position + i)];
    }

    // [1 2 5 2 5]
    // shifto n volte i geni complementari per riempire il vettore...
    for(int i = 0; i < n; ++i){
        genes[pbc(position + i)] = genes[pbc(position + m + i)];
    }

    // [1 2 5 4 3]
    for(int i = 0; i < m; ++i){
        genes[pbc(position + n + i)] = block[i];
    }
}

void Individual::permutate(int position_1, int position_2, int m){
    if(position_1 == 0) {
        position_1 = 1;
    }
    if(position_2 == 0) {
        position_2 = 1;
    }
    if(position_1 > (n_genes - 1)) {
        position_1 = pbc(position_1);
    }
    if(position_2 > (n_genes - 1)) {
        position_2 = pbc(position_2);
    }
    if(m > n_genes/2) {
        m = 1;
    }
    if(m > abs(position_2 - position_1)){
        m = abs(position_1-position_2);
    }
    if(position_1 > position_2){
        int appo = position_2;
        position_2 = position_1;
        position_1 = appo;
    }

    int block[m];

    for(int i = 0; i < m; ++i){
        block[i] = genes[pbc(position_1 + i)];
    }

    for(int i = 0; i < m; ++i){
        genes[pbc(position_1+i)] = genes[pbc(position_2+i)];
    }

    for(int i = 0; i < m; ++i){
        genes[pbc(position_2 + i)] = block[i];
    }
}

void Individual::reverse(int position, int m){
    if(position == 0) {
        position = 1;
    }
    if(position > (n_genes - 1)) {
        position = pbc(position);
    }
    if(m > n_genes - 1) {
        m = 2;
    }

    int block[m];

    // salvo il blocco da invertire
    for(int i = 0; i < m; ++i){
        block[i] = genes[pbc(position + i)];
    }

    // riempio al contrario lo spazio vuoto lasciato dal blocco
    for(int i = 0; i < m; ++i){
        genes[pbc(position + i)] = block[(m - 1) - i];
    }
}

void Individual::crossover(int position, int len, Individual individual_2){
    int block[len];
    if(position == 0) {
        position = 1;
    }
    // salvo il blocco principale
    for(int i = 0; i < len; ++i){
        block[i] = genes[pbc(position + i)];
    } 

    // completa i percorsi con le città mancanti aggiungendole nell'ordine in cui compaiono nella consort 
    int count = 0;
    
    for(int j = 1; j < n_genes; ++j){
        // ciclo sui geni nei blocchi 
        for(int i = 0; i < len; ++i){
            // quando trovo nell'altro genitore il gene presente
            // nel blocco, lo salvo nel buco scavato all'inizio
            if(individual_2.genes[pbc(j)] == block[i] ){
                genes[pbc(position+count)] = individual_2.genes[pbc(j)];
                count++;
            }
        }
    } 
}

void Individual::check(){
    if(this->genes[0] != 1){
        std::cerr << Red("Alert:") << " Individual::check has failed." << std::endl;
        exit(-1);
    }

    int sum = 0;
    for(int i = 0; i < this->n_genes; ++i){
        sum += this->genes[i];
    }

    if(sum != (this->n_genes * (this->n_genes + 1)) / 2){
        std::cerr << Red("Alert:") << " Individual::check has failed." << std::endl;
        exit(-1);
    }
}

int Individual::pbc(int pos){
    //se sfora e non è multiplo della lunghezza restituisco il modulo
    if(pos > (n_genes - 1) && pos % (n_genes - 1) !=0)    pos = pos % (n_genes - 1);
    //se sfora ed è multiplo della lunghezza restituisco la lunghezza
    if(pos > (n_genes - 1) && pos % (n_genes - 1) == 0)    pos = n_genes - 1;
    if(pos < 0)       pos = n_genes + pos % n_genes;
    return pos;
}

bool Individual::equals(Individual& individual_2) {
    if (individual_2.get_n_genes() != this->n_genes) {
        return false;
    }
    for (int i = 0; i < this->n_genes; ++i) {
        if (individual_2.get_gene(i) != this->get_gene(i)) {
            return false;
        }
    }
    return true;
}

Population::Population(Random &rnd){
    Individual individual;
    for(int i = 0; i < n_individuals; ++i){
        individual.empty();
        individual.fill(rnd);
        individuals[i] = individual;
        individuals[i].check();
    }
}
Population::~Population(){}

int Population::get_n_individuals(){
    return n_individuals;
}


void Population::reproduce(Random &rnd){
    
    Individual individual[n_individuals];
    int n_g = individual[0].get_n_genes();
    double expo = expon;
    int h = 0; 
    double len_av = 0;
    
    // Select the best individuals and calculate best length and average length
    for (int i = 0; i < n_individuals; ++i) {
        
        // (prima salvo lunghezze migliori)
        if(i > n_individuals / 2) len_av += this->individuals[i].len;          // Accumulate the length for calculating average
        if(i == n_individuals - 1) this->best_len = this->individuals[i].len;  // Assign the length of the best individual

        // Select an individual according to the power law distribution and save its genetic information
        h = (int)(n_individuals * (1 - pow(rnd.Rannyu(), expo))) - 1; // Select a high index based on power law
        //lim_sup
        if(h >= n_individuals) h = n_individuals - 1;
        //lim inf
        if(h < 0) h = 0;
        individual[i] = this->individuals[h];
    }
    this->best_len_ave = len_av / (double)(n_individuals / 2); 

    // mutazioni
    int pos, m, shift, len;
    int pos1, pos2, index;
    for(int i = 0; i < n_individuals; ++i){

        double prob = rnd.Rannyu();
        if(prob < mutation_probability[0]) {
            pos1 = rnd.Rannyu(1, n_g);
            pos2 = rnd.Rannyu(1, n_g);
            individual[i].swap(pos1, pos2);
        }

        prob = rnd.Rannyu();
        if(prob < mutation_probability[1]){
            pos = rnd.Rannyu(1, n_g);
            m = rnd.Rannyu(1, n_g);
            individual[i].reverse(pos, m);
        }

        prob = rnd.Rannyu();
        if(prob < mutation_probability[2]){
            pos = rnd.Rannyu(1, n_g);
            m = rnd.Rannyu(1, n_g-1);
            shift = rnd.Rannyu(1,n_g);
            individual[i].shift(pos, m, shift);
        }

        prob = rnd.Rannyu();
        if(prob < mutation_probability[3]){
            pos1 = rnd.Rannyu(1, n_g);
            pos2 = rnd.Rannyu(1, n_g);
            m = rnd.Rannyu(1, n_g * 0.6);
            individual[i].permutate(pos1, pos2, m);
        }

        prob = rnd.Rannyu();
        if(prob < mutation_probability[4]){
            pos = rnd.Rannyu(1, n_g);
            len = rnd.Rannyu(1, n_g);
            index = rnd.Rannyu(1, n_individuals);
            individual[i].crossover(pos, len, individual[index]);
        }

        //controllo e salvo per la generazione successiva
        individual[i].check();
        this->individuals[i] = individual[i];
    }
}


Task::Task(std::string circ){
    if(circ == "true"){
        this->task = "circ";
    }
    else if (circ == "false"){
        this->task = "square";
    }
    
}

Task::~Task(){}

// genera coordinate di città in su una circonferenza unitaria
void Task::generate_circular_cities(Random rnd, double radius){
    double phi;

    for(int i = 0; i < this->n_cities; ++i){
        phi = rnd.Rannyu() * 2 * M_PI;
        this->cities_x[i] = radius * cos(phi);
        this->cities_y[i] = radius * sin(phi);
    }
}

/// genera coordinate di città in un quadrato di lato 2
void Task::generate_squared_cities(Random rnd){
    for(int i = 0; i < this->n_cities; ++i){
        this->cities_x[i] = rnd.Rannyu(-1, 1);
        this->cities_y[i] = rnd.Rannyu(-1, 1);
    }
}

void Task::generate_cities(Random rnd){
    if(task == "circ"){
        generate_circular_cities(rnd);
    }
    else if(task == "square"){
        generate_squared_cities(rnd);
    }
    else{
        std::cerr << "Parameter <circ> must be 'true or 'false'. Exit." << std::endl;
        exit(-1);
    }
}


// stampa le coordinate delle città nell'ordine indicato da un individuo
void Task::print_cities(int generation, Individual chr){

    std::ofstream stream;
    stream.open(std::string(ROOT_PATH) + "/Data/09.1_" + task + "_city_coord_" + intToStringWithLeadingZeros(generation + 1, 3) + ".dat");

    int seq;
    for(int i = 0; i < n_cities; ++i){
        seq = chr.get_gene(i) - 1; // nell'ordine indicato dal individuo
        stream << cities_x[seq]  << " " << cities_y[seq] << std::endl;
    }
    stream.close();
}

// stampa la lunghezza media della migliore metà di individui in una popolazione
void Task::print_bests_len_ave(int generation, int part, Population populatio){
    std::ofstream str;
    str.open(std::string(ROOT_PATH) + "/Data/09.1_" + task + "_best_len_average.dat", std::ios_base::app);
    str << generation + 1 << " " << populatio.best_len_ave << std::endl;
    str.close();
}

// Stampa la lunghezza del miglior individuo di una popolazione
void Task::print_best_len(int generation, Population populatio){
    std::ofstream str;
    str.open(std::string(ROOT_PATH) + "/Data/09.1_" + task + "_best_len.dat", std::ios_base::app);
    str << generation  + 1 << " " << populatio.best_len << std::endl;
    str.close();
}

// prende un individuo e con l'ordine dei geni calcola la distanza tra le città e la fitness
void Task::eval_fitness(Individual &individual){
    double fitness = 0;
    double accu = 0;
    double acculen = 0;
    double len2 = 0; // lunghezza senza sqrt (per fitness)
    int ng = individual.get_n_genes();

    // somma dalla prima all'ultima
    for(int i = 0; i < ng - 1; ++i){

        len2 = pow( cities_x[individual.get_gene(i + 1) - 1] - cities_x[individual.get_gene(i) - 1], 2) 
            + pow( cities_y[individual.get_gene(i + 1) - 1] - cities_y[individual.get_gene(i) - 1], 2);
        accu += len2;
        acculen += sqrt(len2);
    }
    // distanza tra la prima e l'ultima
    len2 = pow(cities_x[individual.get_gene(0) - 1] - cities_x[individual.get_gene(ng - 1) - 1], 2) 
            + pow(cities_y[individual.get_gene(0) - 1] - cities_y[individual.get_gene(ng - 1) - 1], 2);
    accu += len2;
    acculen += sqrt(len2);

    fitness = 1.0 / accu;
    individual.fitness = fitness;
    individual.len = acculen;
}

// prende una popolazione e applica eval_fitness su tutti i suoi individui
void Task::eval(Population &population){ 
    for(int i = 0; i < population.get_n_individuals(); ++i){
        eval_fitness(population.individuals[i]);
    }
}

// prende una popolazione e mette in ordine i suoi individui dal peggiore al migliore, sulla base della fitness
void Task::sort_population(Population &population){
    quickSort(population, 0, population.get_n_individuals()-1);
}