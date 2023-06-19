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

Population::Population(){}
Population::~Population(){}

void Population::set_configuration(Random &rnd){
   Chromosome chromosome;
   for(int i = 0; i < n_individuals; ++i){
      chromosome.empty();
      chromosome.fill(rnd);
      chromosomes[i] = chromosome;
      chromosomes[i].check();
   }
}

int Population::get_n_individuals(){
   return n_individuals;
}

/// Metodo complesso: travasa la popolazione in una generazione successiva,
/// generata per mutazioni dei migliori individui della generazione attuale.
void Population::mutate(Random &rnd){
    // variabili
    Chromosome chromosome[n_individuals];
    int ng = chromosome[0].get_n_genes();
    //const double *mp = chr_appo[0].mutation_probability;
    //int a,b,c;
    double expo = expon; // <<<<< questo parametro è scelto arbitrariamente
    int h = 0; 
    double len_av = 0;
    
    // metto la parte migliore della generazione attuale nel cromosoma di appoggio
    for (int i = 0; i < n_individuals; ++i) {
        
        // (prima salvo lunghezze migliori)
        if(i>n_individuals/2) len_av += this->chromosomes[i].len;          // media migliori
        if(i == n_individuals-1) this->best_len = this->chromosomes[i].len; // il migliore

        // estraggo uno dei migliori secondo legge di potenza, poi lo salvo
        h = (int)(n_individuals*(1-pow(rnd.Rannyu(),expo)))-1; // estrae un indice alto
        if(h>=n_individuals) h = n_individuals-1;                       // lim_sup
        if(h<0) h = 0;                                // lim_inf
        chromosome[i] = this->chromosomes[h];
        //chr_appo[i].Check();
    }
    this->best_len_ave = len_av/(double)(n_individuals/2); 

    // lo faccio evolvere e salvo per la successiva generazione
    int pos, m, shift, len;
        int pos1, pos2, index;
    for(int i = 0; i < n_individuals; ++i){

        double prob = rnd.Rannyu();
        if(prob < mutation_probability[0]) {
            //std::cout << "--PP--" << std::endl;
            pos1 = rnd.Rannyu(1,ng);
            pos2 = rnd.Rannyu(1,ng);
            chromosome[i].swap(pos1, pos2);
        }

        prob = rnd.Rannyu();
        if(prob < mutation_probability[1]){
            //std::cout << "--Inv--" << std::endl;
            pos = rnd.Rannyu(1,ng);
            m = rnd.Rannyu(1,ng);
            chromosome[i].reverse(pos, m);
        }

        prob = rnd.Rannyu();
        if(prob < mutation_probability[2]){
            //std::cout << "--Shift--" << std::endl;
            pos = rnd.Rannyu(1,ng);
            m = rnd.Rannyu(1,ng-1);
            shift = rnd.Rannyu(1,ng);
            chromosome[i].shift(pos, m, shift);
        }

        prob = rnd.Rannyu();
        if(prob < mutation_probability[3]){
            //std::cout << "--MP--" << std::endl;
            pos1 = rnd.Rannyu(1, ng);
            pos2 = rnd.Rannyu(1, ng);
            m = rnd.Rannyu(1, ng*0.6);
            chromosome[i].permutate(pos1, pos2, m);
            //std::cout << a <<" "<< b <<" "<< c << std::endl;
            //chr_appo[i].Print();
        }

        prob = rnd.Rannyu();
        if(prob < mutation_probability[4]){
            //std::cout << "--CR--" << std::endl;
            pos = rnd.Rannyu(1,ng); // pos
            len = rnd.Rannyu(1,ng); // len
            index = rnd.Rannyu(1,n_individuals); // altro individuo
            //chr_appo[i].Print();
            chromosome[i].crossover(pos, len, chromosome[index]);
            //chr_appo[i].Print();
        }

        // controllo e salvo per la generazione successiva
        chromosome[i].check();
        this->chromosomes[i] = chromosome[i];
    }
   //delete chr_appo;
}

//========================================================//
//       CHROMOSOME                                       //
//========================================================//

Chromosome::Chromosome(){} // costruttore
Chromosome::~Chromosome(){} // distruttore

Chromosome& Chromosome::operator= (const Chromosome& chromosome)
{
    // COPIES
    // vettore di geni
    for(int i = 0; i < n_genes; ++i) genes[i] = chromosome.genes[i];
    // fitness
    fitness = chromosome.fitness;
    //lughezza
    len = chromosome.len;

    // return the existing object so we can chain this operator
    return *this;
}

// methods

// get e set

void Chromosome::set_gen(int vec[n_genes]){
   for(int i = 0; i< n_genes; ++i) genes[i] = vec[i];
}

void Chromosome::set_fitness(double fitness){ 
   this->fitness = fitness;
}

void Chromosome::set_len(double len){ 
   this->len = len;
}

double Chromosome::get_fitness(){
   return fitness;
}

double Chromosome::get_len(){ 
   return len;
}

int Chromosome::get_n_genes(){
   return n_genes;
}

int Chromosome::get_gene(int i){
   return genes[i];
}


// altri metodi

/// stampa i geni del cromosoma
void Chromosome::print(){

   for(int i = 0; i < n_genes; ++i) std::cout << genes[i] << " ";
   std::cout << std::endl;
}

/// riempie il cromosoma con geni casuali
void Chromosome::fill(Random &rnd){
   //if(Gen[0] != 0) {std::cerr << "Can't fill full chromosome!" << std::endl; abort();}
   if(genes[0] != 0) {std::cerr << "Can't fill full chromosome: I clean it." << std::endl; this->empty();}
   
   int r = (int)(rnd.Rannyu() * n_genes);
   genes[0] = 1;
   // prendo interi in ordine
   for(int i = 2; i < n_genes + 1; ++i){
      
      while(genes[r]!=0){ // finchè la posizione è vuota ...
         r = (int)(rnd.Rannyu() * n_genes);
      }
      genes[r] = i; // ...la riempio...
   }
}

/// svuota il cromosoma
void Chromosome::empty(){
   for(int i = 0; i < n_genes; ++i){
      genes[i] = 0; 
   }
}


// accessori alle mutazioni

/// Applica pbc quando necessario durante le mutazioni
int Chromosome::pbc(int pos){
   if(pos > (n_genes - 1) && pos % (n_genes - 1) !=0)    pos = pos % (n_genes - 1); // non Ng, perchè devo saltare la posizione 0
   if(pos > (n_genes - 1) && pos % (n_genes - 1) == 0)    pos = n_genes - 1;
   if(pos < 0)       pos = n_genes + pos % n_genes;
   return pos;
}

/// Verifica indirettamente (può fallire ma è molto improbabile)
/// che il cromosoma contenga tutti i geni e il primo sia fisso.
/// Verifica che la somma totale dei geni sia corretta.
void Chromosome::check(){
if(this->genes[0] != 1){ // controllo che il primo gene sia sempre 1
        std::cerr << "Check: something went wrong. First gene is no more 1." << std::endl;
        exit(-1);
    }

    int sum = 0;
    for(int i = 0; i < this->n_genes; ++i){
        sum += this->genes[i];
    }

    if(sum != (this->n_genes * (this->n_genes + 1)) / 2){
        std::cerr << "Chromosome::check: something went wrong. Genes were not conserved through mutations." << std::endl;
        exit(-1);
    }
}


// Mutazioni

/// 1) permuta due geni
void Chromosome::swap(int position_1, int position_2){
   if(position_1 == 0) {
      // std::cerr << "Chromosome::PairPermut : can't move first gene: I move the second one." << std::endl;
      position_1 = 1;
   }
   if(position_2 == 0) {
      // std::cerr << "Chromosome::PairPermut : can't move first gene: I move the second one." << std::endl;
      position_2 = 1;
   }

   int g1 = genes[pbc(position_1)];
   int g2 = genes[pbc(position_2)];
   
   genes[pbc(position_1)] = g2;
   genes[pbc(position_2)] = g1;
}

/// 2) in posizione POS inverte M geni 
void Chromosome::reverse(int position, int m){
   if(position == 0) {
      // std::cerr << "Chromosome::Inversion : can't move first gene: I move the second one." << std::endl;
      position = 1;
   }
   if(position > (n_genes - 1)) {
      // std::cerr << "Chromosome::Inversion : position out of Chromosome dimension. I apply PBC." << std::endl;
      position = pbc(position);
   }
   if(m > n_genes - 1) {
      // std::cerr << "Chromosome::Inversion : too many genes to swap. I set m = 2." << std::endl; 
      m = 2;
   }

   int block[m];

   // metto da parte il blocco da swappare
   for(int i = 0; i < m; ++i)
   {
      block[i] = genes[pbc(position + i)];
   }

   // riempio al contrario il buco lasciato dal blocco
   for(int i = 0; i < m; ++i)
   {
      genes[pbc(position + i)] = block[(m - 1) - i];
   }

   //delete[] block;
}

/// 3) a partire da una posizione POS, sposta M geni adiacenti in avanti di N posizioni, 
///    (eccetto il primo gene e col vincolo m < Ng-1)
void Chromosome::shift(int position, int m, int n){

    if(position == 0) {
        std::cerr << "Chromosome::Shift : can't move first gene: I move the second one." << std::endl;
        position = 1;
    }
    if(position > (n_genes - 1)) {
        std::cerr << "Chromosome::Shift : position out of Chromosome dimension. I apply PBC." << std::endl;
        position = pbc(position);
    }

    if(m >= n_genes - 1) {
        std::cerr << "Chromosome::Shift : too many genes to move. I set m = 1." << std::endl; 
        m = 1;
    }

    int block[m];

    // metto da parte il blocco da shiftare
    for(int i = 0; i < m; ++i){
        block[i] = genes[pbc(position + i)];
    }

    // shifto n volte i geni complementari per riempire il buco...
    for(int i = 0; i < n; ++i){
        genes[pbc(position + i)] = genes[pbc(position + m + i)];
    }

    // ... e rimetto il blocco
    for(int i = 0; i < m; ++i){
        genes[pbc(position + n + i)] = block[i];
    }

    //delete[] block;
}

/// 5) scambia M geni in posizione POS1 con altrettanti in POS2
void Chromosome::permutate(int position_1, int position_2, int m){
    if(position_1 == 0) {
        std::cerr << "Chromosome::MPermut : can't move first gene: I move the second one." << std::endl;
        position_1 = 1;
    }
    if(position_2 == 0) {
        std::cerr << "Chromosome::MPermut : can't move first gene: I move the second one." << std::endl;
        position_2 = 1;
    }
    if(position_1 > (n_genes - 1)) {
        std::cerr << "Chromosome::MPermut : position out of Chromosome dimension. I apply PBC." << std::endl;
        position_1 = pbc(position_1);
    }
    if(position_2 > (n_genes - 1)) {
        std::cerr << "Chromosome::MPermut : position out of Chromosome dimension. I apply PBC." << std::endl;
        position_2 = pbc(position_2);
    }
    if(m > n_genes/2) {
        //std::cerr << "Chromosome::MPermut : too many genes to swap. I set m = 1." << std::endl; 
        m = 1;
    }
    if(m > abs(position_2 - position_1)){ //<< casino
        //std::cerr << "Chromosome::MPermut : m > abs(pos1-pos2). I set m = abs(pos1-pos2)." << std::endl; 
        m = abs(position_1-position_2);
    }
    if(position_1 > position_2){
        int appo = position_2;
        position_2 = position_1;
        position_1 = appo;
    }

    int block[m];

    // metto da parte il blocco da swappare
    for(int i = 0; i < m; ++i){
        block[i] = genes[pbc(position_1 + i)];
    } 

    // riempio il buco lasciato dal blocco...
    for(int i = 0; i < m; ++i){
        genes[pbc(position_1+i)] = genes[pbc(position_2+i)];
    }

    // ... e rimetto il blocco in pos2
    for(int i = 0; i < m; ++i){
        genes[pbc(position_2 + i)] = block[i];
    }

    //delete[] block;
}

/// 6) fa il crossover di len geni a partire dalla posizione pos, con un cromosoma parent2
void Chromosome::crossover(int position, int len, Chromosome chromosome_2){

    if(position == 0) {
        std::cerr << "Chromosome::Crossover : can't move first gene: I move the second one." << std::endl;
        position = 1;
    }

    int block[len];
    //int block2[len];

    // 1. cut their paths at the same position:
    //    metto da parte i blocchi da swappare
    //    (conservando la prima parte)
    for(int i = 0; i < len; ++i){
        block[i] = genes[pbc(position + i)];
        //block2[j] = parent2.Gen[Pbc(pos+j)];
    } 

    // 2. complete the paths with the missing cities adding them in the **order** 
    //    in which they appear in the consort (vale anche se sono separati!):
    int count = 0;
    //int count2 = 0;
    
    for(int j = 1; j < n_genes; ++j){
        // ciclo sui geni nei blocchi 
        for(int i = 0; i < len; ++i){
            // quando trovo nell'altro genitore il gene presente
            // nel blocco, lo salvo nel buco scavato all'inizio
            if(chromosome_2.genes[pbc(j)] == block[i] ){
                genes[pbc(position+count)] = chromosome_2.genes[pbc(j)];
                count++;
            }
            //if(Gen[Pbc(j)] == block2[i] ){ // se il crossover modifica solo il cromosoma parent1, questo ciclo non serve
            //   parent2.Gen[Pbc(pos+count2)] = Gen[Pbc(j)];
            //   count2++;
            //}
        }
    } 
}


/// Costruttore
Task::Task(std::string circ){
    if(circ == "true"){
        this->task = "circ";
    }
    else if (circ == "false"){
        this->task = "square";
    }
    
}
/// Distruttore
Task::~Task(){}

/// genera coordinate di città in su una circonferenza unitaria
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
        std::cerr << "Parameter <circ> must be 'true0 or 'false'. Exit." << std::endl;
        exit(-1);
    }
}

void Task::load_cities(std::string filename){
    std::ifstream coord;
    coord.open(filename);
    
    double pos;

    if(coord.is_open()){
        for(int i = 0; i < n_cities; ++i){
            coord >> pos; // load x
            cities_x[i] = pos;
            coord >> pos; // load y
            cities_y[i] = pos;
        }
    }else{
        std::cout << "Unable to read "+ filename + ". Exit." << std::endl; 
        exit(-1);
    }

    coord.close();
}

/// stampa le coordinate delle città nell'ordine indicato da un cromosoma
void Task::print_cities(int generation, Chromosome chr){

    std::ofstream stream;
    stream.open(std::string(ROOT_PATH) + "/Data/09.1_" + task + "_city_coord_" + intToStringWithLeadingZeros(generation + 1, 3) + ".dat");

    int seq;
    for(int i = 0; i < n_cities; ++i){
        // stampo 
        seq = chr.get_gene(i) - 1; // nell'ordine indicato dal cromosoma
        stream << cities_x[seq]  << " " << cities_y[seq] << std::endl;
    }
    stream.close();
}

/// stampa la lunghezza media della migliore metà di individui in una popolazione
void Task::print_bests_len_ave(int generation, int part, Population pop){
    std::ofstream str;
    str.open(std::string(ROOT_PATH) + "/Data/09.1_" + task + "_best_len_average.dat", std::ios_base::app);

    // ho una variabile BestsLenAve in pop che tiene 
    // conto della media nella metà migliore 
    // (viene riempita quando faccio le mutazioni)

    str << generation + 1 << " " << pop.best_len_ave << std::endl;
    str.close();
}

/// Stampa la lunghezza del miglior individuo di una popolazione
void Task::print_best_len(int generation, Population pop){
    std::ofstream str;
    str.open(std::string(ROOT_PATH) + "/Data/09.1_" + task + "_best_len.dat", std::ios_base::app);

    // ho una variabile BestLen in pop che tiene 
    // conto della lunghezza migliore 
    // (viene riempita quando faccio le mutazioni)

    str << generation  + 1 << " " << pop.best_len << std::endl;
    str.close();
}

/// prende un cromosoma e con l'ordine dei geni calcola la distanza tra le città e la fitness
void Task::eval_fitness(Chromosome &chromosome){
    double fitn = 0;
    double accu = 0;
    double acculen = 0;
    double len2 = 0;
    int ng = chromosome.get_n_genes();

    // segmenti dalla prima all'ultima
    for(int i = 0; i < ng - 1; ++i){

        len2 = pow( cities_x[chromosome.get_gene(i + 1) - 1] - cities_x[chromosome.get_gene(i) - 1], 2) 
            + pow( cities_y[chromosome.get_gene(i + 1) - 1] - cities_y[chromosome.get_gene(i) - 1], 2);
        accu += len2;
        acculen += sqrt(len2);
    }
    // distanza tra la prima e l'ultima
    len2 = pow(cities_x[chromosome.get_gene(0) - 1] - cities_x[chromosome.get_gene(ng - 1) - 1], 2) 
            + pow(cities_y[chromosome.get_gene(0) - 1] - cities_y[chromosome.get_gene(ng - 1) - 1], 2);
    accu += len2;
    acculen += sqrt(len2);

    fitn = 1.0 / accu;
    chromosome.fitness = fitn;
    chromosome.len = acculen;
}

/// prende una popolazione e applica EvalFitness su tutti i suoi individui
void Task::eval(Population &population){ 
    for(int i = 0; i < population.get_n_individuals(); ++i){
        eval_fitness(population.chromosomes[i]);
    }
}

/// prende una popolazione e mette in ordine i suoi individui dal peggiore al migliore
void Task::sort_population(Population *population){
    quickSort(population, 0, population->get_n_individuals()-1);
}


//========================================================//
//       OTHER FUNCTIONS                                  //
//========================================================//

/// funzione necessaria in quicksort
int partition(Population *population, int low, int high){
    
   // prendo la ftn dell'elemento più a destra
   double pivot = population->chromosomes[high].fitness;

   // puntatore all'elemento più grande (ipotizzo sia a sx) 
   int i = (low-1);

   // comparo ogni ftn del vettore con la ftn pivot
   for (int j = low; j < high; j++) {
      if (population->chromosomes[j].fitness <= pivot) {
         // se ftn è minore, scambio l'intero chr con quello puntato da i
         i++;

         // swap
         Chromosome t = population->chromosomes[i];
         population->chromosomes[i] = population->chromosomes[j];
         population->chromosomes[j] = t;
      }
   }

   // scambio pivot con l'elemento puntato da i
   Chromosome tt = population->chromosomes[i+1];
   population->chromosomes[i+1] = population->chromosomes[high];
   population->chromosomes[high] = tt;

   // return the partition point
   return (i + 1);
}

/// funzione necessaria in SortPop
void quickSort(Population *population, int low, int high) {
  if (low < high) {
    int pi = partition(population, low, high);

    quickSort(population, low, pi - 1); 
    quickSort(population, pi + 1, high); 
   }
}

//------------------------Progress bar --------------------------

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
	std::cout << Green(perc.substr(0, 3 * ++s));
	std::cout << Gray(perc.substr(3 * s, 3 * 10));
	std::cout << " " << int(prog * 100.0 / NGeneration) << " %\r";
	std::cout.flush();
}
//------------------------------------------------------------------

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
