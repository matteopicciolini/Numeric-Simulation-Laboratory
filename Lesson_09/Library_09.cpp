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

Chromosome::Chromosome(){}
Chromosome::~Chromosome(){}

void Chromosome::set_gen(int vec[n_genes]){
    for(int i = 0; i < this->n_genes; ++i) this->gen[i] = vec[i];
}

void Chromosome::set_fitness(double fitness){ 
    this->fitness = fitness;
}

void Chromosome::set_len(double len){ 
    this->len = len;
}

double Chromosome::get_fitness(){
    return this->fitness;
}

double Chromosome::get_len(){ 
    return this->len;
}

int Chromosome::get_n_genes(){
    return this->n_genes;
}

int Chromosome::get_gen(int i){
    return this->gen[i];
}

void Chromosome::print(){
    for(int i = 0; i < this->n_genes; ++i) std::cout << gen[i] << " ";
    std::cout << std::endl;
}

void Chromosome::fill(Random &rnd){
    //if(this->gen[0] != 0) {std::cerr << "Can't fill full chromosome: I clean it." << std::endl; this->Empty();}
    int random = (int)(rnd.Rannyu() * this->n_genes); //numero random tra 0 e n_genes
    this->gen[0] = 1; //inserisco nella prima posizione il valore 1
    for(int i = 2; i < this->n_genes + 1; ++i){
        while(this->gen[random] != 0){ // finchè non trovo una posizione vuota in gen, genero un nuovo numero random
            random = (int)(rnd.Rannyu() * this->n_genes);
        }
        this->gen[random] = i; //riempio la posizione
    }
}

void Chromosome::empty(){
    for(int i = 0; i < n_genes; ++i){
        this->gen[i] = 0; 
    }
}

void Chromosome::swap(int position1, int position2){
    if(position1 == 0) {
        std::cerr << "Chromosome::permute : can't move first gene." << std::endl;
        exit(0);
    }
   
    if(position2 == 0) {
        std::cerr << "Chromosome::permute : can't move first gene." << std::endl;
        exit(0);
    }

    int old_value1 = this->gen[this->pbc(position1)];
    int old_value2 = this->gen[this->pbc(position2)];
   
    gen[this->pbc(position1)] = old_value2;
    gen[this->pbc(position2)] = old_value1;
}

void Chromosome::reverse(int position, int m){
    if(position == 0) {
        std::cerr << "Chromosome::reverse: can't move first gene." << std::endl;
        exit(0);
    }
    if(position > this->n_genes-1) {
        std::cerr << "Chromosome::reverse: position out of Chromosome dimension." << std::endl;
        exit(0);
    }
    if(m > this->n_genes-1) {
        std::cerr << "Chromosome::reverse : too many genes to swap." << std::endl; 
        exit(0);
    }

    int block[m];

    // salvo il blocco da invertire
    for(int j = 0; j < m; ++j){
        block[j] = this->gen[this->pbc(position + j)];
    }

    // riempio al contrario il buco lasciato dal blocco
    for(int i = 0; i < m; ++i){
        this->gen[this->pbc(position + i)] = block[(m - 1) - i];
    }
}

void Chromosome::shift(int position, int m, int shift){

    if(position == 0) {
        std::cerr << "Chromosome::shift : can't move first gene." << std::endl;
        exit(0);
    }
    if(position > this->n_genes - 1) {
        std::cerr << "Chromosome::shift : position out of Chromosome dimension." << std::endl;
        exit(0);
    }
    if(m >= this->n_genes - 1) {
        std::cerr << "Chromosome::shift : too many genes to move." << std::endl; 
        exit(0);
    }
    int block[m];
    // [1 4 3 2 5]
    // salvo il blocco da shiftare
    //[4, 3]
    for(int j = 0; j < m; ++j){
        block[j] = this->gen[this->pbc(position + j)];
    }

    // [1 2 5 2 5]
    // shifto n volte i geni complementari per riempire il buco...
    for(int i = 0; i < shift; ++i){
        this->gen[this->pbc(position + i)] = this->gen[this->pbc(position + m + i)];
    }

    this->print();

    // ... e rimetto il blocco
    //[1 2 5 4 3]
    for(int i = 0; i < m; ++i){
        this->gen[this->pbc(position + shift + i)] = block[i];
    }
}



int Chromosome::pbc(int position){
    //se sfora e non è multiplo della lunghezza restituisco il modulo
    if((position > this->n_genes-1) && (position % this->n_genes-1 != 0)){
        position = position % (this->n_genes-1);
    } 
    //se sfora ed è multiplo della lunghezza restituisco la lunghezza
    if((position > this->n_genes-1) && (position % this->n_genes-1 == 0)){
        position = this->n_genes-1;
    }
    if(position < 0) exit(0);
    return position;
}

void Chromosome::permutate(int position1, int position2, int m){
    if(position1 == 0) {
        std::cerr << "Chromosome::permutate: can't move first gene." << std::endl;
        exit(0);
    }
    if(position2 == 0) {
        std::cerr << "Chromosome::permutate: can't move first gene." << std::endl;
        exit(0);
    }
    if(position1 > this->n_genes - 1) {
        std::cerr << "Chromosome::permutate: position out of Chromosome dimension." << std::endl;
        exit(0);
    }
    if(position2 > this->n_genes - 1) {
        std::cerr << "Chromosome::permutate: position out of Chromosome dimension." << std::endl;
        exit(0);
    }
    if(m > this->n_genes / 2) {
        std::cerr << "Chromosome::permutate: too many genes to swap." << std::endl; 
        exit(0);
    }
    if(m > abs(position2 - position1)){ 
        std::cerr << "Chromosome::permutate: m > abs(position1 - position2)." << std::endl; 
    }
    if(position1 > position2){
        int appo = position2;
        position2 = position1;
        position1 = appo;
    }

    int block[m];

    // metto da parte il blocco da swappare
    for(int j = 0; j < m; ++j){
        block[j] = this->gen[this->pbc(position1 + j)];
    } 

    // riempio il buco lasciato dal blocco...
    for(int i = 0; i < m; ++i){
        this->gen[this->pbc(position1 + i)] = this->gen[this->pbc(position2 + i)];
    }

    // ... e rimetto il blocco in pos2
    for(int i = 0; i < m; ++i){
        this->gen[this->pbc(position2 + i)] = block[i];
    }
}