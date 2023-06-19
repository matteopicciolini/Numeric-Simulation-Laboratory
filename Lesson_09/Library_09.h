#pragma once

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <string>
#include <version_config.h>
#include "random.h"

static const int NindPop = 300;      // number of idividuals in population [good:3000]
static const int Ngenes = 34;        // number of genes in a chromosome = number of cities in the problemset
static const int NGeneration = 500; // numero di generazioni consecutive
static const double pmut[6] = {
    0.1, // PairPermut
    0.1, // Invertion
    0.2, // Shift
    0.4, // MPermut
    0.6  // Crossover
    };
static const double expon = 6.0;


class Chromosome {

private:
    static const int n_genes = Ngenes; // regolare da parameters
    int genes[n_genes] = {0};

public:
    double fitness = 0;
    double len = 0;

    // constructors
    Chromosome();
    
    // destructor
    ~Chromosome();

    Chromosome& operator= (const Chromosome& chr);

    // methods

    // set e get
    void set_gen(int vec[n_genes]);
    void set_fitness(double f);
    void set_len(double l);
    double get_fitness();
    double get_len();
    int get_n_genes();
    int get_gene(int i);

    // altri
    void fill(Random &rnd);
    void empty();
    void print();

    // mutations
    void swap(int n, int m);
    void shift(int pos, int m, int n);
    void permutate(int pos1, int pos2, int m);
    void reverse(int pos, int m);
    void crossover(int pos, int len, Chromosome parent2);
    // accessori alle mutazioni
    int pbc(int pos);
    void check();
};


class Population {

private:
    static const int n_individuals = NindPop; // regolare da parameters
    const double *mutation_probability = pmut; // basic mutation probability
public:
    double best_len;
    double best_len_ave;
    Chromosome chromosomes[n_individuals];

    // constructors
    Population();

    // destructor
    ~Population();

    // methods
    int get_n_individuals();
    void set_configuration(Random &rnd);
    void mutate(Random &rnd);
};


class Task {

private:
    static const int n_cities = Ngenes; // regolare da parameters
    double cities_x[n_cities];
    double cities_y[n_cities];
    std::string task;

public:
    Task(std::string circ); // costruttore
    ~Task(); // distruttore

    void generate_circular_cities(Random rnd, double radius = 1.);
    void generate_squared_cities(Random rnd);
    void generate_cities(Random rnd);
    void print_cities(int generation, Chromosome chr);
    void print_bests_len_ave(int generation, int part, Population pop);
    void print_best_len(int generation, Population pop);
    void eval_fitness(Chromosome &chr);
    void eval(Population &pop);
    void sort_population(Population *pop);
    void load_cities(std::string filename);
    
};

void Random_Start(Random &random_generator);
int partition(Population *pop, int low, int high);
void quickSort(Population *pop, int low, int high);
void Delete_old_files(std::string pattern);
void Usage(int argc, char* argv[], std::string &circ, std::string &print_city);
void Progress_bar(int& s, int& prog);
std::string Red(std::string red);
std::string Gray(std::string gray);
std::string Green(std::string green);
std::string intToStringWithLeadingZeros(int value, int width);