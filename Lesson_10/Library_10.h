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
#include <mpi.h>
#include <chrono>
#include <thread>


//COSTANTI
static const int NindPop = 5000; //Numero di individui per popolazione
static const int Ngenes = 50; // numero di città
static const int NGeneration = 500; // numero di generazioni
static const double Mutation_Probability[5] = {
    0.1, // swap
    0.1, // reverse
    0.2, // shift
    0.2, // permutate
    0.5  // crossover
    };
static const double expon = 6.0;

// CLASSI

class Individual {

private:
    static const int n_genes = Ngenes;
    int genes[n_genes] = {0};

public:
    double fitness = 0;
    double len = 0;

    Individual();
    ~Individual();

    Individual& operator= (const Individual& chr);


    void set_gene(int vec[n_genes]);
    void set_fitness(double f);
    void set_len(double l);
    double get_fitness();
    double get_len();
    int get_n_genes();
    int get_gene(int i);

    void fill(Random &rnd);
    void empty();
    void print();
    int pbc(int pos);
    void check();

    // mutations
    void swap(int n, int m);
    void shift(int pos, int m, int n);
    void permutate(int pos1, int pos2, int m);
    void reverse(int pos, int m);
    void crossover(int pos, int len, Individual parent2);
};


class Population {

private:
    static const int n_individuals = NindPop;
    const double *mutation_probability = Mutation_Probability;
public:
    double best_len;
    double best_len_ave;
    Individual individuals[n_individuals];


    Population();
    ~Population();

    int get_n_individuals();
    void set_configuration(Random &rnd);
    void reproduce(Random &rnd);
};


class Task {

private:
    static const int n_cities = Ngenes;
    double cities_x[n_cities];
    double cities_y[n_cities];

public:
    Task();
    ~Task();

    void generate_circular_cities(Random rnd, double radius = 1.);
    void generate_squared_cities(Random rnd);

    void print_cities(int generation, Individual chromosome, int rank, std::string migr);
    void print_bests_len_ave(int generation, int part, Population population, int rank, std::string migr);
    void print_best_len(int generation, Population population, int rank, std::string migr);

    void eval_fitness(Individual &chromosome);
    void eval(Population &population);
    void sort_population(Population *population);
    void load_cities(std::string filename);
};

// FUNZIONI
void swap(Individual& a, Individual& b);
void Random_Start(Random &random_generator, int primes_selector = 0);
int partition(Population *population, int low, int high);
void quickSort(Population *population, int low, int high);
void Delete_old_files(std::string pattern);
void Usage(int argc, char* argv[], std::string &print_city, std::string &migrate, int &migration_freq);
void Progress_bar(int& s, int& prog, int rank);
std::string Red(std::string red);
std::string Gray(std::string gray);
std::string Green(std::string green);
std::string intToStringWithLeadingZeros(int value, int width);