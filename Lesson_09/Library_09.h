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

//COSTANTI
static const int NindPop = 1000; //Numero di individui per popolazione
static const int Ngenes = 34; // numero di città
static const int NGeneration = 500; // numero di generazioni
static const double Mutation_Probability[5] = {
    0.11, // swap
    0.1, // reverse
    0.05, // shift
    0.05, // permutate
    0.8  // crossover
};
static const double expon = 6.0;

// CLASSI
class Individual{

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

    /**
     * @brief Fills the individual's genes with unique values using a random number generator.
     * If the first gene is already filled, it displays an alert and empties the individual.
     * 
     * @param rnd The random number generator used to generate random numbers.
     */
    void fill(Random &rnd);
    /**
     * @brief Empties the individual's genes by setting them to zero.
     */
    void empty();
    void print();
    int pbc(int pos);
    void check();
    bool equals(Individual& individual_2);

    // mutations
    /**
     * @brief Swaps the values of two positions in the individual's genes array.
     *        If either position is 0, it is treated as position 1.
     *
     * @param position_1 The first position to swap.
     * @param position_2 The second position to swap.
     */
    void swap(int position_1, int position_2);

    /**
     * @brief Shifts a block of genes within the individual's genes array.
     *        If the position is 0, it is treated as position 1.
     *        If the position exceeds the maximum index, it is adjusted using periodic boundary conditions.
     *        If the block size (m) is equal to or greater than n_genes - 1, m is set to 1.
     *
     * @param position The starting position of the block to shift.
     * @param m The size of the block to shift.
     * @param n The number of times to shift the block.
     */
    void shift(int position, int m, int n);

    /**
     * @brief Permutes a block of genes between two positions in the individual's genes array.
     *        If either position is 0, it is treated as position 1.
     *        If a position exceeds the maximum index, it is adjusted using periodic boundary conditions.
     *        If m is greater than n_genes/2, m is set to 1.
     *        If m is greater than the distance between position_1 and position_2, m is set to the distance.
     *        If position_1 is greater than position_2, they are swapped.
     *
     * @param position_1 The first position defining the block to permute.
     * @param position_2 The second position defining the block to permute.
     * @param m The size of the block to permute.
     */
    void permutate(int position_1, int position_2, int m);

    /**
     * @brief Reverses a block of genes within the individual's genes array.
     *        If the position is 0, it is treated as position 1.
     *        If the position exceeds the maximum index, it is adjusted using periodic boundary conditions.
     *        If m is greater than n_genes - 1, m is set to 2.
     *
     * @param position The starting position of the block to reverse.
     * @param m The size of the block to reverse.
     */
    void reverse(int position, int m);

    /**
     * @brief Performs crossover between the individual's genes and genes from another individual.
     *        The crossover is performed starting from the specified position and for a specified length.
     *        If the position is 0, it is treated as position 1.
     *
     * @param position The starting position of the crossover.
     * @param len The length of the crossover block.
     * @param individual_2 The other individual to perform crossover with.
     */
    void crossover(int position, int len, Individual individual_2);
};


class Population {

private:
    static const int n_individuals = NindPop;
    const double *mutation_probability = Mutation_Probability;
public:
    double best_len;
    double best_len_ave;
    Individual individuals[n_individuals];


    Population(Random &rnd);
    ~Population();

    int get_n_individuals();

    /**
     * @brief Reproduces the population by generating new individuals through mutation and crossover.
     *
     * This function performs reproduction on the population by generating new individuals based on the existing individuals.
     * It selects a portion of the best individuals according to a power law distribution, and saves their genetic information.
     * The best length and average length of the population are also calculated during this process.
     * Mutation operations (swap, reverse, shift, permutate) and crossover are then applied to the new individuals to introduce genetic variations.
     * The mutated individuals are checked for validity and assigned to replace the original individuals in the population.
     *
     * @param rnd A reference to the random number generator object.
     */
    void reproduce(Random &rnd);
};


class Task {

private:
    static const int n_cities = Ngenes;
    double cities_x[n_cities];
    double cities_y[n_cities];
    std::string task;

public:
    Task(std::string circ);
    ~Task();

    void generate_circular_cities(Random rnd, double radius = 1.);
    void generate_squared_cities(Random rnd);
    void generate_cities(Random rnd);

    void print_cities(int generation, Individual individual);
    void print_bests_len_ave(int generation, int part, Population pop);
    void print_best_len(int generation, Population pop);

    /**
     * @brief Calculates the fitness of an individual based on the order of its genes, 
     * considering the distances between cities.
     * The fitness is determined by the inverse of the total distance traveled. The shorter the distance, the higher the fitness.
     *
     * @param individual The individual for which to calculate the fitness.
     */
    void eval_fitness(Individual &individual);

    /**
     * @brief Evaluates the fitness of all individuals in a population by applying eval_fitness function to each individual.
     *
     * @param population The population to evaluate.
     */
    void eval(Population &population);

    /**
     * @brief Sorts the individuals in a population in descending order based on their fitness.
     *
     * @param population The population to evaluate.
     */
    void sort_population(Population &population);
};

// FUNZIONI
void swap(Individual& a, Individual& b);
void Random_Start(Random &random_generator);
int partition(Population &population, int low, int high);
void quickSort(Population &population, int low, int high);
void Delete_old_files(std::string pattern);
void Usage(int argc, char* argv[], std::string &circ, std::string &print_city);
void Progress_bar(int& s, int& prog);
std::string Red(std::string red);
std::string Gray(std::string gray);
std::string Green(std::string green);
std::string intToStringWithLeadingZeros(int value, int width);