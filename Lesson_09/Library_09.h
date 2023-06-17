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

static const int N_genes = 5;
static const double Mutation_probability[6] = {
    0.1, // PairPermut
    0.1, // Invertion
    0.2, // Shift
    0.2, // Shift2
    0.4, // MPermut
    0.6  // Crossover
};

void Random_Start(Random &random_generator);


class Chromosome {

private:
  static const int n_genes = N_genes; // regolare da parameters
  int gen[n_genes] = {0};

public:
  double fitness = 0;
  double len = 0;
  const double *mutation_probability = Mutation_probability; // basic mutation probability

  // constructors
  Chromosome();

  // destructor
  ~Chromosome();

  Chromosome& operator= (const Chromosome& chr);

  // methods

  // set e get
  void set_gen(int vec[n_genes]);
  void set_fitness(double fitness);
  void set_len(double len);
  double get_fitness();
  double get_len();
  int get_n_genes();
  int get_gen(int i);

  // altri
  void fill(Random &rnd);
  void empty();
  void print();

  // mutations
  void swap(int position1, int position2);
  void shift(int position, int m, int shift);
  void Shift2(int shift);
  void permutate(int pos1, int pos2, int m);
  void reverse(int position, int m);
  void crossover(int pos, int len, Chromosome parent2);
  // accessori alle mutazioni
  int pbc(int pos);
};