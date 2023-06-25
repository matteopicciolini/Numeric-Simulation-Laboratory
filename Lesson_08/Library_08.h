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

void Random_Start(Random &random_generator, int primes_selector = 0);

double EvalPotential(double x, double a = -2.5, double b = 1.);
double EvalWaveFunctionSquared(double x, double mu, double sigma);
double EvalWaveFunctionSecondDerivative(double x, double mu, double sigma);
double EvalWaveFunction(double x, double mu, double sigma);


void Equilibrate(int nblocks, int L, double &initialPosition, Random &rnd, double delta, double mu, double sigma);
void Metropolis(double &initialPosition, Random &rnd, double delta, double &accepted, double &attempted, double mu, double sigma);

//double Error(double AV, double AV2, int n);
double Error(double sum, double sum2, int iblk);

//------------------------Progress bar --------------------------
std::string Green(std::string green);
std::string Gray(std::string gray);
std::string Red(std::string red);
void Progress_bar(int& prog, int N, std::string perc);
//------------------------------------------------------------------