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

//Data blocking
unsigned int M = 100000;
unsigned int N = 100;
unsigned int L = M / N;
double sum_prog = 0, sum_prog_2 = 0, blk_ave = 0;

//Parameters
int i = 0, counter = 0;
double mu = 1;
double sigma = 0.5;
double step = 1.2;
double mu_old;
double sigma_old;
double T_max = 2;

//Variables
double x;
double x_old = 0;
double energy_old;
double energy_new;
double sigma_energy_old;
double sigma_energy_new;
double E = 0;

Random rnd;

//Probabilities
double A = 0;

//Files
ofstream out_1, out_2, out_3, out_4;

//Functions
void Initialization();
void Evolution();
void Move();
void A_0();
void Energy();
void Reset();
double psi(double);
void eval_energy();
void Random_Start(Random &random_generator);