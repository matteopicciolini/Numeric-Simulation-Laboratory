#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <string>
#include "random.h"
#include <version_config.h>

Random rnd;  

int nstep = 1e4;
int nblk = 100;
int steps = nstep * nblk;
double L = 1.5;
bool gauss;
int np_rend = 10000;

int mod = 1;
int contastep = 0;

std::string prob_str;


double d = 0;// distance fron the center
double var_prog_d = 0;
double mean_prog_d = 0;

double r[3]; // initial position
double y[3]; // position before the step
double x[3]; // position after the step

double dr_g[3];

double q, A;
int attempted_gs = 0, accepted_gs = 0;


double de = 0;
double var_prog_de = 0;
double mean_prog_de = 0;

double re[3];
double ye[3];
double xe[3];

double dr_e[3];

double qe, Ae;
int attempted_es = 0, accepted_es = 0;


//functions
void Input();

/**
 * @brief Computes the probability of finding the electron in the ground state at a certain position in space.
 * @param x an array of 3 elements representing the cartesian coordinates of the electron
 * @return the probability of finding the electron in the ground state at the specified position
*/
double prob_ground(double*);

/**
 * @brief Computes the probability of finding the electron in the excited state at a certain position in space.
 * @param x an array of 3 elements representing the cartesian coordinates of the electron
 * @return the probability of finding the electron in the excited state at the specified position
*/
double prob_excited(double*);

/**
 * @brief This function resets all the variables used in the Monte Carlo simulation to their initial values.
 * The distances 'd' and 'de' are set to zero. The acceptance rates for the ground state ('attempted_gs' and
 * 'accepted_gs') and excited state ('attempted_es' and 'accepted_es') are reset to zero.
 * @param void
*/
void ResetAllZero();

/**
 * @brief This function saves the positions of the particles in the ground state and excited state to output files.
 * It first opens two output files, one for each state, and appends the current positions to the end of each file.
 * The file names are generated based on the 'pattern' parameter and the ROOT_PATH environment variable.
 * @param void
 */
void SavePos();

/**
 * @brief This function performs a Monte Carlo move by proposing a new step for each particle in the system,
 * calculating the acceptance probability for the ground state and excited state, and updating the
 * positions of the particles based on the acceptance rate.
 * @param void
*/
void Move();

/**
 * @brief The distance of each particle from the origin is calculated using the Euclidean distance formula, and the total
 * distance of all particles in each state is accumulated in the 'd' and 'de' variables, respectively.
 * @param void
 */
void Accumulate();

/**
 * @brief This function takes an integer argument 'i' and prints the acceptance rates for the ground and excited states of a the i-th block.
 * The output is printed to the console in a formatted manner.
 * @param blknum The block number for which the acceptance rates are being printed.
 * @return void
*/
void PrintRate(int blknum);

/**
 * @brief This function calculates the block averages for the distances of the particles from the origin in both the ground
 * state and excited state. The accumulated distances are divided by the number of Monte Carlo steps (nstep) to obtain
 * the average distance for each state. The mean_prog_d and mean_prog_de variables are then updated with the average
 * distances for each state, and the var_prog_d and var_prog_de variables are updated with the square of the average
 * distances for each state.
 * @param void
 */
void BlockAverages();

/**
 * @brief This function saves the distances of the particles from the origin in the ground state and excited state to output files,
 * along with the mean and standard error of the mean for each distance. It first opens two output files, one for each state,
 * and appends the current distances and statistical measures to the end of each file. The file names are generated based on
 * the 'pattern' parameter and the ROOT_PATH environment variable. The 'i' parameter is the current block index, which is used
 * to calculate the standard error of the mean.
 * @param blknum The current block index
 */
void SaveDist(int blknum);

/**
 * @brief Initializes the random number generator by reading the prime numbers and the seed from files.
 * @param random_generator a reference to a Random object representing the random number generator to be initialized
 */
void Random_Start(Random &random_generator);

/**
 * @brief Computes the error on the mean of a quantity for the current block.
 * @param sum the sum of the quantity for the current block
 * @param sum2 the sum of the squares of the quantity for the current block
 * @param iblk the current block number
 * @return the error on the mean of the quantity for the current block, or 0 if the block is the first one
*/
double Error(double sum, double sum2, int iblk);


/**
 * @brief This function proposes a new step for the Monte Carlo simulation. 
 * If 'gauss' is false, it generates a random displacement vector for the
 * particles in the system within the simulation box. Otherwise, it uses
 * a Gaussian distribution to generate the displacement vector.
 * @param void
*/
void ProposeStep();


std::string perc = "▪▪▪▪▪▪▪▪▪▪";

std::string Green(std::string green){
	return std::string("\033[1;32m") + green + "\033[0m";
}

std::string Gray(std::string gray){
	return std::string("\033[1;90m") + gray + "\033[0m";
}

std::string Red(std::string red){
	return std::string("\033[1;31m") + red + "\033[0m";
}

std::string pattern;