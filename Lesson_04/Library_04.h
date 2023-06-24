/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

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

int seed[4];
Random rnd;
int nconf = 1;

//parameters, observables
const int m_props = 1000;
int n_props, iv, ik, it, ie, iw;
double vtail, ptail, bin_size, nbins, sd;
double walker[m_props];

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_pres, stima_kin, stima_etot, stima_temp;
double err_pot, err_press, err_kin, err_etot, err_temp, err_gdir;

//configuration
const int m_part = 108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// thermodynamical state
int npart;
double beta, temp, energy, vol, rho, box, rcut;

// simulation
int iNVET, nstep, nblk, restart;
double delta;

//pigreco
const double pi = 3.1415927;

//functions

/**
 * @brief Reads the input parameters and initializes the simulation.
 * This function reads the input parameters from the appropriate files and initializes the simulation.
 * It sets up the necessary variables, arrays, and random number generator seed.
 * The input files used are "Primes" for random number generator primes,
 * "seed.in" or "seed.out" for the random number generator seed,
 * and "config.in" for the initial configuration.
 * The function also calculates and prints some initial values for measured properties.
 */
void Input(void);

/**
 * @brief Resets the block averages and other variables for a new block.
 * This function resets the block averages, as well as other variables 
 * such as the number of attempted and accepted moves,
 * to prepare for a new block of measurements.
 * If it is the first block (iblk == 1), it also resets the global averages.
 * 
 * @param iblk The index of the current block.
 * 
*/
void Reset(int);


/**
 * @brief Updates the block averages with the current measurements.
 * This function updates the block averages with the current measurements.
 * It adds the values of the measured properties in the walker array to the blk_av array,
 * and increments the blk_norm variable by 1.
*/
void Accumulate(void);

/**
 * @brief Computes the averages of the measured properties for the current block.
 * 
 * This function computes the averages of the measured properties for the current block
 * by dividing the accumulated values in the `blk_av` array by the `blk_norm` variable.
 * It also computes the errors on the averages using the blocking method.
 * The computed averages and errors are stored in the `ave` and `err` arrays, respectively.
 */
void Averages(int);

/**
 * @brief Performs a Monte Carlo or Molecular Dynamics move.
 * This function performs a Monte Carlo move or a Molecular Dynamics move
 * depending on the value of the flag iNVET.
 * In the Molecular Dynamics move, the Verlet integration scheme is used to update the particle positions and velocities.
 * Forces acting on each particle are computed using the Lennard-Jones potential.
 * 
 */
void Move(void);

/**
 * @brief Writes the final configuration to a file.
 * 
 * This function writes the final configuration of the system to a file named "config.final".
 * The configuration includes the positions of all particles.
 */
void ConfFinal(void);
void ConfXYZ(int);

/**
 * @brief Measures the properties of the system.
 * This function measures the properties of the system, including the potential energy, kinetic energy,
 * temperature, total energy, and pressure.
 * The measured values are stored in the walker array.
*/
void Measure(void);

/**
 * @brief Computes the potential energy of a particle.
 * Given the position of a particle ip and the coordinates of all other particles,
 * this function computes the potential energy of the particle using the Lennard-Jones potential.
 * It returns the computed potential energy.
 * @param xx The x-coordinate of the particle.
 * @param yy The y-coordinate of the particle.
 * @param zz The z-coordinate of the particle.
 * @param ip The index of the particle.
 * @return The potential energy of the particle.
*/
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);

/**
 * @brief Computes the force acting on a particle.
 * Given the index of a particle ip and the direction idir (0 for x, 1 for y, 2 for z),
 * this function computes the force acting on the particle in the specified direction.
 * The force is computed as the negative gradient of the potential energy.
 * It returns the computed force.
 * @param ip The index of the particle.
 * @param idir The direction of the force (0 for x, 1 for y, 2 for z).
 * @return The force acting on the particle in the specified direction.
*/

double Force(int, int);
void Usage(int argc, char* argv[]);
void Delete_old_files();


//Usage and filename
std::string phase;
std::string eq;
std::string random_lib_path = std::string(ROOT_PATH) + "/random-library/";
std::string input_path, output_path, eq_str, input_form_eq_path;
std::string pattern;


//------------------------Progress bar --------------------------
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

void Progress_bar(int& s, int& iblk, int& prog){
	std::cout << Red("Block number " + std::to_string(iblk) + ". Progress: ");
	std::cout << Green(perc.substr(0, 3 * ++s));
	std::cout << Gray(perc.substr(3 * s, 3 * 10));
	std::cout << " " << int(prog * 100.0 / nstep) << " %\r";
	std::cout.flush();
}
//------------------------------------------------------------------


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/