#pragma once

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <string>
#include "random.h"
#include <version_config.h>

//Random numbers
int seed[4];
Random rnd;

//parameters, observables
const int m_props = 1000;
int n_props, iu, ic, im, ix, ig;
double nbins;
double walker[m_props];

// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_u, stima_c, stima_m, stima_x, stima_g;
double err_u,err_c, err_m, err_x, err_g;

//configuration
const int m_spin = 50;
double s[m_spin];

// thermodynamical state
int nspin;
double Beta, temp, J, h;

// simulation
int nstep, nblk, metro;

//functions
void Usage(int argc, char* argv[]);
void Input(void);

/**
 * @brief Deletes old files with a specific pattern in the "Data" directory. 
 * The pattern is defined by the input arguments and consists of a prefix followed by three floating point numbers.
 * The function prompts the user to confirm the deletion before proceeding.
 * @param void
 */
void Delete_old_files();

/**
 * @brief Reset the block averages for a given block index. If the block index is 1, reset also the global averages.
 * @param iblk The block index to reset the averages for.
*/
void Reset(int);

/**
 * @brief Updates the block averages of the physical observables.
 * The current values of the observables in walker are added to the corresponding values in blk_av.
 * The normalization factor blk_norm is incremented by 1. 
 */
void Accumulate(void);

/**
 * @brief Computes and prints the block averages and their errors for the current block.
 * If the magnetic field h is zero, the energy, heat capacity, and magnetic susceptibility are computed and written to output files.
 * If the magnetic field is non-zero, the magnetization is computed and written to an output file.
 * The results are printed to the console.
 * The output files are written in append mode.
 * @param iblk The index of the current block.
 */
void Averages(int);

/**
 * @brief Performs a single Monte Carlo move on the lattice.
 * The move consists of selecting a spin at random and attempting to flip it.
 * The type of move is determined by the metro parameter:
 *      metro = 1: Metropolis algorithm is used to accept or reject the move based on the Boltzmann factor.
 *      metro = 0: Glauber algorithm is used to accept or reject the move based on the transition probability.
 * The number of attempted and accepted moves are updated accordingly.
 * @param metro The type of move to use (1 for Metropolis, 0 for Glauber).
 */
void Move(int);

/**
 * @brief Outputs the final configuration of the system to a file and saves the random number generator seed.
 * The final configuration is written to a file named "config.final" in the directory specified by random_lib_path.
 * The output file is written in text mode.
 */
void ConfFinal(void);

/**
 * @brief Computes the physical observables of the current state of the system.
 * The energy u and magnetization m are computed by cycling over all spins.
 * The computed values are stored in walker for later use in the simulation.
 * The energy squared u*u and magnetization squared m*m are also computed and stored in walker.
 */
void Measure(void);

/**
 * Computes the energy contribution of a spin flip according to the Boltzmann factor.
 * The energy of the system with the spin at position ip flipped to sm is computed using the Ising model Hamiltonian.
 * The energy contribution is computed using the Boltzmann factor.
 * @param sm The new spin value (-1 or +1) after the flip.
 * @param ip The index of the spin to be flipped.
 * @return The energy contribution of the spin flip according to the Boltzmann factor.
*/
double Boltzmann(int, int);

/**
 * @brief Implements periodic boundary conditions for a given index.
 * If the index is larger than or equal to the number of spins nspin, it is shifted back by nspin.
 * If the index is negative, it is shifted forward by nspin.
 * @param i The index to apply periodic boundary conditions to.
 * @return int The shifted index, after applying periodic boundary conditions.
 */
int Pbc(int);

/**
 * @brief Computes the statistical error of a block average.
 * The error is estimated using the blocking method, assuming that the data is uncorrelated.
 * @param sum The sum of the values in the block.
 * @param sum2 The sum of the squares of the values in the block.
 * @param iblk The index of the block.
 * @return The estimated error for the block average.
*/
double Error(double, double, int);

/**
 * @brief Writes the simulation results to output files for plotting.
 * The output files depend on the simulation parameters and the chosen pattern.
 * Depending on the value of the magnetic field h, different observables are written:
 *      If h = 0, the energy, heat capacity, and magnetic susceptibility are written.
 *      If h != 0, the magnetization and magnetic field are written.
 * The output files are written in append mode.
 * 
 */
void Results();