#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

#include <random.h>
#include <version_config.h>

//Start with Random Library
/**
 * @brief Initializes the random number generator by reading the prime numbers and the seed from files.
 * @param random_generator a reference to a Random object representing the random number generator to be initialized
 */
void Random_Start(Random &random_generator);

//Exercise_02.1
/**
 * @brief Computes the error of a set of n measurements by calculating the standard deviation of the sample mean.
 * 
 * @param AV a double representing the sample mean of the measurements
 * @param AV2 a double representing the squared sample mean of the measurements
 * @param n an integer representing the number of measurements
 * @return double representing the standard deviation of the sample mean
 */
double error(double AV, double AV2, int n);

/**
 * @brief This method performs a Monte Carlo simulation using the provided 
 * characteristic functions. It calculates the average and its error for a 
 * series of blocks. The random numbers are generated using the provided 
 * random generator and the characteristic function pointers. 
 * The calculated averages and errors are written to the provided output file.

 * 
 * @param n_iteration an integer specifying the total number of iterations.
 * @param n_block an integer specifying the number of blocks.
 * @param random_generator a reference to a Random object used for generating random numbers.
 * @param file_output a reference to an ofstream object representing the output file.
 * @param p_to_function a pointer to a characteristic function that takes a reference to a double and a Random object.
 * @param p_to_function_2 a pointer to a transformation function that takes a reference to a double, a reference to a double, and a reference to an integer.
 */
void blocks_averege(int n_iteration, int n_block, Random &random_generator, 
                    std::ofstream &file_output, void (*p_to_function)(double&, Random &), 
                    void (*p_to_function_2)(double&, double&, int&));

/**
 * @brief This method calculates the sum of the values that the integrand assumes
 * at random x in [0; 1). Integrand is 0.5 * M_PI * cos(M_PI * 0.5 * x).
 * This functions in used as a pointer function for "blocks_average" function.
 *
 * @param sum accumulator for summing values.
 * @param random_generator the Random object for random number generation.
 */
void characteristic_function_uniform_1(double& sum, Random & random_generator);

/**
 * @brief This method calculates the average of the values that the integrand assumes
 * at random x in [0; 1). Integrand is 0.5 * M_PI * cos(M_PI * 0.5 * x).
 * This functions in used as a pointer function for "blocks_average" function.
 *
 * @param sum accumulator for summing values.
 * @param random_generator the Random object for random number generation.
 */
void characteristic_function_uniform_2(double& ave_i, double& sum, int& L);

/**
 * @brief This method calculates the sum of the values that the integrand assumes
 * at random x in [0; 1). Integrand is 0.5 * M_PI * cos(M_PI * 0.5 * x). This function
 * implement importance sampling algorithm.
 * This functions in used as a pointer function for "blocks_average" function.
 *
 * @param sum accumulator for summing values.
 * @param random_generator the Random object for random number generation.
 */
void characteristic_function_sampling_1(double& sum, Random & random_generator);

//Exercise_02.2
/**
 * Simulates a random walk with discrete steps.
 *
 * @param random_generator a Random object for random number generation.
 * @param cartesian_position Cartesian coordinates of the walker's position.
 * @param step_lenght length of each step.
 */
void random_walk_discrete_step(Random &random_generator, std::vector<double> &cartesian_position, double step_lenght);

/**
 * Simulates a random walk with continuous steps.
 *
 * @param random_generator a Random object for random number generation.
 * @param cartesian_position Cartesian coordinates of the walker's position.
 * @param step_lenght length of each step.
 */
void random_walk_continue_step(Random &random_generator, std::vector<double> &cartesian_position, double step_lenght);