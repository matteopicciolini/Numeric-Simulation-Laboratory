#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

#include <random.h>
#include <version_config.h>

//Start with Random Library
/**
 * @brief Initializes the random number generator by reading the prime numbers and the seed from files.
 * @param random_generator a reference to a Random object representing the random number generator to be initialized
 */
void Random_Start(Random &random_generator);

//Exercise_01.1a
/**
 * @brief Computes the error of a set of n measurements by calculating the standard deviation of the sample mean.
 * 
 * @param AV a double representing the sample mean of the measurements
 * @param AV2 a double representing the squared sample mean of the measurements
 * @param n an integer representing the number of measurements
 * @return double representing the standard deviation of the sample mean
 */
double error(double AV, double AV2, int n);

//  Method 1
/**
 * @brief This method calculates the average and squared average of a 
 * random variable using a Monte Carlo method. It performs a specified 
 * number of iterations divided into blocks. For each block, it generates 
 * random numbers using the provided random generator and calculates the sum. 
 * The average and squared average of each block are stored in the given vectors.
 * 
 * @param n_iteration an integer specifying the total number of iterations.
 * @param n_block an integer specifying the number of blocks to divide the iterations into.
 * @param ave a reference to a vector of doubles that will store the calculated averages.
 * @param ave_2 a reference to a vector of doubles that will store the calculated squared averages.
 * @param random_generator a reference to a Random object used for generating random numbers.
 */

void average_calc(int n_iteration, int n_block, std::vector<double> &ave, std::vector<double> &ave_2, Random &random_generator);

/**
 * @brief This method calculates the progressive average and its error for a series of blocks. 
 * It takes the calculated averages and squared averages stored in the given vectors and 
 * calculates the progressive sum for each block. 
 * The progressive average and its error are then written to the provided output file.
 * 
 * @param n_iteration an integer specifying the total number of iterations.
 * @param n_block an integer specifying the number of blocks.
 * @param ave a reference to a vector of doubles containing the calculated averages.
 * @param ave_2 a reference to a vector of doubles containing the calculated squared averages.
 * @param file_output a reference to an ofstream object representing the output file.
*/
void blocks_averege(int n_iteration, int n_block, std::vector<double> &ave, std::vector<double> &ave_2, std::ofstream &file_output);


// Method 2
/**
 * @brief This method calculates the sum of random numbers generated using 
 * the provided random generator. 
 * It is a simple characteristic function used in the Monte Carlo simulation.
 * @param sum a reference to a double representing the sum of random numbers (updated by the method).
 * @param random_generator a reference to a Random object used for generating random numbers.

*/
void characteristic_function_simple(double &sum, Random &random_generator);

/**
 * @brief This method performs a Monte Carlo simulation using a provided characteristic
 *  function. It calculates the average and its error for a series of blocks. 
 * The random numbers are generated using the provided random generator and a 
 * characteristic function pointer. 
 * The calculated averages and errors are written to the provided output file.
 * 
 * @param n_iteration an integer specifying the total number of iterations.
 * @param n_block an integer specifying the number of blocks.
 * @param random_generator a reference to a Random object used for generating random numbers.
 * @param file_output a reference to an ofstream object representing the output file.
 * @param p_to_function a pointer to a characteristic function that takes a reference to a double and a Random object.
*/
void blocks_averege(int n_iteration, int n_block, Random &random_generator, std::ofstream &file_output, void (*p_to_function)(double&, Random &));

/**
 * @brief This method calculates the sum of squared deviations from a random 
 * number generated using the provided random generator. 
 * It is a characteristic function used in the Monte Carlo simulation.
 * 
 * @param sum a reference to a double representing the sum of squared deviations (updated by the method).
 * @param random_generator a reference to a Random object used for generating random numbers.
 */
void characteristic_function_sigma(double &sum, Random &random_generator);

//Exercise_01.2

/**
 * @brief This method calculates a random number following an exponential 
 * distribution given a lambda value and a random number in the range [0, 1). 
 * It is used as a transformation function in the Monte Carlo simulation.
 * 
 * @param lambda a double specifying the lambda value of the exponential distribution.
 * @param random a double representing a random number in the range [0, 1).
 * @return double representing a random number following an exponential distribution.
 */
double exponential_distribution(double lambda, double random);

/**
 * @brief This method calculates a random number following a Lorentz distribution given
 *  a mean, gamma, and a random number in the range [0, 1). It is used as a 
 * transformation function in the Monte Carlo simulation.
 * 
 * @param mean a double specifying the mean of the Lorentz distribution.
 * @param gamma a double specifying the gamma value of the Lorentz distribution.
 * @param random a double representing a random number in the range [0, 1).
 * @return double a double representing a random number following a Lorentz distribution.
 */
double Lorentz_distribution(double mean, double gamma, double random);

//Exercise_01.3

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
void blocks_averege(int n_iteration, int n_block, Random &random_generator, std::ofstream &file_output, void (*p_to_function)(double&, Random &), void (*p_to_function_2)(double&, double&, int&));

/**
 * @brief This method performs a Buffon's needle experiment using the provided 
 * random generator. It calculates the number of successful hints (when the needle 
 * crosses a line) given the current number of hints. 
 * The number of successful hints is updated accordingly.
 * 
 * @param nhint a reference to a double representing the number of hints (updated by the method).
 * @param random_generator a reference to a Random object used for generating random numbers.
 */
void characteristic_function_buffon(double& nhint, Random & random_generator);

/**
 * @brief This method calculates the average distance between successful hints 
 * in a Buffon's needle experiment given the total distance (sum of distances) 
 * and the number of hints. 
 * It updates the average distance and resets the number of hints.
 * 
 * @param ave_i a reference to a double representing the average distance (updated by the method).
 * @param n_hint a reference to a double representing the number of hints (reset by the method).
 * @param L a reference to an integer representing the total distance.
 */
void characteristic_function_buffon_2(double& ave_i, double& n_hint, int& L);










