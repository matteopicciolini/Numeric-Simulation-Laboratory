#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>

#include <random.h>
#include <version_config.h>

//Start with Random Library
/**
 * @brief Initializes the random number generator by reading the prime numbers and the seed from files.
 * @param random_generator a reference to a Random object representing the random number generator to be initialized
 */
void Random_Start(Random &random_generator);


/**
 * @brief Computes the error of a set of n measurements by calculating the standard deviation of the sample mean.
 * 
 * @param AV a double representing the sample mean of the measurements
 * @param AV2 a double representing the squared sample mean of the measurements
 * @param n an integer representing the number of measurements
 * @return double representing the standard deviation of the sample mean
 */
double error(double AV, double AV2, int n);