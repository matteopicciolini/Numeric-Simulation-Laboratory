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
void average_calc(int n_iteration, int n_block, std::vector<double> &ave, std::vector<double> &ave_2, Random &random_generator);
void blocks_averege(int n_iteration, int n_block, std::vector<double> &ave, std::vector<double> &ave_2, std::ofstream &file_output);

//  Method 2
void characteristic_function_simple(double &sum, Random &random_generator);
void blocks_averege(int n_iteration, int n_block, Random &random_generator, std::ofstream &file_output, void (*p_to_function)(double&, Random &));

//Exercise_01.1b
void characteristic_function_sigma(double &sum, Random &random_generator);

//Exercise_01.2
double exponential_distribution(double lambda, double random);
double Lorentz_distribution(double mean, double gamma, double random);

//Exercise_01.3
void blocks_averege(int n_iteration, int n_block, Random &random_generator, std::ofstream &file_output, void (*p_to_function)(double&, Random &), void (*p_to_function_2)(double&, double&, int&));
void characteristic_function_buffon(double& nhint, Random & random_generator);
void characteristic_function_buffon_2(double& ave_i, double& n_hint, int& L);