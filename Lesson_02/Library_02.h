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
void blocks_averege(int n_iteration, int n_block, Random &random_generator, 
                    std::ofstream &file_output, void (*p_to_function)(double&, Random &), 
                    void (*p_to_function_2)(double&, double&, int&));
void characteristic_function_uniform_1(double& sum, Random & random_generator);
void characteristic_function_uniform_2(double& ave_i, double& sum, int& L);
void characteristic_function_sampling_1(double& sum, Random & random_generator);

//Exercise_02.2
void random_walk_discrete_step(Random &random_generator, std::vector<double> &cartesian_position, double step_lenght);
void random_walk_continue_step(Random &random_generator, std::vector<double> &cartesian_position, double step_lenght);
void random_walk_fill_matrix_with_step(int n_iteration, int n_blocks, int n_step, Random &random_generator, std::vector<std::vector<double>> &matrix, double step_lenght, std::string simulation_type);
void random_walk_block_averege(int n_blocks, int n_step, std::vector<std::vector<double>> &matrix, std::ofstream &file_output);
