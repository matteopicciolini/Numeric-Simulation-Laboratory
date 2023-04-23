#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

#include <random.h>
#include <version_config.h>

//Start with Random Library
void Random_Start(Random &random_generator);

//Exercise_01.a
double error(double AV, double AV2, int n);