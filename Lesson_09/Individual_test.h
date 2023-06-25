#pragma once

#include "Library_09.h"
#include <cassert>

class Individual_test{
private:
    Individual chromosome;
public:
    Individual_test();
    ~Individual_test(){};

    void swap_test();
    void reverse_test();
    void permutate_test();
    void shift_test();
};