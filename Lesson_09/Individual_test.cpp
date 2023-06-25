#include "Individual_test.h"

Individual_test::Individual_test(){
    int gen[34] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
    
    chromosome.set_gene(gen);
    swap_test();

    chromosome.set_gene(gen);
    reverse_test();

    chromosome.set_gene(gen);
    permutate_test();

    chromosome.set_gene(gen);
    shift_test();
}

void Individual_test::swap_test(){
    int gen_swap[34] = {1, 3, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
    Individual correct_chromosome;
    correct_chromosome.set_gene(gen_swap);
    chromosome.swap(1, 2);
    assert(chromosome.equals(correct_chromosome));
};


void Individual_test::reverse_test(){
    int gen_reverse[34] = {1, 2, 7, 6, 5, 4, 3, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
    Individual correct_chromosome;
    correct_chromosome.set_gene(gen_reverse);
    chromosome.reverse(2, 5);
    assert(chromosome.equals(correct_chromosome));
};

void Individual_test::permutate_test(){
    int gen_permut[34] = {1, 2, 7, 8, 5, 6, 3, 4, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
    Individual correct_chromosome;
    correct_chromosome.set_gene(gen_permut);
    chromosome.permutate(2, 6, 2);
    assert(chromosome.equals(correct_chromosome));
};

void Individual_test::shift_test(){
    int gen_shift[34] = {1, 2, 7, 8, 9, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
    Individual correct_chromosome;
    correct_chromosome.set_gene(gen_shift);
    chromosome.shift(2, 4, 3);
    assert(chromosome.equals(correct_chromosome));
};
