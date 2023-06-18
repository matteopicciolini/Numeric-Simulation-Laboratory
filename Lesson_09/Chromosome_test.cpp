#include "Chromosome_test.h"

void Chromosome_Test::runTests() {
    this->mutations_tests();
}

void Chromosome_Test::mutations_tests() {
    Chromosome chromosome(10);
    Chromosome verified_chromosome(10);

    // SET TEST
    int gen[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    chromosome.set_gen(gen);
    verified_chromosome.set_gen(gen);
    assert(chromosome.equals(verified_chromosome));
    
    // SWAP TEST
    int gen_swap[] = {1, 2, 6, 4, 5, 3, 7, 8, 9, 10};
    verified_chromosome.set_gen(gen_swap);
    chromosome.swap(2, 5);
    assert(chromosome.equals(verified_chromosome));
    
    // REVERSE TEST
    int gen_reverse[] = {1, 2, 7, 3, 5, 4, 6, 8, 9, 10};
    verified_chromosome.set_gen(gen_reverse);
    chromosome.reverse(2, 5);
    assert(chromosome.equals(verified_chromosome));

    // PERMUTATE TEST
    int gen_permut[] = {1, 2, 6, 8, 5, 4, 7, 3, 9, 10};
    verified_chromosome.set_gen(gen_permut);
    chromosome.permutate(2, 6, 2);
    assert(chromosome.equals(verified_chromosome));

    // INVERSION TEST
    int gen_shift[] = {1, 2, 7, 3, 9, 6, 8, 5, 4, 10};
    verified_chromosome.set_gen(gen_shift);
    chromosome.shift(2, 4, 3);
    assert(chromosome.equals(verified_chromosome));
}
