#include "Library_09.h"
#include "Chromosome_test.h"

int main (int argc, char* argv[]){
    //mutation test
    Chromosome_Test chromosome_test;
    chromosome_test.runTests();

    //Random Start
	Random rnd;
	Random_Start(rnd);

    Chromosome chromosome;
    chromosome.fill(rnd);
    chromosome.print();
    
    return 0;
}