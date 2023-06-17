#include "Library_09.h"

int main (int argc, char* argv[]){
    //Random Start
	Random rnd;
	Random_Start(rnd);

    Chromosome chromosome;
    chromosome.fill(rnd);
    chromosome.print();

    chromosome.swap(2, 3);
    chromosome.print();

    chromosome.reverse(1, 3);
    chromosome.print();
    
    chromosome.shift(1, 2, 2);
    chromosome.print();

    return 0;
}