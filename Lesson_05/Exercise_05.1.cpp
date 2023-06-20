#include "Library_05.h"

int main (int argc, char* argv[]){

    //Usage and phase choosing
    //Usage(argc, argv);
    //std::cout << std::endl;

    Input();
    std::cout << std::endl;
    
    //equilibration
    for(int i = 0; i < 600000; ++i){
        Move();
    }

    for(int i = 0; i < nblk; ++i){     
        ResetAllZero();
        for (int j = 0; j < nstep; ++j){           
            SavePos();    
            Move();
            Accumulate();
        }            
        PrintRate(i);
        BlockAverages();
        SaveDist(i);

    }

    return 0;
}



void Input(){
  std::ifstream ReadInput;

  std::cout << "Hydrogen atom" << std::endl;

  Random_Start(rnd);
  
  ReadInput.open("input-output/input.dat");
    
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> L;
  ReadInput >> r[0] >> r[1] >> r[2];
  ReadInput >> re[0] >> re[1] >> re[2];
  ReadInput >> gauss;

  if (!gauss){
    prob_str = "Unif";
  }
  else if (gauss){
    prob_str = "Gauss";
  }

  steps = nstep * nblk;

  std::cout << "Total number of steps = " << steps << std::endl;
  std::cout << "Number of blocks = " << nblk << std::endl;
  std::cout << "Number of steps in one block = " << nstep << std::endl << std::endl;
  std::cout << "Step lenght = " << L << std::endl;
  std::cout << "Initial position (GS) = (" << r[0] << " " << r[1] << " " << r[2] << ")" << std::endl;
  std::cout << "Initial position (ExSt) = (" << re[0] << " " << re[1] << " " << re[2] << ")" << std::endl;
  std::cout << std::endl;

  Delete_old_files();
  ReadInput.close();
}

void Delete_old_files(){
    pattern = "05.1_" + prob_str + "_" + std::to_string(r[0]).substr(0, 4) + "_" + std::to_string(r[1]).substr(0, 4) + "_" + std::to_string(r[2]).substr(0, 4);
    std::filesystem::path directory_path = std::string(ROOT_PATH) + "/Data";
    std::cout << "This program will delete files in Data with pattern '" + pattern + "'. Press <enter> to confirm." << std::endl;
    std::cin.ignore();
    for (auto& file : std::filesystem::directory_iterator(directory_path)) {
            if (std::filesystem::is_regular_file(file) && file.path().filename().string().find(pattern) != std::string::npos) {
                    std::filesystem::remove(file.path());
                    std::cout << "Deleted file: " << file.path() << std::endl;
            }
    }
    std::cout << std::endl;    
}

void ResetAllZero(){
    //distances
    d = 0;
    de = 0; 

    //acceptance rate
    attempted_gs = 0;
    attempted_es = 0;
    accepted_gs = 0;
    accepted_es = 0;
}

void SavePos(){

    std::ofstream WritePosGS;
    WritePosGS.open(std::string(ROOT_PATH) + "/Data/" + pattern + "_Ground_State.dat", std::ios::app);
    std::ofstream WritePosES;
    WritePosES.open(std::string(ROOT_PATH) + "/Data/" + pattern + "_Excited_State.dat", std::ios::app);

    for(int i = 0; i < 3; ++i){
        y[i] = r[i];
        ye[i] = re[i];
        WritePosGS << y[i] << " ";
        WritePosES << ye[i] << " ";
    }
    WritePosGS << std::endl;
    WritePosES << std::endl;

    WritePosGS.close();
    WritePosES.close();
}

void ProposeStep(){
    if(!gauss){
        for(int i = 0; i < 3; ++i){
            dr_g[i] = rnd.Rannyu(-L, L);
            dr_e[i] = rnd.Rannyu(-L, L);
        }
    }

    else if(gauss){
        for(int i = 0; i < 3; ++i){
            dr_g[i] = rnd.Gauss(0, L);
            dr_e[i] = rnd.Gauss(0, L);
        }
    }
}

void Move(){

    ProposeStep();

    for(int i = 0; i < 3; ++i){
        x[i] = y[i] + dr_g[i];
        xe[i] = ye[i] + dr_e[i];
    }

    q = prob_ground(x) / prob_ground(y);
    A = std::min(1., q);
    qe = prob_excited(xe) / prob_excited(ye);
    Ae = std::min(1., qe);
    
    if(rnd.Rannyu() < A){
        for(int i = 0; i < 3; ++i){
            r[i] = x[i];
        }
        accepted_gs++;
    }
    attempted_gs++;
    
    if(rnd.Rannyu() < Ae){
        for(int i = 0; i < 3; ++i){
            re[i] = xe[i];
        }
        accepted_es++;
    }
    attempted_es++;
}

void Accumulate(){
    d += sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    de += sqrt(re[0] * re[0] + re[1] * re[1] + re[2] * re[2]);
}

void PrintRate(int i){
    std::cout << "Block number " << i+1 << std::endl;
    std::cout << "                    (Ground)  (Excited)  " << std::endl;
    std::cout << "Acceptance rates:   " << static_cast<double>(accepted_gs) / attempted_gs << "    " << static_cast<double>(accepted_es) / attempted_es << std::endl;
    std::cout << "-----------------------------------" << std::endl;
}

void BlockAverages(){
    d /= nstep;   
    mean_prog_d += d;
    var_prog_d += d * d;

    de /= nstep; 
    mean_prog_de += de;
    var_prog_de += de * de;
}

void SaveDist(int blknum){

    std::ofstream WriteResultGS;
    WriteResultGS.open(std::string(ROOT_PATH) + "/Data/" + pattern + "_Dist_Ground_State.dat", std::ios::app);
    std::ofstream WriteResultES;
    WriteResultES.open(std::string(ROOT_PATH) + "/Data/" + pattern + "_Dist_Excited_State.dat", std::ios::app);

    WriteResultGS << blknum << " " << d << " " << mean_prog_d / (blknum + 1) << " " << Error(mean_prog_d / (blknum + 1), var_prog_d / (blknum + 1), blknum) << " " << std::endl;
    WriteResultES << blknum <<  " " << de << " " << mean_prog_de / (blknum + 1) << " " << Error(mean_prog_de / (blknum + 1), var_prog_de / (blknum + 1), blknum) << " " << std::endl;   

    WriteResultGS.close();
    WriteResultES.close();
}


double prob_ground(double x[3]){
    double d = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    double psi = pow(M_E, -d) / sqrt(M_PI);
    return psi * psi;
}

double prob_excited(double x[3]){
    double d = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    double costheta = x[2]/d;
    double psi = 1. / 8. * sqrt(2. / M_PI) * d * pow(M_E, -d / 2) * costheta;
    return psi * psi;
}

double Error(double sum, double sum2, int iblk){
    if(iblk == 0){
      return 0;
   }
   else{
      return sqrt((sum2 - sum * sum) / iblk);
   }
}

void Random_Start(Random &random_generator){
	
	int seed[4];
	int p1, p2;
	std::ifstream Primes(std::string(ROOT_PATH) + "/random-library/Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else {std::cerr << "PROBLEM: Unable to open Primes" << std::endl;}
	
	Primes.close();
	
	std::ifstream input(std::string(ROOT_PATH) + "/random-library/seed.in");
	std::string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				random_generator.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else {std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;}
}