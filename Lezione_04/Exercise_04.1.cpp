#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include "Library_04.h"
#include <version_config.h>

std::string random_lib_path = std::string(ROOT_PATH) + "/random-library/";
std::string input_path, output_path, eq_str, input_form_eq_path;


int main(int argc, char* argv[]){

	//Usage and phase choosing
	if(argc == 3){
		phase = static_cast<std::string>(argv[1]);
		eq = static_cast<std::string>(argv[2]);
		if ((phase != "solid" && phase != "liquid" && phase != "gas") && (eq != "false" && eq != "true")){ 
    	std::cerr << "Wrong usage. Usage: './Exercise_04.1 <phase> <equilibration>' with <phase> = {solid, liquid, gas} and <equilibration> = {true, false}" << std::endl;
			return -1;
  	}
	}
	else{
		std::cerr << "Wrong usage. Usage: './Exercise_04.1 <phase> <equilibration>' with <phase> = {solid, liquid, gas} and <equilibration> = {true, false}" << std::endl;
		return -1;
	}
  
	if(eq == "true"){
		std::cout << "You choose EQUILIBRATION mode" << std::endl;
		input_path = "input/eq/" + phase + "/input." + phase;
		output_path = "input/eq/" + phase + "/";
		eq_str = "eq_";
	}
	else if (eq == "false"){
		std::cout << "You choose SIMULATION mode" << std::endl;
		input_form_eq_path = "input/eq/" + phase + "/";
		input_path = "input/input." + phase;
		output_path = "";
		eq_str = "";
	}
	else{}

	//delete old files
	std::string pattern = "04.1_" + eq_str + phase;
  std::filesystem::path directory_path = std::string(ROOT_PATH) + "/Data";
	std::string confirm = "";
	std::cout << "This program will delete files in Data with pattern '" + pattern + "'. Press <enter> to confirm." << std::endl;
	std::cin.ignore();
	for (auto& file : std::filesystem::directory_iterator(directory_path)) {
			if (std::filesystem::is_regular_file(file) && file.path().filename().string().find(pattern) != std::string::npos) {
					std::filesystem::remove(file.path());
					std::cout << "Deleted file: " << file.path() << std::endl;
			}
	}
	std::cout << std::endl;

	//Inizialization
  Input(); 
	std::cout << std::endl;

	//------------------------Progress bar declaration--------------------------
	int stp2 = nstep/10;
	int s = 0;
	//--------------------------------------------------------------------------

	//Start simulation
  int nconf = 1;
  for(int iblk = 1; iblk <= nblk; iblk++){
		s = 0;

    Reset(iblk);   //Reset block averages

    for(int istep = 1; istep <= nstep; istep++){

      Move();
      Measure();
      Accumulate(); //Update block averages

			//Write actual configuration in XYZ format
      /*if(istep%10 == 0){
        ConfXYZ(nconf);
        nconf += 1;
      }*/

			//----------------------------Progress Bar------------------------------
			if (istep % stp2 == 0){
				Progress_bar(s, iblk, istep);
			}
			//----------------------------------------------------------------------
    }
		std::cout << 	Green("Block number " + std::to_string(iblk) + " completed successfully. Progress: ▪▪▪▪▪▪▪▪▪▪ 100%") << std::endl;
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  rnd.SaveSeed(output_path + "seed.out"); //Save seed

  return 0;
}





void Input(void){

  std::ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  std::cout << "Classic Lennard-Jones fluid        " << std::endl;
  std::cout << "MD(NVE) / MC(NVT) simulation       " << std::endl << std::endl;
  std::cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << std::endl << std::endl;
  std::cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << std::endl << std::endl;
  std::cout << "The program uses Lennard-Jones units " << std::endl;

//Read seed for random numbers
  int p1, p2;
  Primes.open(random_lib_path + "Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

	//Read input informations

  ReadInput.open(input_path);
  ReadInput >> iNVET;
  ReadInput >> restart;

  if(restart) Seed.open(input_form_eq_path + "seed.out");
  else Seed.open(random_lib_path + "seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  Seed.close();

  ReadInput >> temp;
  beta = 1.0 / temp;
  std::cout << "Temperature = " << temp << std::endl;

  ReadInput >> npart;
  std::cout << "Number of particles = " << npart << std::endl;

  ReadInput >> rho;
  std::cout << "Density of particles = " << rho << std::endl;
  vol = (double)npart/rho;
  box = pow(vol, 1.0 / 3.0);
  std::cout << "Volume of the simulation box = " << vol << std::endl;
  std::cout << "Edge of the simulation box = " << box << std::endl;

  ReadInput >> rcut;
  std::cout << "Cutoff of the interatomic potential = " << rcut << std::endl << std::endl;
    
  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  std::cout << "The program perform Metropolis moves with uniform translations" << std::endl;
  std::cout << "Moves parameter = " << delta << std::endl;
  std::cout << "Number of blocks = " << nblk << std::endl;
  std::cout << "Number of steps in one block = " << nstep << std::endl << std::endl;
  ReadInput.close();

	//Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  iw = 4; //Total pressure
  n_props = 5; //Number of observables

	//Read initial configuration
  std::cout << "Read initial configuration" << std::endl << std::endl;
  if(restart){
    ReadConf.open(input_form_eq_path + "config.out", std::ios::out);
    ReadVelocity.open(input_form_eq_path + "velocity.out", std::ios::out);
    for (int i = 0; i < npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else{
    ReadConf.open("input/eq/config.in");
    std::cout << "Prepare velocities with center of mass velocity equal to zero " << std::endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < npart; ++i)
    {
      vx[i] = rnd.Gauss(0., sqrt(temp));
      vy[i] = rnd.Gauss(0., sqrt(temp));
      vz[i] = rnd.Gauss(0., sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim = 0; idim < 3; ++idim) sumv[idim] /= (double)npart; //media sulla singola direzione
    double sumv2 = 0.0, fs;
    for (int i = 0; i < npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    } 
    //deviazione dalla media
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    std::cout << "velocity scale factor: " << fs << std::endl << std::endl;
    for (int i = 0; i < npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i = 0; i < npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

  for (int i = 0; i < npart; ++i)
  {
    if(iNVET)
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
//Evaluate properties of the initial configuration
  Measure();

//Print initial values for measured properties
  std::cout << "Initial potential energy = " << walker[iv] / (double)npart << std::endl;
  std::cout << "Initial temperature      = " << walker[it] << std::endl;
  std::cout << "Initial kinetic energy   = " << walker[ik] / (double)npart << std::endl;
  std::cout << "Initial pressure         = " << walker[iw] << std::endl;
  std::cout << "Initial total energy     = " << walker[ie] / (double)npart << std::endl;

  return;
}


void Move()
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.Rannyu()*npart);

    //Old
      energy_old = Boltzmann(x[o],y[o],z[o],o);

    //New
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o], y[o], z[o],o);

    //Metropolis test
      p = exp(beta * (energy_old - energy_new));
      if(p >= rnd.Rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i = 0; i < npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i = 0; i < npart; ++i){ //Verlet integration scheme

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i]) / (2.0 * delta);
      vy[i] = Pbc(ynew - yold[i]) / (2.0 * delta);
      vz[i] = Pbc(znew - zold[i]) / (2.0 * delta);

      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;
      attempted = attempted + 1.0;
    }
  }
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene = 0.0;
  double dx, dy, dz, dr;

  for (int i = 0; i < npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx * dx + dy * dy + dz * dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0 * ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f = 0.0;
  double dvec[3], dr;

  for (int i = 0; i < npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure() //Properties measurement
{
  double v = 0.0, kin = 0.0, w = 0.0;
  double vij, w_ij;
  double dx, dy, dz, dr;

//cycle over pairs of particles
  for (int i = 0; i < npart - 1; ++i)
  {
    for (int j = i + 1; j < npart; ++j)
    {
// distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx * dx + dy * dy + dz * dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        w_ij = 1./pow(dr, 12) - 0.5/pow(dr, 6.0);
        v += vij;
        w += w_ij;
      }
    }          
  }

  for (int i = 0; i < npart; ++i) kin += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

  walker[iv] = 4.0 * v; // Potential energy
  walker[ik] = kin; // Kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
  walker[ie] = 4.0 * v + kin;  // Total energy;
  walker[iw] = 48.0/(3.0 * vol) * w + rho * walker[it];  //Pressure

  return;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i = 0; i < n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i = 0; i < n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i = 0; i < n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   std::ofstream Epot, Ekin, Etot, Temp, Press;
   const int wd = 12;
    
    std::cout << "Block number " << iblk << std::endl;
    std::cout << "Acceptance rate " << accepted/attempted << std::endl << std::endl;
    
    Epot.open(std::string(ROOT_PATH) + "/Data/04.1_" + eq_str + phase + "_epot.dat", std::ios::app);
    Ekin.open(std::string(ROOT_PATH) + "/Data/04.1_" + eq_str + phase + "_ekin.dat", std::ios::app);
    Temp.open(std::string(ROOT_PATH) + "/Data/04.1_" + eq_str + phase + "_temp.dat", std::ios::app);
    Etot.open(std::string(ROOT_PATH) + "/Data/04.1_" + eq_str + phase + "_etot.dat", std::ios::app);
    Press.open(std::string(ROOT_PATH) + "/Data/04.1_" + eq_str + phase + "_press.dat", std::ios::app);
    
    stima_pot = blk_av[iv] / blk_norm / (double)npart; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot * stima_pot;
    err_pot = Error(glob_av[iv], glob_av2[iv], iblk);
    
    stima_kin = blk_av[ik] / blk_norm / (double)npart; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin * stima_kin;
    err_kin = Error(glob_av[ik], glob_av2[ik], iblk);

    stima_etot = blk_av[ie] / blk_norm / (double)npart; //Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot * stima_etot;
    err_etot = Error(glob_av[ie], glob_av2[ie], iblk);

    stima_temp = blk_av[it] / blk_norm; //Temperature
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp * stima_temp;
    err_temp = Error(glob_av[it], glob_av2[it], iblk);
    
    stima_pres = blk_av[iw] / blk_norm; //Press
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres * stima_pres;
    err_press = Error(glob_av[iw], glob_av2[iw], iblk);

    //Potential energy per particle
    Epot << std::setw(wd) << iblk <<  std::setw(wd) << stima_pot << std::setw(wd) << glob_av[iv]/(double)iblk << std::setw(wd) << err_pot << std::endl;
    //Kinetic energy
    Ekin << std::setw(wd) << iblk <<  std::setw(wd) << stima_kin << std::setw(wd) << glob_av[ik]/(double)iblk << std::setw(wd) << err_kin << std::endl;
    //Total energy
    Etot << std::setw(wd) << iblk <<  std::setw(wd) << stima_etot << std::setw(wd) << glob_av[ie]/(double)iblk << std::setw(wd) << err_etot << std::endl;
    //Temperature
    Temp << std::setw(wd) << iblk <<  std::setw(wd) << stima_temp << std::setw(wd) << glob_av[it]/(double)iblk << std::setw(wd) << err_temp << std::endl;
    //Pressure
    Press << std::setw(wd) << iblk <<  std::setw(wd) << stima_pres << std::setw(wd) << glob_av[iw]/(double)iblk << std::setw(wd) << err_press << std::endl;

    std::cout << "----------------------------" << std::endl << std::endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Press.close();
}


void ConfFinal(void)
{
  std::ofstream WriteConf, WriteVelocity, WriteSeed;

  std::cout << "Print final configuration to file config.out" << std::endl << std::endl;
  WriteConf.open(output_path + "config.out");
  WriteVelocity.open(output_path + "velocity.out");
  for (int i = 0; i < npart; ++i)
  {
    WriteConf << x[i] / box << "   " <<  y[i] / box << "   " << z[i] / box << std::endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << std::endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  std::ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + std::to_string(nconf) + ".xyz");
  WriteXYZ << npart << std::endl;
  WriteXYZ << "This is only a comment!" << std::endl;
  for (int i = 0; i < npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << std::endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum / (double)iblk, 2)) / (double)iblk);
}