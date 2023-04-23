#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include "random.h"
#include <version_config.h>

Random rnd;  

int nstep = 1e4;                      // numero di step per blocco
int nblk = 100;                   // numero di blocchi
int steps = nstep * nblk;          // numero di step totale
double L = 1.5;                   // lato del cubo visitabile / sigma per caso gaussiano
bool gauss;
int np_rend = 10000;
//int mod = (int)steps/np_rend;
int mod = 1;
int contastep = 0;

std::string prob_str;

// variabili di appoggio per Ground State
double d = 0;           // distanza dal centro
double var_prog_d = 0;  // appoggio per errori statistici
double mean_prog_d = 0;

double r[3]; // posizione iniziale
double y[3];                 // posizione prima del passo
double x[3];                 // posizione dopo il passo

double dr_g[3];

double q, A;
int attempted_gs = 0, accepted_gs = 0;

// variabili di appoggio per Excited State
double de = 0;          // distanza dal centro
double var_prog_de = 0; // appoggio per errori statistici
double mean_prog_de = 0;

double re[3]; // posizione iniziale
double ye[3];                 // posizione prima del passo
double xe[3];                 // posizione dopo il passo

double dr_e[3];

double qe, Ae;
int attempted_es = 0, accepted_es = 0;

//functions
void Input();
double prob_gs(double*);
double prob_exc(double*);
double min(double, double);
void ResetAll();
void SavePos();
void Move();
void Accumulate();
void PrintAccRate(int blknum);
void BlockAverages();
void SaveDist(int blknum);

void Random_Start(Random &random_generator);
double Error(double sum, double sum2, int iblk);

std::string perc = "▪▪▪▪▪▪▪▪▪▪";

std::string Green(std::string green){
	return std::string("\033[1;32m") + green + "\033[0m";
}

std::string Gray(std::string gray){
	return std::string("\033[1;90m") + gray + "\033[0m";
}

std::string Red(std::string red){
	return std::string("\033[1;31m") + red + "\033[0m";
}

std::string pattern;