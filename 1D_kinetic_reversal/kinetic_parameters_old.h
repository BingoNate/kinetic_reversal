#ifndef KINETIC_PARAMETERS_H
#define KINETIC_PARAMETERS_H

#include "NR_kinetic_old.h"

const int M = 30;   // Splitting the 2*PI angles to M pieces.
const int SST = 1;  // Actually no meaning.
const int NDFR = 0;
const int GP = 0;        // sample GP + 1 points.
const double DRV = 0.0;  // No meaning.
const double PI = 3.14159265358979323846;
const double dta = PI / M;  // Small theta pieces.

/* Need to input. */
int nowtime;
int TOTALTIME;  // Total simulation time.
int INTER1;
int INTER2;
int Nx;  // Actually use.
int Ny;
int N;      // Square system size. Equal to Nx.
double LX;  // Real system size.
double LY;

double DENS;  // GLobal rescaled density.
double Nosi;  // defaut = 0.001.
double dh;    // delta_x = delta_y; grid size default = 3.0.
double dtt;   // time step. default = 0.01.

double v0;        // self-propelled velocity for drift terms.
double reversal;  // reversal rate of particle.
double D_para;    // Set to small, i.e. 0.10.
double D_perp;
double DNN;  // D_para - D_perp; vary from 0 to 2.8.
double DPP;  // D_para + D_perp; Keep 2.8 fixed.
double DFR;  // Dr: rotaional diffusion. default = 2.0.
double DFQ;  //= 0.50 * (DPP - DNN);// Actually is D_perp.
double DFP;  //= DNN;        // D_para - D_perp.
double dxx;  // Grid size.
double dyy;
/*********************/
double UXU[M + 1][M + 1];
double usy[M + 1];  // Store the sin(\theta).
double ucx[M + 1];
double rusy[M + 1];
double rucx[M + 1];
long static RSD[1] = {static_cast<long>(-time(NULL))};

#ifdef _MSC_VER
#define Subdirectory "Data\\"
#else
#define Subdirectory "Data//"
#endif  // _MSC_VER

#endif  //! KINETIC_PARAMETERS_H
