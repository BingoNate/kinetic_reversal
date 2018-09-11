#ifndef KINETIC_PARAMETERS_H
#define KINETIC_PARAMETERS_H

#include "nr3.h"

const int N = 40;   // Square system size.
const int M = 30;   // Splitting the 2*PI angles to M pieces.
const int SST = 1;  // Actually no meaning.
const int NDFR = 0;
const double PI = 3.14159265358979323846;
const double dta = PI / M;  // Small theta pieces.
const double DRV = 0.0;     // No meaning.

int TOTALTIME;  // Total simulation time.
int INTER1;
int INTER2;
int Nx;  // Actually use.
int Ny;

double v0;        // self-propelled velocity for drift terms.
double reversal;  // reversal rate of particle.
double D_para;
double D_perp;
double DNN;  // D_para - D_perp; vary from 0 to 2.8.
double DPP;  // D_para + D_perp; Keep 2.8 fixed.
double DFR;  // Dr: rotaional diffusion. default = 4.0.
double DFQ;  //= 0.50 * (DPP - DNN);// Actually is D_perp.
double DFP;  //= DNN;        // D_para - D_perp.
double dxx;  // Grid size.
double dyy;

long seed = static_cast<long>(time(NULL));

static int set_q;
static double DENS;  // GLobal rescaled density.
static double Nosi;  // defaut = 0.001.
static double dh;    // delta_x = delta_y; grid size default = 3.0.
static double dtt;   // time step. default = 0.01.

static double UXU[M + 1][M + 1];
static double usy[M + 1];  // Store the sin(\theta).
static double ucx[M + 1];
static double rusy[M + 1];
static double rucx[M + 1];

static int LX;  // System size.
static int LY;

#ifdef _MSC_VER
#define Subdirectory "Data\\"
#else
#define Subdirectory "Data//"
#endif  // _MSC_VER

#endif  // !KINETIC_PARAMETERS_H
