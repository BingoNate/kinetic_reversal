#ifndef INITIAL_H
#define INITIAL_H

#include "NR_kinetic_old.h"
#include "kinetic_parameters_old.h"

const int gap_ini = Ny / 4;

#define PREP_check(Fpre, Fnum, R_W)              \
  {                                              \
    FILE *fp;                                    \
    errno_t err;                                 \
    {                                            \
      char PRT[100];                             \
      sprintf_s(PRT, "%s%d.dat", Fpre, Fnum);    \
      if ((err = fopen_s(&fp, PRT, R_W)) != 0) { \
        nowtime = Fnum;                          \
        yesfile++;                               \
        { fclose(fp); }                          \
      }                                          \
    }                                            \
  }

void make_subdirectory() {
  char subfile[300];
  sprintf_s(subfile, "mkdir %s", Subdirectory);
  system(subfile);
}

void AjustPara() {
  dtt *= 1.0;  // Rescale delta t.
  dyy = 0;
  N = Nx;
  LX = dxx * Nx;
  LY = dyy * Ny;
  dh = dxx;
  DPP = D_para + D_perp;  // Set anisotropic diffusion constants. Isotropic is
                          // set D_para = D_perp.
  DNN = D_para - D_perp;  // Best set may be set D_para = v0^2 / (2 * reversal).
  DFQ = D_perp;
  DFP = DNN;
}

void ReadParameter() {
  PREP("aParameter", 0, "r")
  fscanf_s(fp, "%d=Total_time(TOTALTIME)", &TOTALTIME);
  fscanf_s(fp, "%d=INTER1(INTER1)", &INTER1);
  fscanf_s(fp, "%d=INTER2(INTER2)", &INTER2);
  fscanf_s(fp, "%d=System_size(Nx)", &Nx);
  fscanf_s(fp, "%d=System_size(Ny)", &Ny);
  fscanf_s(fp, "%lf=Density(DENS)", &DENS);
  fscanf_s(fp, "%lf=Nosi(Nosi)", &Nosi);
  fscanf_s(fp, "%lf=Velocity(v0)", &v0);
  fscanf_s(fp, "%lf=Reversal(reversal)", &reversal);
  fscanf_s(fp, "%lf=D_para(D_para)", &D_para);
  fscanf_s(fp, "%lf=D_perp(D_perp)", &D_perp);
  fscanf_s(fp, "%lf=DFR(DFR)", &DFR);
  fscanf_s(fp, "%lf=Dxx(dxx)", &dxx);
  fscanf_s(fp, "%lf=Dyy(dyy)", &dyy);
  fscanf_s(fp, "%lf=Dtt(dtt)", &dtt);
  CLSP;
  AjustPara();
}

void ini_lattice(double *x1, double *x2) {
  int i;

  for (i = 1; i <= Nx; ++i) {
    x1[i] = dxx * i;
  }

  for (i = 1; i <= Ny; ++i) {
    x2[i] = dyy * i;
  }
}

void ini_potential() {
  int i, j;

  for (i = 1; i <= M; ++i) {
    usy[i] = sin(i * dta - 0.5 * dta);
    ucx[i] = cos(i * dta - 0.5 * dta);
  }

  for (i = 1; i <= M; ++i) {
    for (j = 1; j <= M; ++j) {
      // -|uu|:
      UXU[i][j] = 1 - fabs(ucx[i] * ucx[j] + usy[i] * usy[j]);
      // -|uxu|:
      // UXU[i][j] = fabs(usy[i] * ucx[j] - ucx[i] * usy[j]);
    }
  }
}

void ini_zero_distri(double ***ftu, double ***fmtu) {
  int k, i, j;

  /* Initial two point distributions: f_T and f_M */
  for (k = 1; k <= M; ++k) {
    for (i = 1; i <= Nx; ++i) {
      for (j = 1; j <= Ny; ++j) {
        ftu[k][i][j] = 0;
        fmtu[k][i][j] = 0;
      }
    }
  }
}

void check_exist_file(double ***ftu, double ***fmtu) {
  int k, i, j;
  int yesfile = 0;
  int a, b, c;
  char subfile[300];
  sprintf_s(subfile, "%s%s", Subdirectory, "Adistribution_");

  for (i = 1; i <= TOTALTIME; i += INTER1) {
    PREP_check(subfile, i, "r")
  }

  if (yesfile != 0) {
    sprintf_s(subfile, "%s%s", Subdirectory, "Adistribution_");
    PREP(subfile, nowtime, "r")
    printf("Reading distribution: T = %d", nowtime);

    for (k = 1; k <= M; ++k)
      for (i = 1; i <= Nx; ++i)
        for (j = 1; j <= Ny; ++j) {
          fscanf_s(fp, "%d,%d,%d,%lf,%lf", &a, &b, &c, &ftu[k][i][j],
                   &fmtu[k][i][j]);
        }

    CLSP
  }

  // stay above initial condition
  ftu[1][1][SST] += 0.000;
  fmtu[1][1][SST] += 0.000;
}

void transfer_tmp_distri(double ***ftu, double ***fmtu, double ***pftu,
                         double ***pfmtu) {
  int k, i, j;

  for (k = 1; k <= M; ++k)
    for (i = 1; i <= Nx; ++i)
      for (j = 1; j <= Ny; ++j) {
        pftu[k][i][j] = ftu[k][i][j];    // transfer value
        pfmtu[k][i][j] = fmtu[k][i][j];  // transfer value
      }
}

/*
  Precondition: band is along x-axis. (Smallest x = 3.)
  We just need to read 1-line distribution vertical to x.
  Then we copy the total system distribution with it.
  And we can add some small noise to the distributions.
*/
void read_band(double ***ftu, double ***fmtu) {
  /*It is a distribution file.*/
  int k, i, j;
  int a, b, c;
  int shfty = 5;
  double ***tu, ***tmu;
  tu = f3tensor(1, M, 1, Nx, 1, Ny);
  tmu = f3tensor(1, M, 1, Nx, 1, Ny);
  PREP("Band", 0, "r")

  for (k = 1; k <= M; ++k) {
    for (i = 1; i <= 10; ++i) {
      for (j = 1; j <= Ny; ++j) {
        fscanf_s(fp, "%d,%d,%d,%lf,%lf", &a, &b, &c, &tu[k][i][j],
                 &tmu[k][i][j]);
      }
    }
  }

  CLSP;

  /*Set distributions with band state along y-axis.*/
  for (k = 1; k <= M; ++k) {
    for (i = 1; i <= Nx; ++i) {
      for (j = 1; j <= Ny; ++j) {
        ftu[k][i][j] = tu[k][1][j];
        fmtu[k][i][j] = tmu[k][1][j];
      }
    }
  }

  /*Add some small noise to distributions.*/
  for (k = 1; k <= M; ++k) {
    for (i = 1; i <= Nx; ++i) {
      for (j = 1; j <= Ny; ++j) {
        ftu[k][i][j] *= (1.0 + Nosi * (ran1(RSD) - 0.50));
        fmtu[k][i][j] *= (1.0 + Nosi * (ran1(RSD) - 0.50));
      }
    }
  }

  free_f3tensor(tu, 1, M, 1, 10, 1, Ny);
  free_f3tensor(tmu, 1, M, 1, 10, 1, Ny);
}

/*
  Set ordered system with orientation along with x-axis.
  Initialize to generate band state.
*/
void pre_set_band(double ***ftu, double ***fmtu) {
  int k, i, j;

  for (k = 1; k <= M; ++k) {
    for (i = 1; i <= Nx; ++i) {
      for (j = 1; j <= Ny; ++j) {
        if (i >= 40 && i < 60) {
          if (M / 2 == k) {
            ftu[k][i][j] = DENS * (PI * 3.0 / 2.0) / (2.0 * PI) * M * Nx / 20 +
                           Nosi * (ran1(RSD) - 0.5);
          } else {
            ftu[k][i][j] = 0;
          }
        } else {
          ftu[k][i][j] = 0;
        }
        fmtu[k][i][j] = 0.50 * ftu[k][i][j];
      }
    }
  }
}

/**************************************Different Initial
 * Conditions************************************************************/
void initial_continue(double *x1, double *x2, double ***ftu, double ***pftu,
                      double ***fmtu, double ***pfmtu) {
  /*Initial system real size.*/
  ini_lattice(x1, x2);
  /*Initial system potential kernel.*/
  ini_potential();
  /*Set distributions to zero.*/
  ini_zero_distri(ftu, fmtu);
  /*Check files in Data and set current time-step distributions. Reset the
   * time.*/
  check_exist_file(ftu, fmtu);
  /*Transfer temperory distributions.*/
  transfer_tmp_distri(ftu, fmtu, pftu, pfmtu);
}

void initial_read_band(double *x1, double *x2, double ***ftu, double ***pftu,
                       double ***fmtu, double ***pfmtu) {
  /*Initial system real size.*/
  ini_lattice(x1, x2);
  /*Initial system potential kernel.*/
  ini_potential();
  /*Set distributions to zero.*/
  ini_zero_distri(ftu, fmtu);
  /*Check files in Data and set current time-step distributions. Reset the
   * time.*/
  read_band(ftu, fmtu);
  /*Transfer temperory distributions.*/
  transfer_tmp_distri(ftu, fmtu, pftu, pfmtu);
}

void initial(double *x1, double *x2, double ***ftu, double ***pftu,
             double ***fmtu, double ***pfmtu) {
  /*Initial system real size.*/
  ini_lattice(x1, x2);
  /*Initial system potential kernel.*/
  ini_potential();
  /*Set distributions to zero.*/
  ini_zero_distri(ftu, fmtu);
  /*Check files in Data and set current time-step distributions. Reset the
   * time.*/
  pre_set_band(ftu, fmtu);
  /*Transfer temperory distributions.*/
  transfer_tmp_distri(ftu, fmtu, pftu, pfmtu);
}

#endif  // INITIAL_H
