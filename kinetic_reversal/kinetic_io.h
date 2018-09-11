#ifndef KINETIC_IO_H
#define KINETIC_IO_H

#include "kinetic_parameters.h"

void make_dir() {
  char subfile[100];
  sprintf_s(subfile, "mkdir %s", Subdirectory);
  system(subfile);
}

void adjust_para() {
  DPP = D_para + D_perp;
  DNN = D_para - D_perp;
  // Set anisotropic diffusion constants. Isotropic is set D_para = D_perp.
  // Best set may be set D_para = v0^2 / (2 * reversal).
  DFQ = D_perp;
  DFP = D_para - D_perp;
  dh = dxx;
  LX = static_cast<int>(dxx * Nx);
  LY = static_cast<int>(dyy * Ny);
  dh = dxx;
  dtt *= 3.0 / 4.0;  // Only 1 equation. Rescale delta t.
}

void read_para() {
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
  CLSP

      set_q = 0;
  adjust_para();
}

void initial_set(int set_q, vector<double>& x1, vector<double>& x2,
                 Mat3DDoub& ftu, Mat3DDoub& pftu, Mat3DDoub& fmtu,
                 Mat3DDoub& pfmtu) {
  Ran myran(seed);
  if (set_q == 0) {
    for (int i = 0; i < M; i++) {
      double dtheta = i * dta - 0.5 * dta;
      usy[i] = sin(dtheta);
      ucx[i] = cos(dtheta);
    }
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < M; j++) {
        UXU[i][j] = fabs(usy[i] * ucx[j] - ucx[i] * usy[j]);
      }
    }
    for (int i = 0; i < Nx; i++) {
      x1[i] = dxx * i;
    }
    for (int i = 0; i < Ny; i++) {
      x2[i] = dxx * i;
    }

    for (int k = 0; k < M; k++) {
      for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
          //*exp(-0.01*(pow(x1[i]-x1[N/2],2)+pow(x2[j]-x2[N/2],2)))*exp(-(k-M/2)*(k-M/2)*0.01)
          ftu[k][i][j] = DENS * (PI * 3.0 / 2.0 + Nosi * (myran.doub() - 0.5)) /
                         (2.0 * PI);  //+0.5*exp(-pow(x2[j]-x2[N/2],2));
          fmtu[k][i][j] = ftu[k][i][j];
        }
      }
    }
  } else {
    for (int k = 0; k < M; k++) {
      for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
          //*exp(-0.01*(pow(x1[i]-x1[N/2],2)+pow(x2[j]-x2[N/2],2)))*exp(-(k-M/2)*(k-M/2)*0.01)
          ftu[k][i][j] = DENS * (PI * 3.0 / 2.0 + Nosi * (myran.doub() - 0.5)) /
                         (2.0 * PI);  //+0.5*exp(-pow(x2[j]-x2[N/2],2));
          fmtu[k][i][j] = ftu[k][i][j];
        }
      }
    }
  }
  // Transfer.
  for (int k = 0; k < M; k++)
    for (int i = 0; i < Nx; i++)
      for (int j = 0; j < Ny; j++) {
        pftu[k][i][j] = ftu[k][i][j];
        pfmtu[k][i][j] = fmtu[k][i][j];
      }
}

void outscreen(int mt, const Mat3DDoub& ftu, const Mat3DDoub& fmtu) {
  double sum = 0;
  double pm = 0;
  printf_s("T = %d\n", mt);

  for (int k = 0; k < M; k++)
    for (int i = 0; i < Nx; i++)
      for (int j = 0; j < Ny; j++) {
        sum += ftu[k][i][j];
        pm += fmtu[k][i][j];
      }

  printf_s("Distribution: ft = %f, fm = %f\n", sum, pm);  // Number conservation
}

void output(int m_step, const Mat3DDoub& ftu, const Mat3DDoub& fmtu,
            Mat3DDoub& kru, const vector<double>& x1,
            const vector<double>& x2) {
  int i, j, k;

  if (m_step % INTER1 == 0) {
    char subfile[300];
    sprintf_s(subfile, "%s%s", Subdirectory, "a_");
    PREP(subfile, m_step, "w")  // distribution.

    for (k = 0; k < M; k++)
      for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
          fprintf_s(fp, "%d,%d,%d,%e,%e\n", k, i, j, ftu[k][i][j], fmtu[k][i][j]);
        }
      }

    CLSP;
    for (i = 0; i < Nx; i++) {
      for (j = 0; j < Ny; j++) {
        kru[1][i][j] = 0;
      }
    }

    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++)
        for (k = 0; k < M; k++) {
          kru[1][i][j] += ftu[k][i][j];  // Total distribution on all pi angles.
        }

    sprintf_s(subfile, "%s%s", Subdirectory, "d_");  // Density.
    PREP(subfile, m_step, "w")

    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++)
        fprintf_s(fp, "%d,%d,%f\n", i, j, kru[1][i][j] * 2.0 * PI / M);

    CLSP;

    for (i = 0; i < Nx; i++) {
      for (j = 0; j < Ny; j++) {
        kru[2][i][j] = 0;
        kru[3][i][j] = 0;
        kru[4][i][j] = 0;
      }
    }
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++)
        for (k = 0; k < M; k++) {
          kru[2][i][j] += ftu[k][i][j] * ucx[k];
          kru[3][i][j] += ftu[k][i][j] * usy[k];
        }

    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        kru[2][i][j] /= kru[1][i][j];
        kru[3][i][j] /= kru[1][i][j];
      }

    double Sxx = 0, Sxy = 0;
    double Syy = 0;

    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        kru[2][i][j] = 0;
        kru[3][i][j] = 0;
        kru[4][i][j] = 0;
        kru[5][i][j] = 0;
        kru[6][i][j] = 0;
      }

    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++)
        for (k = 0; k < M; k++) {
          kru[2][i][j] += ftu[k][i][j] * (2 * ucx[k] * ucx[k] - 1);
          kru[3][i][j] += ftu[k][i][j] * (2 * usy[k] * usy[k] - 1);
          kru[4][i][j] += ftu[k][i][j] * 2 * usy[k] * ucx[k];
        }

    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        kru[2][i][j] /= kru[1][i][j];
        kru[3][i][j] /= kru[1][i][j];
        kru[4][i][j] /= kru[1][i][j];
      }

    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        Sxx += kru[2][i][j];
        Syy += kru[3][i][j];
        Sxy += kru[4][i][j];
      }

    Sxx /= (1.0 * Nx * Ny);
    Sxy /= (1.0 * Nx * Ny);
    Syy /= (1.0 * Nx * Ny);

    double orderS = 0, angleS = 0;
    orderS = sqrt(Sxy * Sxy - Sxx * Syy);
    angleS = Sxy / (orderS - Sxx);
    angleS = atan(1.0 / angleS);
    printf_s("Nematic Order = %f, Angle = %f\n", orderS, angleS);
    PREP("order_tran", 0, "a+")
    fprintf_s(fp, "%d,%e,%e\n", m_step, orderS, angleS);
    CLSP
  }
}

#endif  // !KINETIC_IO_H
