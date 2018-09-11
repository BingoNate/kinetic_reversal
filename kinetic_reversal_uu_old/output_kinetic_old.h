#ifndef OUTPUT_KINETIC_H
#define OUTPUT_KINETIC_H

#include "NR_kinetic_old.h"
#include "kinetic_parameters_old.h"

void outscreen(int m, double ***ftu, double ***fmtu) {
  int k, i, j;
  double sum = 0;
  double psum = 0;
  double rho = 0;
  double sump1 = 0, sumq1 = 0;
  double sump2 = 0, sumq2 = 0;
  double aver_polar = 0, aver_nematic = 0;
  double **p1, **q1;
  double **p2, **q2;
  p1 = matrix(1, Nx, 1, Ny);
  q1 = matrix(1, Nx, 1, Ny);
  p2 = matrix(1, Nx, 1, Ny);
  q2 = matrix(1, Nx, 1, Ny);

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++) {
      p1[i][j] = 0;
      q1[i][j] = 0;
      p2[i][j] = 0;
      q2[i][j] = 0;

      for (k = 1; k <= M; k++) 
      {
        sum += ftu[k][i][j];
        psum += fmtu[k][i][j];
        p1[i][j] += ftu[k][i][j] * (2 * SQR(ucx[k]) - 1) * dta;
        q1[i][j] += ftu[k][i][j] * (2 * ucx[k] * usy[k]) * dta;
        p2[i][j] += ftu[k][i][j] * ucx[k] * dta;
        q2[i][j] += ftu[k][i][j] * usy[k] * dta;
      }

      p1[i][j] *= 2.0;
      q1[i][j] *= 2.0;
      p2[i][j] *= 2.0;
      q2[i][j] *= 2.0;
      sump1 += p1[i][j];
      sumq1 += q1[i][j];
      sump2 += p2[i][j];
      sumq2 += q2[i][j];
    }

  rho = sum * 2.0 * dta;
  aver_nematic =
      sqrt(SQR(sump1) + SQR(sumq1)) / (2.0 * rho);  // rho = sum * 2.0 * dta.
  aver_polar = sqrt(SQR(sump2) + SQR(sumq2)) / rho;
  printf("T = %d:\n", m);
  printf("ft = %.6f, fm = %.6f\nNematic = %.6f, Polar = %.6f\n\n", sum, psum,
         aver_nematic, aver_polar);  
  free_matrix(p1, 1, Nx, 1, Ny);
  free_matrix(q1, 1, Nx, 1, Ny);
  free_matrix(p2, 1, Nx, 1, Ny);
  free_matrix(q2, 1, Nx, 1, Ny);
}

void output(int m, double ***ftu, double ***fmtu, double ***kru, double *x1,
            double *x2) {
  int i, j, k;
  char subfile[300];
  sprintf_s(subfile, "%s%s", Subdirectory, "Adistribution_");
  PREP(subfile, m, "w")  // distributions on every lattice.

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++) {
      for (j = 1; j <= Ny; j++) {
        fprintf_s(fp, "%d,%d,%d,%e,%e\n", k, i, j, ftu[k][i][j], fmtu[k][i][j]);
      }
    }

  CLSP

      for (i = 1; i <= Nx; i++) for (j = 1; j <= Ny; j++) kru[1][i][j] = 0;

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++)
      for (k = 1; k <= M; k++) {
        kru[1][i][j] +=
            ftu[k][i][j];  // Total distribution on all lattice grid.
      }

  sprintf_s(subfile, "%s%s", Subdirectory, "density_");  // Density.
  PREP(subfile, m, "w")

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++)
      fprintf_s(fp, "%d,%d,%f\n", i, j, kru[1][i][j] * 2.0 * dta);

  CLSP

      for (i = 1; i <= Nx; i++) for (j = 1; j <= Ny; j++) {
    kru[2][i][j] = 0;
    kru[3][i][j] = 0;
    kru[4][i][j] = 0;
  }

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++)
      for (k = 1; k <= M; k++) {
        kru[2][i][j] += ftu[k][i][j] * ucx[k];  // Polar.
        kru[3][i][j] += ftu[k][i][j] * usy[k];
      }

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++) {
      kru[2][i][j] /= kru[1][i][j];
      kru[3][i][j] /= kru[1][i][j];
    }
  /*
  sprintf_s(subfile, "%s%s", Subdirectory, "polar_");
  PREP(subfile, m, "a+")  // Every lattice's polar elements.

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++)
      fprintf_s(fp, "%f,%f,%f,%f,%f,%f\n", x1[i], x2[j],
                x1[i] + kru[2][i][j] * dxx, x2[j] + kru[3][i][j] * dyy,
                kru[2][i][j], kru[3][i][j]);

  CLSP;
*/
  double Sxx = 0;
  double Sxy = 0;
  double Syy = 0;
  double Pox = 0, Poy = 0;

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++) {
      kru[2][i][j] = 0;
      kru[3][i][j] = 0;
      kru[4][i][j] = 0;
      kru[5][i][j] = 0;
      kru[6][i][j] = 0;
      kru[7][i][j] = 0;
      kru[8][i][j] = 0;
    }

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++)
      for (k = 1; k <= M; k++) {
        kru[2][i][j] += ftu[k][i][j] * (2 * ucx[k] * ucx[k] - 1);
        kru[3][i][j] += ftu[k][i][j] * (2 * usy[k] * usy[k] - 1);
        kru[4][i][j] += ftu[k][i][j] * 2 * usy[k] * ucx[k];
        kru[7][i][j] += ftu[k][i][j] * ucx[k];
        kru[8][i][j] += ftu[k][i][j] * usy[k];
      }

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++) {
      kru[2][i][j] /= kru[1][i][j];
      kru[3][i][j] /= kru[1][i][j];
      kru[4][i][j] /= kru[1][i][j];
      kru[7][i][j] /= kru[1][i][j];
      kru[8][i][j] /= kru[1][i][j];
    }

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++) {
      Sxx += kru[2][i][j];
      Syy += kru[3][i][j];
      Sxy += kru[4][i][j];
      Pox += kru[7][i][j];
      Poy += kru[8][i][j];
    }

  Sxx /= (1.0 * Nx * Ny);
  Sxy /= (1.0 * Nx * Ny);
  Syy /= (1.0 * Nx * Ny);
  Pox /= (1.0 * Nx * Ny);
  Poy /= (1.0 * Nx * Ny);
  double orderS = 0, angleS = 0;
  double orderP = 0, angleP = 0;
  orderS = sqrt(Sxy * Sxy - Sxx * Syy);
  orderP = sqrt(SQR(Pox) + SQR(Poy));
  angleS = Sxy / (orderS - Sxx);
  angleS = atan(1 / angleS);
  angleP = atan2(Poy, Pox);
  PREP("GlobalOrder", 0, "a+")  // global order parameter time-series.
  fprintf_s(fp, "%d,%e,%e,%e,%e\n", m, orderS, orderP, angleS, angleP);
  CLSP
      /*
      for (i = 1; i <= Nx; i++) for (j = 1; j <= Ny; j++) {
    kru[5][i][j] =
        sqrt(kru[4][i][j] * kru[4][i][j] -
             kru[2][i][j] * kru[3][i][j]);  // eigenvalue of order parameter
    kru[6][i][j] = kru[4][i][j] /
                   (kru[5][i][j] - kru[2][i][j]);  // angle of order parameter
    kru[6][i][j] = atan(1 / kru[6][i][j]);         // angle of order parameter
  }

  sprintf_s(subfile, "%s%s", Subdirectory, "orders_");
  PREP(subfile, m, "w")  // Order parameters.

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++)
      fprintf_s(fp, "%f,%f,%f,%f\n", x1[i], x2[j], kru[6][i][j],
                kru[5][i][j]);  // x y angle magnitude

  CLSP
*/
      sprintf_s(subfile, "%s%s", Subdirectory, "flow_");  // Particle flow.
  PREP(subfile, m, "w")

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++) {
      double jax, jay;
      double jb, jt;
      int ai, bi, aj, bj;
      ai = i - 1;
      bi = i + 1;
      aj = j - 1;
      bj = j + 1;

      if (ai <= 0) ai += Nx;

      if (bi > Nx) bi -= Nx;

      if (aj <= 0) aj += Ny;

      if (bj > Ny) bj -= Ny;

      jax = 0;
      jay = 0;

      for (k = 1; k <= M; k++) {
        jax += -DFP * ucx[k] * ucx[k] * (ftu[k][bi][j] - ftu[k][ai][j]) /
                   (2 * dxx) -
               DFP * ucx[k] * usy[k] * (ftu[k][i][bj] - ftu[k][i][aj]) /
                   (2 * dyy) -
               DFQ * (ftu[k][bi][j] - ftu[k][ai][j]) / (2 * dxx);
        jay += -DFP * ucx[k] * usy[k] * (ftu[k][bi][j] - ftu[k][ai][j]) /
                   (2 * dxx) -
               DFP * usy[k] * usy[k] * (ftu[k][i][bj] - ftu[k][i][aj]) /
                   (2 * dyy) -
               DFQ * (ftu[k][i][bj] - ftu[k][i][aj]) / (2 * dyy);
      }

      jb = sqrt(jax * jax + jay * jay);
      jt = jay / jb;
      jt = asin(jt);

      if (jax < 0) jt = PI - jt;

      fprintf_s(fp, "%f,%f,%f,%f\n", x1[i], x2[j], jt, jb);
    }

  CLSP
}

#endif  //! OUTPUT_KINETIC_H
