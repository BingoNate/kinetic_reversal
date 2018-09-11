#include "NR_kinetic_old.h"
#include "initial.h"
#include "kinetic_parameters_old.h"
#include "output_kinetic_old.h"

int int_kp_simple(double ***ftu, double ***kru, double *x1, double *x2) {
  int i, j;
  int k1, k2;

  for (k1 = 1; k1 <= M; k1++) {
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        kru[k1][i][j] = 0;
      }
  }

  for (k1 = 1; k1 <= M; k1++) {
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        for (k2 = 1; k2 <= M; k2++) {
          if ((k2 == k1) || fabs(k2 - k1) == M) {
            continue;
          }
          kru[k1][i][j] += UXU[k1][k2] * ftu[k2][i][j];
        }
      }
  }

  for (k1 = 1; k1 <= M; k1++) {
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        kru[k1][i][j] = kru[k1][i][j] * dta * 2;
      }
  }

  return 1;
}

/*
  This routine caculates the differential in angular space.
*/
int differ1(double ***ftu, double ***pftu, double ***fmtu, double ***pfmtu,
            double ***kru, double ***nftu) {
  int i, j, k;
  int ai, bi, aj, bj, ak, bk;
  double maxd = DFR;
  double ***a, ***b, ***c;
  double **alpha, **beta;
  double *rr, *rx;
  a = f3tensor(1, Nx, 1, Ny, 1, M);
  b = f3tensor(1, Nx, 1, Ny, 1, M);
  c = f3tensor(1, Nx, 1, Ny, 1, M);
  alpha = matrix(1, Nx, 1, Ny);
  beta = matrix(1, Nx, 1, Ny);
  rr = vector(1, M);
  rx = vector(1, M);

  /******************************************************************************************************/
  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++)
      for (k = 1; k <= M; k++) {
        ak = k - 1;
        bk = k + 1;

        if (ak < 1) ak += M;

        if (bk > M) bk -= M;

        b[i][j][k] = 2 * maxd * dtt / (3 * dta * dta) -
                     DFR * (kru[bk][i][j] - 2 * kru[k][i][j] + kru[ak][i][j]) *
                         dtt / (6 * dta * dta) +
                     1;
      }

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++)
      for (k = 2; k <= M; k++) {
        ak = k - 1;
        bk = k + 1;

        if (ak < 1) ak += M;

        if (bk > M) bk -= M;

        a[i][j][k] =
            -maxd * dtt / (3 * dta * dta) -
            DFR * (kru[ak][i][j] - kru[k][i][j]) * dtt / (6 * dta * dta);  //?
      }

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++)
      for (k = 1; k <= M - 1; k++) {
        ak = k - 1;
        bk = k + 1;

        if (ak < 1) ak += M;

        if (bk > M) bk -= M;

        c[i][j][k] =
            -maxd * dtt / (3 * dta * dta) -
            DFR * (kru[bk][i][j] - kru[k][i][j]) * dtt / (6 * dta * dta);
      }

  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++) {
      k = M;
      ak = k - 1;
      bk = k + 1;

      if (ak < 1) ak += M;

      if (bk > M) bk -= M;

      alpha[i][j] =
          -maxd * dtt / (3 * dta * dta) -
          DFR * (kru[bk][i][j] - kru[k][i][j]) * dtt / (6 * dta * dta);
      k = 1;
      ak = k - 1;
      bk = k + 1;

      if (ak < 1) ak += M;

      if (bk > M) bk -= M;

      beta[i][j] = -maxd * dtt / (3 * dta * dta) -
                   DFR * (kru[ak][i][j] - kru[k][i][j]) * dtt / (6 * dta * dta);
    }

  /******************************************************************************************************/

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = 0;
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;

        if (ai <= 0) ai += Nx;

        if (bi > Nx) bi -= Nx;

        if (aj <= 0) aj += Ny;

        if (bj > Ny) bj -= Ny;

        if (ak < 1) ak += M;

        if (bk > M) bk -= M;

        nftu[k][i][j] += (DFP *
                          (ucx[k] * ucx[k] *
                               ((ftu[k][bi][j] - ftu[k][i][j]) -
                                (ftu[k][i][j] - ftu[k][ai][j])) +
                           usy[k] * usy[k] *
                               ((ftu[k][i][bj] - ftu[k][i][j]) -
                                (ftu[k][i][j] - ftu[k][i][aj])) +
                           0.5 * usy[k] * ucx[k] *
                               ((ftu[k][bi][bj] - ftu[k][bi][aj]) -
                                (ftu[k][ai][bj] - ftu[k][ai][aj]))) /
                          (dh * dh));  // DFP=D_para-D_perp
        nftu[k][i][j] += (DFQ *
                          (((ftu[k][bi][j] - ftu[k][i][j]) -
                            (ftu[k][i][j] - ftu[k][ai][j])) +
                           ((ftu[k][i][bj] - ftu[k][i][j]) -
                            (ftu[k][i][j] - ftu[k][i][aj]))) /
                          (dh * dh));  // DFQ=D_perp
        nftu[k][i][j] +=
            -v0 * (ucx[k] * (fmtu[k][bi][j] - fmtu[k][ai][j]) / (2.0 * dh) +
                   usy[k] * (fmtu[k][i][bj] - fmtu[k][i][aj]) /
                       (2.0 * dh));  // drift term.
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = dtt * nftu[k][i][j] / 3 + ftu[k][i][j];
      }

  /***************Updating to tmp variable*********************/
  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++) {
      for (k = 1; k <= M; k++) {
        rr[k] = nftu[k][i][j];
      }

      cyclic(a[i][j], b[i][j], c[i][j], alpha[i][j], beta[i][j], rr, rx, M);

      for (k = 1; k <= M; k++) {
        pftu[k][i][j] = rx[k];
      }
    }

  /*******************************************************************/
  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = 0;
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;

        if (ai <= 0) ai += Nx;

        if (bi > Nx) bi -= Nx;

        if (aj <= 0) aj += Ny;

        if (bj > Ny) bj -= Ny;

        if (ak < 1) ak += M;

        if (bk > M) bk -= M;

        nftu[k][i][j] += (DFP *
                          (ucx[k] * ucx[k] *
                               ((fmtu[k][bi][j] - fmtu[k][i][j]) -
                                (fmtu[k][i][j] - fmtu[k][ai][j])) +
                           usy[k] * usy[k] *
                               ((fmtu[k][i][bj] - fmtu[k][i][j]) -
                                (fmtu[k][i][j] - fmtu[k][i][aj])) +
                           0.5 * usy[k] * ucx[k] *
                               ((fmtu[k][bi][bj] - fmtu[k][bi][aj]) -
                                (fmtu[k][ai][bj] - fmtu[k][ai][aj]))) /
                          (dh * dh));  // DFP=D_para-D_perp
        nftu[k][i][j] += (DFQ *
                          (((fmtu[k][bi][j] - fmtu[k][i][j]) -
                            (fmtu[k][i][j] - fmtu[k][ai][j])) +
                           ((fmtu[k][i][bj] - fmtu[k][i][j]) -
                            (fmtu[k][i][j] - fmtu[k][i][aj]))) /
                          (dh * dh));  // DFQ=D_perp
        nftu[k][i][j] +=
            -v0 * (ucx[k] * (ftu[k][bi][j] - ftu[k][ai][j]) / (2.0 * dh) +
                   usy[k] * (ftu[k][i][bj] - ftu[k][i][aj]) /
                       (2.0 * dh));                        // drift term.
        nftu[k][i][j] += -2.0 * reversal * fmtu[k][i][j];  // reversal term.
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = dtt * nftu[k][i][j] / 3 + fmtu[k][i][j];
      }

  /***************Updating to tmp variable*********************/
  for (i = 1; i <= Nx; i++)
    for (j = 1; j <= Ny; j++) {
      for (k = 1; k <= M; k++) {
        rr[k] = nftu[k][i][j];
      }

      cyclic(a[i][j], b[i][j], c[i][j], alpha[i][j], beta[i][j], rr, rx, M);

      for (k = 1; k <= M; k++) {
        pfmtu[k][i][j] = rx[k];
      }
    }

  /*************************************************************/
  free_f3tensor(a, 1, Nx, 1, Ny, 1, M);
  free_f3tensor(b, 1, Nx, 1, Ny, 1, M);
  free_f3tensor(c, 1, Nx, 1, Ny, 1, M);
  free_matrix(alpha, 1, Nx, 1, Ny);
  free_matrix(beta, 1, Nx, 1, Ny);
  free_vector(rr, 1, M);
  free_vector(rx, 1, M);
  return 0;
}

int differ2(double ***ftu, double ***pftu, double ***fmtu, double ***pfmtu,
            double ***kru, double ***nftu) {
  int i, j, k;
  int ai, bi, aj, bj, ak, bk;
  double maxd = DFP;
  double ***a, ***b, ***c;
  double **alpha, **beta;
  double *rr, *rx;
  a = f3tensor(1, M, 1, Ny, 1, Nx);
  b = f3tensor(1, M, 1, Ny, 1, Nx);
  c = f3tensor(1, M, 1, Ny, 1, Nx);
  alpha = matrix(1, M, 1, Ny);
  beta = matrix(1, M, 1, Ny);
  rr = vector(1, Nx);
  rx = vector(1, Nx);

  /********************************************************************************/
  for (k = 1; k <= M; k++)
    for (j = 1; j <= Ny; j++)
      for (i = 1; i <= Nx; i++) {
        b[k][j][i] =
            (ucx[k] * ucx[k] * maxd + DFQ) * 2 * dtt / (3 * dh * dh) + 1;
      }

  for (k = 1; k <= M; k++)
    for (j = 1; j <= Ny; j++)
      for (i = 2; i <= Nx; i++) {
        a[k][j][i] = -(ucx[k] * ucx[k] * maxd + DFQ) * dtt / (3 * dh * dh);
      }

  for (k = 1; k <= M; k++)
    for (j = 1; j <= Ny; j++)
      for (i = 1; i <= Nx - 1; i++) {
        c[k][j][i] = -(ucx[k] * ucx[k] * maxd + DFQ) * dtt / (3 * dh * dh);
      }

  for (k = 1; k <= M; k++)
    for (j = 1; j <= Ny; j++) {
      alpha[k][j] = -(ucx[k] * ucx[k] * maxd + DFQ) * dtt / (3 * dh * dh);
      beta[k][j] = -(ucx[k] * ucx[k] * maxd + DFQ) * dtt / (3 * dh * dh);
    }

  /********************************************************************************/

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = 0;
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;

        if (ai <= 0) ai += Nx;

        if (bi > Nx) bi -= Nx;

        if (aj <= 0) aj += Ny;

        if (bj > Ny) bj -= Ny;

        if (ak < 1) ak += M;

        if (bk > M) bk -= M;

        nftu[k][i][j] += (DFP *
                          (usy[k] * usy[k] *
                               ((ftu[k][i][bj] - ftu[k][i][j]) -
                                (ftu[k][i][j] - ftu[k][i][aj])) +
                           0.5 * usy[k] * ucx[k] *
                               ((ftu[k][bi][bj] - ftu[k][bi][aj]) -
                                (ftu[k][ai][bj] - ftu[k][ai][aj]))) /
                          (dh * dh));
        nftu[k][i][j] += (DFQ *
                          (((ftu[k][i][bj] - ftu[k][i][j]) -
                            (ftu[k][i][j] - ftu[k][i][aj]))) /
                          (dh * dh));
        nftu[k][i][j] += ((DFR * ((ftu[bk][i][j] - ftu[k][i][j]) -
                                  (ftu[k][i][j] - ftu[ak][i][j]))) /
                          (dta * dta));
        nftu[k][i][j] += ((DFR * ((ftu[bk][i][j] + ftu[k][i][j]) *
                                      (kru[bk][i][j] - kru[k][i][j]) -
                                  (ftu[k][i][j] + ftu[ak][i][j]) *
                                      (kru[k][i][j] - kru[ak][i][j]))) /
                          (2 * dta * dta));
        nftu[k][i][j] +=
            -v0 * (ucx[k] * (fmtu[k][bi][j] - fmtu[k][ai][j]) / (2.0 * dh) +
                   usy[k] * (fmtu[k][i][bj] - fmtu[k][i][aj]) /
                       (2.0 * dh));  // drift term. Here Modified.
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = dtt * nftu[k][i][j] / 3 + ftu[k][i][j];
      }

  /***************Updating to tmp variable*********************/
  for (k = 1; k <= M; k++)
    for (j = 1; j <= Ny; j++) {
      for (i = 1; i <= Nx; i++) {
        rr[i] = nftu[k][i][j];
      }

      cyclic(a[k][j], b[k][j], c[k][j], alpha[k][j], beta[k][j], rr, rx, Nx);

      for (i = 1; i <= Nx; i++) {
        pftu[k][i][j] = rx[i];
      }
    }

  /********************************************************************************/
  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = 0;
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;

        if (ai <= 0) ai += Nx;

        if (bi > Nx) bi -= Nx;

        if (aj <= 0) aj += Ny;

        if (bj > Ny) bj -= Ny;

        if (ak < 1) ak += M;

        if (bk > M) bk -= M;

        nftu[k][i][j] += (DFP *
                          (usy[k] * usy[k] *
                               ((fmtu[k][i][bj] - fmtu[k][i][j]) -
                                (fmtu[k][i][j] - fmtu[k][i][aj])) +
                           0.5 * usy[k] * ucx[k] *
                               ((fmtu[k][bi][bj] - fmtu[k][bi][aj]) -
                                (fmtu[k][ai][bj] - fmtu[k][ai][aj]))) /
                          (dh * dh));
        nftu[k][i][j] += (DFQ *
                          (((fmtu[k][i][bj] - fmtu[k][i][j]) -
                            (fmtu[k][i][j] - fmtu[k][i][aj]))) /
                          (dh * dh));
        nftu[k][i][j] += ((DFR * ((fmtu[bk][i][j] - fmtu[k][i][j]) -
                                  (fmtu[k][i][j] - fmtu[ak][i][j]))) /
                          (dta * dta));
        nftu[k][i][j] += ((DFR * ((fmtu[bk][i][j] + fmtu[k][i][j]) *
                                      (kru[bk][i][j] - kru[k][i][j]) -
                                  (fmtu[k][i][j] + fmtu[ak][i][j]) *
                                      (kru[k][i][j] - kru[ak][i][j]))) /
                          (2 * dta * dta));
        nftu[k][i][j] +=
            -v0 * (ucx[k] * (ftu[k][bi][j] - ftu[k][ai][j]) / (2.0 * dh) +
                   usy[k] * (ftu[k][i][bj] - ftu[k][i][aj]) /
                       (2.0 * dh));                        // drift term.
        nftu[k][i][j] += -2.0 * reversal * fmtu[k][i][j];  // drift term.
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = dtt * nftu[k][i][j] / 3 + fmtu[k][i][j];
      }

  /***************Updating to tmp variable*********************/
  for (k = 1; k <= M; k++)
    for (j = 1; j <= Ny; j++) {
      for (i = 1; i <= Nx; i++) {
        rr[i] = nftu[k][i][j];
      }

      cyclic(a[k][j], b[k][j], c[k][j], alpha[k][j], beta[k][j], rr, rx, Nx);

      for (i = 1; i <= Nx; i++) {
        pfmtu[k][i][j] = rx[i];
      }
    }

  /*************************************************************/
  free_f3tensor(a, 1, M, 1, Ny, 1, Nx);
  free_f3tensor(b, 1, M, 1, Ny, 1, Nx);
  free_f3tensor(c, 1, M, 1, Ny, 1, Nx);
  free_matrix(alpha, 1, M, 1, Ny);
  free_matrix(beta, 1, M, 1, Ny);
  free_vector(rr, 1, Nx);
  free_vector(rx, 1, Nx);
  return 0;
}

int differ3(double ***ftu, double ***pftu, double ***fmtu, double ***pfmtu,
            double ***kru, double ***nftu) {
  int i, j, k;
  int ai, bi, aj, bj, ak, bk;
  double maxd = DFP;
  double ***a, ***b, ***c;
  double **alpha, **beta;
  double *rr, *rx;
  a = f3tensor(1, M, 1, Nx, 1, Ny);
  b = f3tensor(1, M, 1, Nx, 1, Ny);
  c = f3tensor(1, M, 1, Nx, 1, Ny);
  alpha = matrix(1, M, 1, Nx);
  beta = matrix(1, M, 1, Nx);
  rr = vector(1, Ny);
  rx = vector(1, Ny);

  /***************************************************************************************/
  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        b[k][i][j] =
            (usy[k] * usy[k] * maxd + DFQ) * 2 * dtt / (3 * dh * dh) + 1;
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 2; j <= Ny; j++) {
        a[k][i][j] = -(usy[k] * usy[k] * maxd + DFQ) * dtt / (3 * dh * dh);
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny - 1; j++) {
        c[k][i][j] = -(usy[k] * usy[k] * maxd + DFQ) * dtt / (3 * dh * dh);
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++) {
      alpha[k][i] = -(usy[k] * usy[k] * maxd + DFQ) * dtt / (3 * dh * dh);
      beta[k][i] = -(usy[k] * usy[k] * maxd + DFQ) * dtt / (3 * dh * dh);
    }

  /***************************************************************************************/
  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = 0;
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;

        if (ai <= 0) ai += Nx;

        if (bi > Nx) bi -= Nx;

        if (aj <= 0) aj += Ny;

        if (bj > Ny) bj -= Ny;

        if (ak < 1) ak += M;

        if (bk > M) bk -= M;

        nftu[k][i][j] += (DFP *
                          (ucx[k] * ucx[k] *
                               ((ftu[k][bi][j] - ftu[k][i][j]) -
                                (ftu[k][i][j] - ftu[k][ai][j])) +
                           0.5 * usy[k] * ucx[k] *
                               ((ftu[k][bi][bj] - ftu[k][bi][aj]) -
                                (ftu[k][ai][bj] - ftu[k][ai][aj]))) /
                          (dh * dh));
        nftu[k][i][j] += (DFQ *
                          (((ftu[k][bi][j] - ftu[k][i][j]) -
                            (ftu[k][i][j] - ftu[k][ai][j]))) /
                          (dh * dh));
        nftu[k][i][j] += ((DFR * ((ftu[bk][i][j] - ftu[k][i][j]) -
                                  (ftu[k][i][j] - ftu[ak][i][j]))) /
                          (dta * dta));
        nftu[k][i][j] += ((DFR * ((ftu[bk][i][j] + ftu[k][i][j]) *
                                      (kru[bk][i][j] - kru[k][i][j]) -
                                  (ftu[k][i][j] + ftu[ak][i][j]) *
                                      (kru[k][i][j] - kru[ak][i][j]))) /
                          (2 * dta * dta));
        nftu[k][i][j] +=
            -v0 * (ucx[k] * (fmtu[k][bi][j] - fmtu[k][ai][j]) / (2.0 * dh) +
                   usy[k] * (fmtu[k][i][bj] - fmtu[k][i][aj]) /
                       (2.0 * dh));  // drift term.
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = dtt * nftu[k][i][j] / 3 + ftu[k][i][j];
      }

  /***************Updating to tmp variable*********************/
  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++) {
      for (j = 1; j <= Ny; j++) {
        rr[j] = nftu[k][i][j];
      }

      cyclic(a[k][i], b[k][i], c[k][i], alpha[k][i], beta[k][i], rr, rx, Ny);

      for (j = 1; j <= Ny; j++) {
        pftu[k][i][j] = rx[j];
      }
    }

  /***************************************************************************************/
  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = 0;
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;

        if (ai <= 0) ai += Nx;

        if (bi > Nx) bi -= Nx;

        if (aj <= 0) aj += Ny;

        if (bj > Ny) bj -= Ny;

        if (ak < 1) ak += M;

        if (bk > M) bk -= M;

        nftu[k][i][j] += (DFP *
                          (ucx[k] * ucx[k] *
                               ((fmtu[k][bi][j] - fmtu[k][i][j]) -
                                (fmtu[k][i][j] - fmtu[k][ai][j])) +
                           0.5 * usy[k] * ucx[k] *
                               ((fmtu[k][bi][bj] - fmtu[k][bi][aj]) -
                                (fmtu[k][ai][bj] - fmtu[k][ai][aj]))) /
                          (dh * dh));
        nftu[k][i][j] += (DFQ *
                          (((fmtu[k][bi][j] - fmtu[k][i][j]) -
                            (fmtu[k][i][j] - fmtu[k][ai][j]))) /
                          (dh * dh));
        nftu[k][i][j] += ((DFR * ((fmtu[bk][i][j] - fmtu[k][i][j]) -
                                  (fmtu[k][i][j] - fmtu[ak][i][j]))) /
                          (dta * dta));
        nftu[k][i][j] += ((DFR * ((fmtu[bk][i][j] + fmtu[k][i][j]) *
                                      (kru[bk][i][j] - kru[k][i][j]) -
                                  (fmtu[k][i][j] + fmtu[ak][i][j]) *
                                      (kru[k][i][j] - kru[ak][i][j]))) /
                          (2 * dta * dta));
        nftu[k][i][j] +=
            -v0 * (ucx[k] * (ftu[k][bi][j] - ftu[k][ai][j]) / (2.0 * dh) +
                   usy[k] * (ftu[k][i][bj] - ftu[k][i][aj]) /
                       (2.0 * dh));                        // drift term.
        nftu[k][i][j] += -2.0 * reversal * fmtu[k][i][j];  // reversal term.
      }

  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        nftu[k][i][j] = dtt * nftu[k][i][j] / 3 + fmtu[k][i][j];
      }

  /***************Updating to tmp variable*********************/
  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++) {
      for (j = 1; j <= Ny; j++) {
        rr[j] = nftu[k][i][j];
      }

      cyclic(a[k][i], b[k][i], c[k][i], alpha[k][i], beta[k][i], rr, rx, Ny);

      for (j = 1; j <= Ny; j++) {
        pfmtu[k][i][j] = rx[j];
      }
    }

  /*************************************************************/
  free_f3tensor(a, 1, M, 1, Nx, 1, Ny);
  free_f3tensor(b, 1, M, 1, Nx, 1, Ny);
  free_f3tensor(c, 1, M, 1, Nx, 1, Ny);
  free_matrix(alpha, 1, M, 1, Nx);
  free_matrix(beta, 1, M, 1, Nx);
  free_vector(rr, 1, Ny);
  free_vector(rx, 1, Ny);
  return 0;
}

void alloc_cyclic_arrays(double ***a, double ***b, double ***c, double **alpha,
                         double **beta, double *rr, double *rx, int dim1,
                         int dim2, int dim3) {
  a = f3tensor(1, dim1, 1, dim2, 1, dim3);
  b = f3tensor(1, dim1, 1, dim2, 1, dim3);
  c = f3tensor(1, dim1, 1, dim2, 1, dim3);
  alpha = matrix(1, dim1, 1, dim2);
  beta = matrix(1, dim1, 1, dim2);
  rr = vector(1, dim3);
  rx = vector(1, dim3);
}

void dealloc_cyclic_arrays(double ***a, double ***b, double ***c,
                           double **alpha, double **beta, double *rr,
                           double *rx, int dim1, int dim2, int dim3) {
  free_f3tensor(a, 1, dim1, 1, dim2, 1, dim3);
  free_f3tensor(b, 1, dim1, 1, dim2, 1, dim3);
  free_f3tensor(c, 1, dim1, 1, dim2, 1, dim3);
  free_matrix(alpha, 1, dim1, 1, dim2);
  free_matrix(beta, 1, dim1, 1, dim2);
  free_vector(rr, 1, dim3);
  free_vector(rx, 1, dim3);
}

void updating(double ***ftu, double ***pftu, double ***fmtu, double ***pfmtu) {
  int i, j, k;

  /*Updating f_T f_M together.*/
  for (k = 1; k <= M; k++)
    for (i = 1; i <= Nx; i++)
      for (j = 1; j <= Ny; j++) {
        ftu[k][i][j] = pftu[k][i][j];
        fmtu[k][i][j] = pfmtu[k][i][j];
        pftu[k][i][j] = 0;
        pfmtu[k][i][j] = 0;
      }
}

int interation(double ***ftu, double ***pftu, double ***fmtu, double ***pfmtu,
               double *x1, double *x2, int m) {
  double ***kru, ***nftu;
  kru = f3tensor(1, M, 1, Nx, 1, Ny);
  nftu = f3tensor(1, M, 1, Nx, 1, Ny);
  /****************************Theta******************************************/
  int_kp_simple(ftu, kru, x1, x2);             // Calculate W(r, u).
  differ1(ftu, pftu, fmtu, pfmtu, kru, nftu);  // differential in theta.
  updating(ftu, pftu, fmtu, pfmtu);
  /***************************Y-direction*************************************/
  int_kp_simple(ftu, kru, x1, x2);
  differ2(ftu, pftu, fmtu, pfmtu, kru, nftu);
  updating(ftu, pftu, fmtu, pfmtu);
  /****************************Theta*******************************************/
  int_kp_simple(ftu, kru, x1, x2);
  differ1(ftu, pftu, fmtu, pfmtu, kru, nftu);
  updating(ftu, pftu, fmtu, pfmtu);
  /****************************X-direction*************************************/
  int_kp_simple(ftu, kru, x1, x2);
  differ3(ftu, pftu, fmtu, pfmtu, kru, nftu);
  updating(ftu, pftu, fmtu, pfmtu);

  /*********************************Output*************************************/

  if ((m - 1) % INTER1 == 0) {
    output(m, ftu, fmtu, kru, x1, x2);
  }

  free_f3tensor(kru, 1, M, 1, Nx, 1, Ny);
  free_f3tensor(nftu, 1, M, 1, Nx, 1, Ny);
  return 0;
}

void simulation() {
  int m;
  double *x1, *x2;
  double ***ftu, ***pftu;
  double ***fmtu, ***pfmtu;
  x1 = vector(1, Nx);
  x2 = vector(1, Ny);
  ftu = f3tensor(1, M, 1, Nx, 1, Ny);
  pftu = f3tensor(1, M, 1, Nx, 1, Ny);
  fmtu = f3tensor(1, M, 1, Nx, 1, Ny);
  pfmtu = f3tensor(1, M, 1, Nx, 1, Ny);
  //initial(x1, x2, ftu, pftu, fmtu, pfmtu);
  initial_read_band(x1, x2, ftu, pftu, fmtu, pfmtu);

  //  initial_continue(x1, x2, ftu, pftu, fmtu, pfmtu);
  for (m = 1; m <= TOTALTIME; m++) {
    nowtime = m;
#ifdef _MSC_VER
    if ((m - 1) % INTER2 == 0) {
      outscreen(m, ftu, fmtu);
    }
#endif  // _MSC_VER

    interation(ftu, pftu, fmtu, pfmtu, x1, x2, nowtime);
  }
  free_vector(x1, 1, Nx);
  free_vector(x2, 1, Ny);
  free_f3tensor(ftu, 1, M, 1, Nx, 1, Ny);
  free_f3tensor(pftu, 1, M, 1, Nx, 1, Ny);
  free_f3tensor(fmtu, 1, M, 1, Nx, 1, Ny);
  free_f3tensor(pfmtu, 1, M, 1, Nx, 1, Ny);
}

int main() {
  make_subdirectory();
  ReadParameter();
  simulation();
  return 0;
}
