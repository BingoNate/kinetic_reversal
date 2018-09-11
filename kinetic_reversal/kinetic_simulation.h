#ifndef KINETIC_SIMULATION_H
#define KINETIC_SIMULATION_H

#include "kinetic_io.h"
#include "kinetic_parameters.h"

void calculate_potential_w(const Mat3DDoub &ftu, Mat3DDoub &kru) {
  int i, j, k1, k2;
  for (k1 = 0; k1 < M; k1++) {
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        kru[k1][i][j] = 0;
      }
  }

  for (k1 = 0; k1 < M; k1++) {
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        for (k2 = 0; k2 < M; k2++) {
          if ((k2 == k1) || fabs(k2 - k1) == M) {
            continue;
          }
          kru[k1][i][j] += UXU[k1][k2] * ftu[k2][i][j];
        }
      }
  }

  for (k1 = 0; k1 < M; k1++) {
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        kru[k1][i][j] = kru[k1][i][j] * dta * 2;  // W(r,u).
      }
  }
}

/*First in theta.*/
int differ_theta(Mat3DDoub &ftu, Mat3DDoub &pftu, Mat3DDoub &fmtu,
                 Mat3DDoub &pfmtu, Mat3DDoub &kru, Mat3DDoub &nftu) {
  int i, j, k;
  int ai, bi;
  int aj, bj;
  int ak, bk;

  Mat3DDoub a(Nx, Ny, M);
  Mat3DDoub b(Nx, Ny, M);
  Mat3DDoub c(Nx, Ny, M);
  MatDoub alpha(Nx, Ny);
  MatDoub beta(Nx, Ny);

  for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++)
      for (k = 0; k < M; k++) {
        ak = k - 1;
        bk = k + 1;
        if (ak < 0) ak += M;
        if (bk == M) bk -= M;

        b[i][j][k] = 2 * DFR * dtt / (3 * dta * dta) -
                     DFR * (kru[bk][i][j] - 2 * kru[k][i][j] + kru[ak][i][j]) *
                         dtt / (6 * dta * dta) +
                     1;
      }

  for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++)
      for (k = 1; k < M; k++) {
        ak = k - 1;
        bk = k + 1;
        if (ak < 0) ak += M;
        if (bk == M) bk -= M;

        a[i][j][k] =
            -DFR * dtt / (3 * dta * dta) -
            DFR * (kru[ak][i][j] - kru[k][i][j]) * dtt / (6 * dta * dta);
      }

  for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++)
      for (k = 0; k < M - 1; k++) {
        ak = k - 1;
        bk = k + 1;
        if (ak < 0) ak += M;
        if (bk == M) bk -= M;

        c[i][j][k] =
            -DFR * dtt / (3 * dta * dta) -
            DFR * (kru[bk][i][j] - kru[k][i][j]) * dtt / (6 * dta * dta);
      }

  for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++) {
      k = M - 1;
      ak = k - 1;
      bk = k + 1;
      if (ak < 0) ak += M;
      if (bk == M) bk -= M;

      alpha[i][j] =
          -DFR * dtt / (3 * dta * dta) -
          DFR * (kru[bk][i][j] - kru[k][i][j]) * dtt / (6 * dta * dta);
      k = 0;
      ak = k - 1;
      bk = k + 1;

      if (ak < 0) ak += M;
      if (bk == M) bk -= M;

      beta[i][j] = -DFR * dtt / (3 * dta * dta) -
                   DFR * (kru[ak][i][j] - kru[k][i][j]) * dtt / (6 * dta * dta);
    }

  /*****f(r,u,t)*****/
  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        //pftu[k][i][j] = ftu[k][i][j];  
        //pfmtu[k][i][j] = fmtu[k][i][j];
        nftu[k][i][j] = 0;
      }

  /*Diffusion terms.*/
  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;

        if (ai < 0) ai += Nx;
        if (bi >= Nx) bi -= Nx;
        if (aj < 0) aj += Ny;
        if (bj >= Ny) bj -= Ny;
        if (ak < 0) ak += M;
        if (bk >= M) bk -= M;

        nftu[k][i][j] +=
            DFP *
            (ucx[k] * ucx[k] *
                 (ftu[k][bi][j] - 2.0 * ftu[k][i][j] + ftu[k][ai][j]) +
             usy[k] * usy[k] *
                 (ftu[k][i][bj] - 2.0 * ftu[k][i][j] + ftu[k][i][aj]) +
             0.5 * usy[k] * ucx[k] *
                 (ftu[k][bi][bj] - ftu[k][bi][aj] - ftu[k][ai][bj] +
                  ftu[k][ai][aj])) /
            (dxx * dyy);  // DFP=D_para-D_perp

        nftu[k][i][j] +=
            DFQ *
            ((ftu[k][bi][j] - 2.0 * ftu[k][i][j] + ftu[k][ai][j]) +
             (ftu[k][i][bj] - 2.0 * ftu[k][i][j] + ftu[k][i][aj])) /
            (dxx * dyy);  // DFQ=D_perp

        nftu[k][i][j] +=
            -v0 * (ucx[k] * (fmtu[k][bi][j] - fmtu[k][ai][j]) / (2.0 * dxx) +
                   usy[k] * (fmtu[k][i][bj] - fmtu[k][i][aj]) /
                       (2.0 * dyy));  // Drift term.
      }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        nftu[k][i][j] =
            dtt * nftu[k][i][j] / 3 + ftu[k][i][j];  // 1/4 timestep.
      }

  vector<double> rr(M);
  vector<double> rx(M);
  vector<double> aa(M);
  vector<double> bb(M);
  vector<double> cc(M);

  for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++)
      for (k = 0; k < M; k++) {
        aa[k] = a[i][j][k];
        bb[k] = b[i][j][k];
        cc[k] = c[i][j][k];
      }

  for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++) {
      for (k = 0; k < M; k++) {
        rr[k] = nftu[k][i][j];
      }

      cyclic(aa, bb, cc, alpha[i][j], beta[i][j], rr, rx);

      for (k = 0; k < M; k++) {
        pftu[k][i][j] = rx[k];  // Solve the equation for f(r,u,t). Updating it.
      }
    }

  /********f_M(r,u,t)******/
  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        nftu[k][i][j] = 0;
      }

  /*Diffusion terms for f_M.*/
  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;
        if (ai < 0) ai += Nx;
        if (bi >= Nx) bi -= Nx;
        if (aj < 0) aj += Ny;
        if (bj >= Ny) bj -= Ny;
        if (ak < 0) ak += M;
        if (bk >= M) bk -= M;

        nftu[k][i][j] +=
            DFP *
            (ucx[k] * ucx[k] *
                 (fmtu[k][bi][j] - 2.0 * fmtu[k][i][j] + fmtu[k][ai][j]) +
             usy[k] * usy[k] *
                 (fmtu[k][i][bj] - 2.0 * fmtu[k][i][j] + fmtu[k][i][aj]) +
             0.5 * usy[k] * ucx[k] *
                 (fmtu[k][bi][bj] - fmtu[k][bi][aj] - fmtu[k][ai][bj] +
                  fmtu[k][ai][aj])) /
            (dxx * dyy);  // DFP=D_para-D_perp

        nftu[k][i][j] +=
            DFQ *
            ((fmtu[k][bi][j] - 2.0 * fmtu[k][i][j] + fmtu[k][ai][j]) +
             (fmtu[k][i][bj] - 2.0 * fmtu[k][i][j] + fmtu[k][i][aj])) /
            (dxx * dyy);  // DFQ=D_perp

        nftu[k][i][j] +=
            -v0 * (ucx[k] * (ftu[k][bi][j] - ftu[k][ai][j]) / (2.0 * dxx) +
                   usy[k] * (ftu[k][i][bj] - ftu[k][i][aj]) /
                       (2.0 * dyy));                        // Drift term.
        nftu[k][i][j] += -2.0 * reversal * fmtu[k][i][j];  // reversal term.
      }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        nftu[k][i][j] =
            dtt * nftu[k][i][j] / 3 + fmtu[k][i][j];  // 1/4 timestep.
      }

  rr.clear();
  rx.clear();
  for (i = 0; i < Nx; i++)
    for (j = 0; j < Ny; j++) {
      for (k = 0; k < M; k++) {
        rr[k] = nftu[k][i][j];
      }
      cyclic(aa, bb, cc, alpha[i][j], beta[i][j], rr, rx);

      for (k = 0; k < M; k++) {
        pfmtu[k][i][j] = rx[k];  // Solve the equation for f_M.
      }
    }
  return 0;
}

int differ_x(Mat3DDoub &ftu, Mat3DDoub &pftu, Mat3DDoub &fmtu, Mat3DDoub &pfmtu,
             Mat3DDoub &kru, Mat3DDoub &nftu) {
  int i, j, k;
  int ai, bi;
  int aj, bj;
  int ak, bk;

  Mat3DDoub a(M, Ny, Nx);
  Mat3DDoub b(M, Ny, Nx);
  Mat3DDoub c(M, Ny, Nx);
  MatDoub alpha(M, Ny);
  MatDoub beta(M, Ny);

  for (k = 0; k < M; k++)
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++) {
        b[k][j][i] =
            (ucx[k] * ucx[k] * DFP + DFQ) * 2 * dtt / (3 * dxx * dyy) + 1;
      }

  for (k = 0; k < M; k++)
    for (j = 0; j < Ny; j++)
      for (i = 1; i < Nx; i++) {
        a[k][j][i] = -(ucx[k] * ucx[k] * DFP + DFQ) * dtt / (3 * dxx * dyy);
      }

  for (k = 0; k < M; k++)
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx - 1; i++) {
        c[k][j][i] = -(ucx[k] * ucx[k] * DFP + DFQ) * dtt / (3 * dxx * dyy);
      }

  for (k = 0; k < M; k++)
    for (j = 0; j < Ny; j++)  // Actually N here.
    {
      alpha[k][j] = -(ucx[k] * ucx[k] * DFP + DFQ) * dtt / (3 * dxx * dyy);
      beta[k][j] = -(ucx[k] * ucx[k] * DFP + DFQ) * dtt / (3 * dxx * dyy);
    }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        pftu[k][i][j] = ftu[k][i][j];
        pfmtu[k][i][j] = fmtu[k][i][j];
        nftu[k][i][j] = 0;
      }

  /***********f(r,u,t)***********/
  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;
        if (ai < 0) ai += Nx;
        if (bi >= Nx) bi -= Nx;
        if (aj < 0) aj += Ny;
        if (bj >= Ny) bj -= Ny;
        if (ak < 0) ak += M;
        if (bk >= M) bk -= M;

        nftu[k][i][j] +=
            DFP *
            (usy[k] * usy[k] *
                 (ftu[k][i][bj] - 2.0 * ftu[k][i][j] + ftu[k][i][aj]) +
             0.5 * usy[k] * ucx[k] *
                 (ftu[k][bi][bj] - ftu[k][bi][aj] - ftu[k][ai][bj] +
                  ftu[k][ai][aj])) /
            (dxx * dyy);
        nftu[k][i][j] +=
            DFQ * (ftu[k][i][bj] - 2.0 * ftu[k][i][j] + ftu[k][i][aj]) /
            (dyy * dyy);
        nftu[k][i][j] +=
            -v0 *
            (usy[k] * (fmtu[k][i][bj] - fmtu[k][i][aj]) / (2.0 * dyy) +
             ucx[k] * (fmtu[k][bi][j] - fmtu[k][ai][j]) /
                 (2.0 * dxx));  // Add drift term here: v0(u.\nabla)f_M(r,u,t).
                                // Only \partial_y.
        nftu[k][i][j] +=
            DFR * (ftu[bk][i][j] - 2.0 * ftu[k][i][j] + ftu[ak][i][j]) /
            (dta * dta);  // Actually, it is just linear term for
                          // \partial_\theta^{2} {f(r,u,t)}.
        nftu[k][i][j] +=
            DFR *
            ((ftu[bk][i][j] + ftu[k][i][j]) * (kru[bk][i][j] - kru[k][i][j]) -
             (ftu[k][i][j] + ftu[ak][i][j]) *
                 (kru[k][i][j] - kru[ak][i][j])) /
            (2 * dta *
             dta);  // This is non-linear term for interaction potential.
      }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        nftu[k][i][j] = dtt * nftu[k][i][j] / 3 + ftu[k][i][j];
      }

  /*y direction.*/
  vector<double> rr(Nx);
  vector<double> rx(Nx);
  vector<double> aa(Nx);
  vector<double> bb(Nx);
  vector<double> cc(Nx);

  for (k = 0; k < M; k++)
    for (j = 0; j < Ny; j++)
      for (i = 0; i < Nx; i++) {
        aa[i] = a[k][j][i];
        bb[i] = b[k][j][i];
        cc[i] = c[k][j][i];
      }

  for (k = 0; k < M; k++)
    for (j = 0; j < Ny; j++) {
      for (i = 0; i < Nx; i++) {
        rr[i] = nftu[k][i][j];
      }

      cyclic(aa, bb, cc, alpha[k][j], beta[k][j], rr, rx);

      for (i = 0; i < Nx; i++) {
        pftu[k][i][j] = rx[i];  // Updating f(r,u,t).
      }
    }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        nftu[k][i][j] = 0;  // Calculate for rr set to 0.
      }

  /*************f_M(r,u,t)************/
  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;
        if (ai < 0) ai += Nx;
        if (bi >= Nx) bi -= Nx;
        if (aj < 0) aj += Ny;
        if (bj >= Ny) bj -= Ny;
        if (ak < 0) ak += M;
        if (bk >= M) bk -= M;

        nftu[k][i][j] +=
            DFP *
            (usy[k] * usy[k] *
                 (fmtu[k][i][bj] - 2.0 * fmtu[k][i][j] + fmtu[k][i][aj]) +
             0.5 * usy[k] * ucx[k] *
                 (fmtu[k][bi][bj] - fmtu[k][bi][aj] - fmtu[k][ai][bj] +
                  fmtu[k][ai][aj])) /
            (dxx * dyy);
        nftu[k][i][j] +=
            DFQ * (fmtu[k][i][bj] - 2.0 * fmtu[k][i][j] + fmtu[k][i][aj]) /
            (dyy * dyy);
        nftu[k][i][j] +=
            -v0 *
            (usy[k] * (ftu[k][i][bj] - ftu[k][i][aj]) / (2.0 * dyy) +
             ucx[k] * (ftu[k][bi][j] - ftu[k][ai][j]) /
                 (2.0 * dxx));  // Add drift term here: v0(u.\nabla)f_M(r,u,t).
                                // Only \partial_y.
        nftu[k][i][j] +=
            -2.0 * reversal * fmtu[k][i][j];  // Add reversal term.
        nftu[k][i][j] +=
            DFR * (fmtu[bk][i][j] - 2.0 * fmtu[k][i][j] + fmtu[ak][i][j]) /
            (dta * dta);  // Actually, it is just linear term for
                          // \partial_\theta^{2} {f(r,u,t)}.
        nftu[k][i][j] +=
            DFR *
            ((fmtu[bk][i][j] + fmtu[k][i][j]) *
                 (kru[bk][i][j] - kru[k][i][j]) -
             (fmtu[k][i][j] + fmtu[ak][i][j]) *
                 (kru[k][i][j] - kru[ak][i][j])) /
            (2 * dta *
             dta);  // This is non-linear term for interaction potential.
      }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        nftu[k][i][j] = dtt * nftu[k][i][j] / 3 + fmtu[k][i][j];
      }

  rr.clear();
  rx.clear();
  for (k = 0; k < M; k++)
    for (j = 0; j < Ny; j++) {
      for (i = 0; i < Nx; i++) {
        rr[i] = nftu[k][i][j];
      }

      cyclic(aa, bb, cc, alpha[k][j], beta[k][j], rr, rx);

      for (i = 0; i < Nx; i++) {
        pfmtu[k][i][j] = rx[i];  // Updating f_M(r,u,t).
      }
    }
  return 0;
}

int differ_y(Mat3DDoub &ftu, Mat3DDoub &pftu, Mat3DDoub &fmtu, Mat3DDoub &pfmtu,
             Mat3DDoub &kru, Mat3DDoub &nftu) {
  int i, j, k;
  int ai, bi;
  int aj, bj;
  int ak, bk;

  Mat3DDoub a(M, Nx, Ny);
  Mat3DDoub b(M, Nx, Ny);
  Mat3DDoub c(M, Nx, Ny);
  MatDoub alpha(M, Nx);
  MatDoub beta(M, Nx);

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        b[k][i][j] =
            (usy[k] * usy[k] * DFP + DFQ) * 2 * dtt / (3 * dxx * dyy) + 1;
      }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 1; j < Ny; j++) {
        a[k][i][j] = -(usy[k] * usy[k] * DFP + DFQ) * dtt / (3 * dxx * dyy);
      }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny - 1; j++) {
        c[k][i][j] = -(usy[k] * usy[k] * DFP + DFQ) * dtt / (3 * dxx * dyy);
      }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)  // Actually N here.
    {
      alpha[k][i] = -(usy[k] * usy[k] * DFP + DFQ) * dtt / (3 * dxx * dyy);
      beta[k][i] = -(usy[k] * usy[k] * DFP + DFQ) * dtt / (3 * dxx * dyy);
    }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        //pftu[k][i][j] = ftu[k][i][j];
        //pfmtu[k][i][j] = fmtu[k][i][j];
        nftu[k][i][j] = 0;
      }

  /**********f(r,u,t)*********/
  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;
        if (ai < 0) ai += Nx;
        if (bi >= Nx) bi -= Nx;
        if (aj < 0) aj += Ny;
        if (bj >= Ny) bj -= Ny;
        if (ak < 0) ak += M;
        if (bk >= M) bk -= M;

        nftu[k][i][j] +=
            DFP *
            (ucx[k] * ucx[k] *
                 (ftu[k][bi][j] - 2.0 * ftu[k][i][j] + ftu[k][ai][j]) +
             0.5 * usy[k] * ucx[k] *
                 (ftu[k][bi][bj] - ftu[k][bi][aj] - ftu[k][ai][bj] +
                  ftu[k][ai][aj])) /
            (dxx * dyy);
        nftu[k][i][j] +=
            DFQ * (ftu[k][bi][j] - 2.0 * ftu[k][i][j] + ftu[k][ai][j]) /
            (dxx * dxx);
        nftu[k][i][j] +=
            -v0 *
            (usy[k] * (fmtu[k][i][bj] - fmtu[k][i][aj]) / (2.0 * dyy) +
             ucx[k] * (fmtu[k][bi][j] - fmtu[k][ai][j]) /
                 (2.0 * dxx));  // Add drift term here: v0(u.\nabla)f_M(r,u,t).
                                // Only \partial_x.
        nftu[k][i][j] +=
            DFR * (ftu[bk][i][j] - 2.0 * ftu[k][i][j] + ftu[ak][i][j]) /
            (dta * dta);
        nftu[k][i][j] +=
            DFR *
            ((ftu[bk][i][j] + ftu[k][i][j]) * (kru[bk][i][j] - kru[k][i][j]) -
             (ftu[k][i][j] + ftu[ak][i][j]) *
                 (kru[k][i][j] - kru[ak][i][j])) /
            (2 * dta * dta);
      }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        nftu[k][i][j] = dtt * nftu[k][i][j] / 3 + ftu[k][i][j];
      }

  vector<double> rr(Ny);
  vector<double> rx(Ny);
  vector<double> aa(Ny);
  vector<double> bb(Ny);
  vector<double> cc(Ny);

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        aa[j] = a[k][i][j];
        bb[j] = b[k][i][j];
        cc[j] = c[k][i][j];
      }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++) {
      for (j = 0; j < Ny; j++) {
        rr[j] = nftu[k][i][j];
      }

      cyclic(aa, bb, cc, alpha[k][i], beta[k][i], rr, rx);

      for (j = 0; j < Ny; j++) {
        pftu[k][i][j] = rx[j];  // Updating f(r,u,t).
      }
    }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        nftu[k][i][j] = 0;  // Set to 0.
      }

  /********f_M(r,u,t)********/
  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        ai = i - 1;
        bi = i + 1;
        aj = j - 1;
        bj = j + 1;
        ak = k - 1;
        bk = k + 1;
        if (ai < 0) ai += Nx;
        if (bi >= Nx) bi -= Nx;
        if (aj < 0) aj += Ny;
        if (bj >= Ny) bj -= Ny;
        if (ak < 0) ak += M;
        if (bk >= M) bk -= M;

        nftu[k][i][j] +=
            DFP *
            (ucx[k] * ucx[k] *
                 (fmtu[k][bi][j] - 2.0 * fmtu[k][i][j] + fmtu[k][ai][j]) +
             0.5 * usy[k] * ucx[k] *
                 (fmtu[k][bi][bj] - fmtu[k][bi][aj] - fmtu[k][ai][bj] +
                  fmtu[k][ai][aj])) /
            (dxx * dyy);
        nftu[k][i][j] +=
            DFQ * (fmtu[k][bi][j] - 2.0 * fmtu[k][i][j] + fmtu[k][ai][j]) /
            (dxx * dxx);
        nftu[k][i][j] +=
            -v0 *
            (usy[k] * (ftu[k][i][bj] - ftu[k][i][aj]) / (2.0 * dyy) +
             ucx[k] * (ftu[k][bi][j] - ftu[k][ai][j]) /
                 (2.0 * dxx));  // Add drift term here: v0(u.\nabla)f_M(r,u,t).
                                // Only \partial_x.
        nftu[k][i][j] +=
            -2.0 * reversal * fmtu[k][i][j];  // Add reversal term.
        nftu[k][i][j] +=
            DFR * (fmtu[bk][i][j] - 2.0 * fmtu[k][i][j] + fmtu[ak][i][j]) /
            (dta * dta);
        nftu[k][i][j] += DFR *
                         ((fmtu[bk][i][j] + fmtu[k][i][j]) *
                              (kru[bk][i][j] - kru[k][i][j]) -
                          (fmtu[k][i][j] + fmtu[ak][i][j]) *
                              (kru[k][i][j] - kru[ak][i][j])) /
                         (2 * dta * dta);
      }

  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++)
      for (j = 0; j < Ny; j++) {
        nftu[k][i][j] = dtt * nftu[k][i][j] / 3 + fmtu[k][i][j];
      }

  rr.clear();
  rx.clear();
  for (k = 0; k < M; k++)
    for (i = 0; i < Nx; i++) {
      for (j = 0; j < Ny; j++) {
        rr[j] = nftu[k][i][j];
      }

      cyclic(aa, bb, cc, alpha[k][i], beta[k][i], rr, rx);

      for (j = 0; j < Ny; j++) {
        pfmtu[k][i][j] = rx[j];  // Updating f_M(r,u,t).
      }
    }
  return 0;
}

void updating(Mat3DDoub &ftu, Mat3DDoub &pftu, Mat3DDoub &fmtu,
  Mat3DDoub &pfmtu)
{
  for (int k = 0; k < M; k++) {
    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
        ftu[k][i][j] = pftu[k][i][j];
        fmtu[k][i][j] = pfmtu[k][i][j];
        pftu[k][i][j] = 0;
        pfmtu[k][i][j] = 0;
      }
    }
  }
}

void evolution(int m_step, Mat3DDoub& ftu, Mat3DDoub& pftu, Mat3DDoub& fmtu, Mat3DDoub& pfmtu,
               vector<double>& x1, vector<double>& x2){
  Mat3DDoub kru(M, Nx, Ny);
  Mat3DDoub nftu(M, Nx, Ny);
  calculate_potential_w(ftu, kru);
  differ_theta(ftu, pftu, fmtu, pfmtu, kru, nftu);
  updating(ftu, pftu, fmtu, pfmtu);
  calculate_potential_w(ftu, kru);
  differ_x(ftu, pftu, fmtu, pfmtu, kru, nftu);
  updating(ftu, pftu, fmtu, pfmtu);
  // calculate_potential_w(ftu, kru);
  //differ_theta(ftu, pftu, fmtu, pfmtu, kru, nftu);
  //updating(ftu, pftu, fmtu, pfmtu);
  // calculate_potential_w(ftu, kru);
  //differ_y(ftu, pftu, fmtu, pfmtu, kru, nftu);
  //updating(ftu, pftu, fmtu, pfmtu);
  // Determine which step to output.
  //output(m_step, ftu, fmtu, kru, x1, x2);
}

void simulation() { 
  vector<double> x1(Nx);
  vector<double> x2(Ny);
  Mat3DDoub ftu(M, Nx, Ny);
  Mat3DDoub pftu(M, Nx, Ny);
  Mat3DDoub fmtu(M, Nx, Ny);
  Mat3DDoub pfmtu(M, Nx, Ny);
  initial_set(set_q, x1, x2, ftu, pftu, fmtu, pfmtu);
  for (int nowtime = 0; nowtime < TOTALTIME; nowtime++) {
    if (nowtime % INTER2 == 0) {
      outscreen(nowtime, ftu, fmtu);
    }
    evolution(nowtime, ftu, pftu, fmtu, pfmtu, x1, x2);
  }
}

#endif  // !KINETIC_SIMULATION_H
