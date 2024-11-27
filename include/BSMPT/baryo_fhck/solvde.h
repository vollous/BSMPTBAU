// https://github.com/NJdevPro/Numerical-Recipes

#pragma once

#include <BSMPT/baryo_fhck/difeq.h>
#include <BSMPT/utility/data_structures.h>
#include <fstream>

struct Solvde
{
  const int itmax;
  const double conv;
  const double slowc;
  const VecDoub &scalv;
  VecInt &indexv;
  const int nb;
  MatDoub &y;
  Difeq &difeq;
  int ne, m;
  VecInt kmax;
  VecDoub ermax;
  Mat3DDoub c;
  MatDoub s;

  inline Solvde(const int itmaxx,
                const double convv,
                const double slowcc,
                VecDoub &scalvv,
                VecInt &indexvv,
                const int nbb,
                MatDoub &yy,
                Difeq &difeqq)
      : itmax(itmaxx)
      , conv(convv)
      , slowc(slowcc)
      , scalv(scalvv)
      , indexv(indexvv)
      , nb(nbb)
      , y(yy)
      , difeq(difeqq)
      , ne(y.rows())
      , m(y.cols())
      , kmax(ne)
      , ermax(ne)
      , c(ne, ne - nb + 1, m + 1)
      , s(MatDoub(ne, 2 * ne + 1, (double)0.))
  {
    int jv, k, nvars = ne * m;
    int k1 = 0, k2 = m;
    int j1 = 0, j2 = nb, j3 = nb, j4 = ne, j5 = j4 + j1, j6 = j4 + j2,
        j7 = j4 + j3, j8 = j4 + j4, j9 = j8 + j1;
    int ic1 = 0, ic2 = ne - nb, ic3 = ic2, ic4 = ne, jc1 = 0, jcf = ic3;
    std::cout << std::setw(8) << "Iter.";
    std::cout << std::setw(10) << "Error" << std::setw(10) << "FAC"
              << std::endl;
    for (int it = 0; it < itmax; it++)
    {
      k = k1;
      difeq.smatrix(k, k1, k2, j9, ic3, ic4, indexv, s, y);
      pinvs(ic3, ic4, j5, j9, jc1, k1);
      for (k = k1 + 1; k < k2; k++)
      {
        int kp = k;
        difeq.smatrix(k, k1, k2, j9, ic1, ic4, indexv, s, y);
        red(ic1, ic4, j1, j2, j3, j4, j9, ic3, jc1, jcf, kp);
        pinvs(ic1, ic4, j3, j9, jc1, k);
      }
      k = k2;
      difeq.smatrix(k, k1, k2, j9, ic1, ic2, indexv, s, y);
      red(ic1, ic2, j5, j6, j7, j8, j9, ic3, jc1, jcf, k2);
      pinvs(ic1, ic2, j7, j9, jcf, k2);
      bksub(jcf, k1, k2);
      double err = 0.0;
      for (int j = 0; j < ne; j++)
      {
        jv          = indexv[j];
        double errj = 0.0, vmax = 0.0;
        int km = 0;
        for (k = k1; k < k2; k++)
        {
          double vz = abs(c[jv][0][k]);
          if (vz > vmax)
          {
            vmax = vz;
            km   = k + 1;
          }
          errj += vz;
        }
        err += errj / scalv[j];
        ermax[j] = c[jv][0][km - 1] / scalv[j];
        kmax[j]  = km;
      }
      err /= nvars;
      double fac = (err > slowc ? slowc / err : 1.0);

      for (int j = 0; j < ne; j++)
      {
        jv = indexv[j];
        for (k = k1 + 1; k < k2 - 1; k++)
          y[j][k] -= fac * c[jv][0][k];
      }
      std::cout << std::setw(6) << it;
      std::cout << std::setw(13) << err;
      std::cout << std::setw(12) << fac << std::endl;

      if (err < conv and it > 3)
      {
        // Save individual iteration of the path
        std::string str = "output_" + std::to_string(1) + ".tsv";
        std::ofstream out(str);
        for (size_t i = 0; i < difeq.z.size(); i++)
        {
          out << difeq.z[i] << "\t" << y[0][i] << "\t" << y[1][i] << "\t"
              << y[2][i] << "\t" << y[3][i] << "\t" << y[4][i] << "\t"
              << y[5][i] << "\t" << y[6][i] << "\t" << y[7][i] << "\n";
        }
        out.close();

        return;
      }
    }
    // throw("Too many iterations in solvde");
  }

  /**
   * @brief pinvs diagonalizes the square subsection of s and stores unreduced
   * coefficients
   *
   * @param ie1
   * @param ie2
   * @param je1
   * @param jsf
   * @param jc1
   * @param k
   */
  inline void pinvs(const size_t ie1,
                    size_t ie2,
                    const size_t je1,
                    const size_t jsf,
                    const int jc1,
                    const int k)
  {
    int jpiv, jp, jcoff, irow, ipiv, icoff;
    size_t je2;
    double pivinv, piv, big;
    const size_t iesize = ie2 - ie1;
    VecInt indxr(iesize);
    VecDoub pscl(iesize);
    je2 = je1 + iesize;
    for (size_t i = ie1; i < ie2; i++)
    {
      big = 0.0;
      for (size_t j = je1; j < je2; j++)
        if (abs(s[i][j]) > big) big = abs(s[i][j]);
      if (big == 0.0) throw("Singular matrix - row all 0, in pinvs");
      pscl[i - ie1]  = 1.0 / big;
      indxr[i - ie1] = 0;
    }
    for (size_t id = 0; id < iesize; id++)
    {
      piv = 0.0;
      for (size_t i = ie1; i < ie2; i++)
      {
        if (indxr[i - ie1] == 0)
        {
          big = 0.0;
          for (size_t j = je1; j < je2; j++)
          {
            if (abs(s[i][j]) > big)
            {
              jp  = j;
              big = abs(s[i][j]);
            }
          }
          if (big * pscl[i - ie1] > piv)
          {
            ipiv = i;
            jpiv = jp;
            piv  = big * pscl[i - ie1];
          }
        }
      }
      if (s[ipiv][jpiv] == 0.0) throw("Singular matrix in routine pinvs");
      indxr[ipiv - ie1] = jpiv + 1;
      pivinv            = 1.0 / s[ipiv][jpiv];
      for (size_t j = je1; j <= jsf; j++)
        s[ipiv][j] *= pivinv;
      s[ipiv][jpiv] = 1.0;
      for (size_t i = ie1; i < ie2; i++)
      {
        if (indxr[i - ie1] != jpiv + 1)
        {
          if (s[i][jpiv] != 0.0)
          {
            double dum = s[i][jpiv];
            for (size_t j = je1; j <= jsf; j++)
              s[i][j] -= dum * s[ipiv][j];
            s[i][jpiv] = 0.0;
          }
        }
      }
    }
    jcoff = jc1 - je2;
    icoff = ie1 - je1;
    for (size_t i = ie1; i < ie2; i++)
    {
      irow = indxr[i - ie1] + icoff;
      for (size_t j = je2; j <= jsf; j++)
        c[irow - 1][j + jcoff][k] = s[i][j];
    }
  }

  inline void bksub(const int jf, const int k1, const int k2)
  {
    int nbf = ne - nb, im = 1;
    for (int k = k2 - 1; k >= k1; k--)
    {
      if (k == k1) im = nbf + 1;
      int kp = k + 1;
      for (int j = 0; j < nbf; j++)
      {
        double xx = c[j][jf][kp];
        for (int i = im - 1; i < ne; i++)
          c[i][jf][k] -= c[i][j][k] * xx;
      }
    }
    for (int k = k1; k < k2; k++)
    {
      int kp = k + 1;
      for (int i = 0; i < nb; i++)
        c[i][0][k] = c[i + nbf][jf][k];
      for (int i = 0; i < nbf; i++)
        c[i + nb][0][k] = c[i][jf][kp];
    }
  }

  /**
   * @brief  red eliminates leading columns of the s matrix using results from
   * prior blocks
   *
   * @param iz1
   * @param iz2
   * @param jz1
   * @param jz2
   * @param jm1
   * @param jm2
   * @param jmf
   * @param ic1
   * @param jc1
   * @param jcf
   * @param kc
   */
  inline void red(const int iz1,
                  const int iz2,
                  const int jz1,
                  const int jz2,
                  const int jm1,
                  const int jm2,
                  const int jmf,
                  const int ic1,
                  const int jc1,
                  const int jcf,
                  const int kc)
  {
    int l, j, i;
    double vx;
    int loff = jc1 - jm1, ic = ic1;
    for (j = jz1; j < jz2; j++)
    {
      for (l = jm1; l < jm2; l++)
      {
        vx = c[ic][l + loff][kc - 1];
        for (i = iz1; i < iz2; i++)
          s[i][l] -= s[i][j] * vx;
      }
      vx = c[ic][jcf][kc - 1];
      for (i = iz1; i < iz2; i++)
        s[i][jmf] -= s[i][j] * vx;
      ic += 1;
    }
  }
};