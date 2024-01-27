/* -------------------------------------------------------------------------*\
 *
 *  NEOS
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of Neos.
 *
 *  Neos is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  Neos is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Neos. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

/**
 * @file Polynomial3D.cpp
 * @brief This file contains the 3D Polynomial Interpolator class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-05-14
 * @copyright Inria
 */

#include <iostream>
#include <math.h>
#include <numeric>
#include <iomanip>
#include <limits>
#include "Polynomial3D.hpp"

namespace neos {

Polynomial3D::Polynomial3D() {
#if USE_EIGEN3
  _solve = fullPiv;
  _matVecProduct = eiMatVecProduct;
#else
  _matVecProduct = matVecProduct;
#endif /* USE_EIGEN3 */
#ifdef USE_LAPACKE
  _solve = inverse;
#endif /* USE_LAPACKE */
}

double Polynomial3D::computeInterpolation4P(const NPoint &xpos, const std::vector<std::array<double,3> > &xref, const std::vector<double> &val_ref) {
  double A[xref.size() * xref.size()];
  double b[xref.size()];
  int size = xref.size();

  std::copy(val_ref.begin(),val_ref.end(), b);

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()]     = xref[i][NPX];
    A[i * xref.size() + 1] = xref[i][NPY];
    A[i * xref.size() + 2] = xref[i][NPZ];
    A[i * xref.size() + 3] = 1.0;
  }
  _solve(size, A);
  _matVecProduct(xref.size(), A, b);

  return b[0] * xpos[NPX] + b[1] * xpos[NPY] + b[2] * xpos[NPZ] + b[3];
}

std::vector<double> Polynomial3D::computeInterpolation4P(const NPoint &xpos, const std::vector<NPoint > &xref) {
  double A[xref.size() * xref.size()];
  std::vector<double> res;

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()]     = xref[i][NPX];
    A[i * xref.size() + 1] = xref[i][NPY];
    A[i * xref.size() + 2] = xref[i][NPZ];
    A[i * xref.size() + 3] = 1.0;
  }
  _solve(xref.size(), A);
  size_t nPoint = xref.size();
  for (size_t i = 0; i < nPoint; i++)
    res.push_back(xpos[NPX] * A[i] + xpos[NPY] * A[i + nPoint] + xpos[NPZ] * A[i + 2*nPoint]  + A[i + 3*nPoint]);
  return res;
}

double Polynomial3D::computeInterpolation8P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref) {
  double A[xref.size() * xref.size()];
  double b[xref.size()];

  std::copy(val_ref.begin(),val_ref.end(), b);

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()]     = xref[i][NPX] * xref[i][NPY] * xref[i][NPZ];
    A[i * xref.size() + 1] = xref[i][NPX] * xref[i][NPY] *      1;
    A[i * xref.size() + 2] = xref[i][NPX] *      1     * xref[i][NPZ];
    A[i * xref.size() + 3] =      1     * xref[i][NPY] * xref[i][NPZ];
    A[i * xref.size() + 4] = xref[i][NPX] *      1     *      1;
    A[i * xref.size() + 5] =      1     * xref[i][NPY] *      1;
    A[i * xref.size() + 6] =      1     *      1     * xref[i][NPZ];
    A[i * xref.size() + 7] =      1     *      1     *      1;
  }
  _solve(xref.size(), A);
  _matVecProduct(xref.size(), A, b);

  return b[0] * xpos[NPX] * xpos[NPY] * xpos[NPZ] +
         b[1] * xpos[NPX] * xpos[NPY] +
         b[2] * xpos[NPX] * xpos[NPZ] +
         b[3] * xpos[NPY] * xpos[NPZ] +
         b[4] * xpos[NPX] +
         b[5] * xpos[NPY] +
         b[6] * xpos[NPZ] +
         b[7];
}

std::vector<double> Polynomial3D::computeInterpolation8P(const NPoint &xpos, const std::vector<NPoint > &xref) {
  double A[xref.size() * xref.size()];
  std::vector<double> res;

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()]     = xref[i][NPX] * xref[i][NPY] * xref[i][NPZ];
    A[i * xref.size() + 1] = xref[i][NPX] * xref[i][NPY] *      1;
    A[i * xref.size() + 2] = xref[i][NPX] *      1     * xref[i][NPZ];
    A[i * xref.size() + 3] =      1     * xref[i][NPY] * xref[i][NPZ];
    A[i * xref.size() + 4] = xref[i][NPX] *      1     *      1;
    A[i * xref.size() + 5] =      1     * xref[i][NPY] *      1;
    A[i * xref.size() + 6] =      1     *      1     * xref[i][NPZ];
    A[i * xref.size() + 7] =      1     *      1     *      1;
  }
  _solve(xref.size(), A);
  size_t nPoint = xref.size();
  for (size_t i = 0; i < nPoint; i++)
    res.push_back(xpos[NPX] * xpos[NPY] * xpos[NPZ] * A[i] + xpos[NPX] * xpos[NPY] * A[i + nPoint]  + xpos[NPX] * xpos[NPZ] * A[i + 2*nPoint] + xpos[NPY] * xpos[NPZ] * A[i + 3*nPoint] + xpos[NPX] * A[i + 4*nPoint] + xpos[NPY] * A[i + 5*nPoint] + xpos[NPZ] * A[i + 6*nPoint] + A[i + 7*nPoint]);
  return res;
}

double Polynomial3D::computeInterpolation7P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref) {
  double A[xref.size() * xref.size()];
  double b[xref.size()];

  std::copy(val_ref.begin(),val_ref.end(), b);

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()    ] = xref[i][NPX] * xref[i][NPY] *      1;
    A[i * xref.size() + 1] = xref[i][NPX] *      1     * xref[i][NPZ];
    A[i * xref.size() + 2] =      1     * xref[i][NPY] * xref[i][NPZ];
    A[i * xref.size() + 3] = xref[i][NPX] *      1     *      1;
    A[i * xref.size() + 4] =      1     * xref[i][NPY] *      1;
    A[i * xref.size() + 5] =      1     *      1     * xref[i][NPZ];
    A[i * xref.size() + 6] =      1     *      1     *      1;
  }
  _solve(xref.size(), A);
  _matVecProduct(xref.size(), A, b);

  return b[0] * xpos[NPX] * xpos[NPY] +
         b[1] * xpos[NPX] * xpos[NPZ] +
         b[2] * xpos[NPY] * xpos[NPZ] +
         b[3] * xpos[NPX] +
         b[4] * xpos[NPY] +
         b[5] * xpos[NPZ] +
         b[6];
}

std::vector<double> Polynomial3D::computeInterpolation7P(const NPoint &xpos, const std::vector<NPoint > &xref) {
  double A[xref.size() * xref.size()];
  std::vector<double> res;

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()] = xref[i][NPX] * xref[i][NPY] *      1;
    A[i * xref.size() + 1] = xref[i][NPX] *      1     * xref[i][NPZ];
    A[i * xref.size() + 2] =      1     * xref[i][NPY] * xref[i][NPZ];
    A[i * xref.size() + 3] = xref[i][NPX] *      1     *      1;
    A[i * xref.size() + 4] =      1     * xref[i][NPY] *      1;
    A[i * xref.size() + 5] =      1     *      1     * xref[i][NPZ];
    A[i * xref.size() + 6] =      1     *      1     *      1;
  }
  _solve(xref.size(), A);
  size_t nPoint = xref.size();
  for (size_t i = 0; i < nPoint; i++)
    res.push_back(xpos[NPX] * xpos[NPY] * A[i] + xpos[NPX] * xpos[NPZ] * A[i + nPoint]  + xpos[NPY] * xpos[NPZ] * A[i + 2*nPoint] + xpos[NPX] * A[i + 3*nPoint] + xpos[NPY] * A[i + 4*nPoint] + xpos[NPZ] * A[i + 5*nPoint] +A[i + 6*nPoint]);
  return res;
}

double Polynomial3D::computeInterpolation10P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref) {
  double A[xref.size() * xref.size()];
  double b[xref.size()];

  std::copy(val_ref.begin(),val_ref.end(), b);

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()     ] =     xref[i][NPX]    *      1            *      1;
    A[i * xref.size() + 1 ] =         1         * xref[i][NPY]        *      1;
    A[i * xref.size() + 2 ] =         1         *      1            * xref[i][NPZ];
    A[i * xref.size() + 3 ] = pow(xref[i][NPX],2) *      1            *      1;
    A[i * xref.size() + 4 ] =         1         * pow(xref[i][NPY],2) *      1;
    A[i * xref.size() + 5 ] =         1         *      1            * pow(xref[i][NPZ],2);
    A[i * xref.size() + 6 ] =     xref[i][NPX]    * xref[i][NPY]        *      1;
    A[i * xref.size() + 7 ] =     xref[i][NPX]    *      1            * xref[i][NPZ];
    A[i * xref.size() + 8 ] =         1         * xref[i][NPY]        * xref[i][NPZ];
    A[i * xref.size() + 9 ] =         1         *      1            *      1;
  }
  _solve(xref.size(), A);
  _matVecProduct(xref.size(), A, b);

  return b[0] * xpos[NPX] +
         b[1] * xpos[NPY] +
         b[2] * xpos[NPZ] +
         b[3] * pow(xpos[NPX],2) +
         b[4] * pow(xpos[NPY],2) +
         b[5] * pow(xpos[NPZ],2) +
         b[6] * xpos[NPX] * xpos[NPY] +
         b[7] * xpos[NPX] * xpos[NPZ] +
         b[8] * xpos[NPY] * xpos[NPZ] +
         b[9];
}

std::vector<double> Polynomial3D::computeInterpolation10P(const NPoint &xpos, const std::vector<NPoint > &xref) {
  double A[xref.size() * xref.size()];
  std::vector<double> res;

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()     ] =     xref[i][NPX]    *      1            *      1;
    A[i * xref.size() + 1 ] =         1         * xref[i][NPY]        *      1;
    A[i * xref.size() + 2 ] =         1         *      1            * xref[i][NPZ];
    A[i * xref.size() + 3 ] = pow(xref[i][NPX],2) *      1            *      1;
    A[i * xref.size() + 4 ] =         1         * pow(xref[i][NPY],2) *      1;
    A[i * xref.size() + 5 ] =         1         *      1            * pow(xref[i][NPZ],2);
    A[i * xref.size() + 6 ] =     xref[i][NPX]    * xref[i][NPY]        *      1;
    A[i * xref.size() + 7 ] =     xref[i][NPX]    *      1            * xref[i][NPZ];
    A[i * xref.size() + 8 ] =         1         * xref[i][NPY]        * xref[i][NPZ];
    A[i * xref.size() + 9 ] =         1         *      1            *      1;
  }
  _solve(xref.size(), A);
  size_t nPoint = xref.size();
  for (size_t i = 0; i < nPoint; i++)
    res.push_back(xpos[NPX] * A[i] + xpos[NPY] * A[i + nPoint] + xpos[NPZ] * A[i + 2*nPoint] + pow(xpos[NPX],2) * A[i + 3*nPoint] + pow(xpos[NPY],2) * A[i + 4*nPoint] + pow(xpos[NPZ],2) * A[i + 5*nPoint] + xpos[NPX] * xpos[NPY] * A[i + 6*nPoint] + xpos[NPX] * xpos[NPZ] * A[i + 7*nPoint] + xpos[NPY] * xpos[NPY] * A[i + 8*nPoint] + A[i + 9*nPoint]);
  return res;
}


double Polynomial3D::computeInterpolation11P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref) {
  double A[xref.size() * xref.size()];
  double b[xref.size()];
  int size = xref.size();
  std::copy(val_ref.begin(),val_ref.end(), b);

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()     ] =     xref[i][NPX]    *      1            *      1;
    A[i * xref.size() + 1 ] =         1         * xref[i][NPY]        *      1;
    A[i * xref.size() + 2 ] =         1         *      1            * xref[i][NPZ];
    A[i * xref.size() + 3 ] = pow(xref[i][NPX],2) *      1            *      1;
    A[i * xref.size() + 4 ] =         1         * pow(xref[i][NPY],2) *      1;
    A[i * xref.size() + 5 ] =         1         *      1            * pow(xref[i][NPZ],2);
    A[i * xref.size() + 6 ] =     xref[i][NPX]    * xref[i][NPY]        *      1;
    A[i * xref.size() + 7 ] =     xref[i][NPX]    *      1            * xref[i][NPZ];
    A[i * xref.size() + 8 ] =         1         * xref[i][NPY]        * xref[i][NPZ];
    A[i * xref.size() + 9 ] =     xref[i][NPX]    * xref[i][NPY]        * xref[i][NPZ];
    A[i * xref.size() + 10] =         1         *      1            *      1;
  }
  _solve(size, A);
  _matVecProduct(xref.size(), A, b);

  return b[0] * xpos[NPX] +
         b[1] * xpos[NPY] +
         b[2] * xpos[NPZ] +
         b[3] * pow(xpos[NPX],2) +
         b[4] * pow(xpos[NPY],2) +
         b[5] * pow(xpos[NPZ],2) +
         b[6] * xpos[NPX] * xpos[NPY] +
         b[7] * xpos[NPX] * xpos[NPZ] +
         b[8] * xpos[NPY] * xpos[NPZ] +
         b[9] * xpos[NPX] * xpos[NPY] * xpos[NPZ] +
         b[10];
}

std::vector<double> Polynomial3D::computeInterpolation11P(const NPoint &xpos, const std::vector<NPoint > &xref) {
  double A[xref.size() * xref.size()];
  std::vector<double> res;


  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()     ] =     xref[i][NPX]    *      1            *      1;
    A[i * xref.size() + 1 ] =         1         * xref[i][NPY]        *      1;
    A[i * xref.size() + 2 ] =         1         *      1            * xref[i][NPZ];
    A[i * xref.size() + 3 ] = pow(xref[i][NPX],2) *      1            *      1;
    A[i * xref.size() + 4 ] =         1         * pow(xref[i][NPY],2) *      1;
    A[i * xref.size() + 5 ] =         1         *      1            * pow(xref[i][NPZ],2);
    A[i * xref.size() + 6 ] =     xref[i][NPX]    * xref[i][NPY]        *      1;
    A[i * xref.size() + 7 ] =     xref[i][NPX]    *      1            * xref[i][NPZ];
    A[i * xref.size() + 8 ] =         1         * xref[i][NPY]        * xref[i][NPZ];
    A[i * xref.size() + 9 ] =     xref[i][NPX]    * xref[i][NPY]        * xref[i][NPZ];
    A[i * xref.size() + 10] =         1         *      1            *      1;
  }
  _solve(11, A);
  size_t nPoint = xref.size();
  for (size_t i = 0; i < nPoint; i++)
    res.push_back(xpos[NPX] * A[i]  + xpos[NPY] * A[i + nPoint] + xpos[NPZ] * A[i + 2*nPoint] + pow(xpos[NPX],2) * A[i + 3*nPoint] + pow(xpos[NPY],2) * A[i + 4*nPoint] + pow(xpos[NPZ],2) * A[i + 5*nPoint] + xpos[NPX] * xpos[NPY] * A[i + 6*nPoint] + xpos[NPX] * xpos[NPZ] * A[i + 7*nPoint] + xpos[NPY] * xpos[NPY] * A[i + 8*nPoint] + xpos[NPX] * xpos[NPY] * xpos[NPZ] * A[i + 9*nPoint]   + A[i + 10*nPoint]);
  return res;
}

double Polynomial3D::computeInterpolation(const NPoint &xpos,
                                          const std::vector<NPoint > &xref,
                                          const std::vector<double> &vref) {
  std::vector<double> weights;
  double val_int = 0.0;
  double s       = 0.0;

  if (xref.size() == 4)
    return computeInterpolation4P(xpos, xref, vref);
  else if (xref.size() == 7)
    return computeInterpolation7P(xpos, xref, vref);
  else if (xref.size() == 8)
    return computeInterpolation8P(xpos, xref, vref);
  else if (xref.size() == 10)
    return computeInterpolation10P(xpos, xref, vref);
  else if (xref.size() == 11)
    return computeInterpolation11P(xpos, xref, vref);
  else {
    std::cout << "Interpolation Polynomial3D not implemented for size " << xref.size() << std::endl;
    exit(-1);
  }

  for (std::size_t i = 0; i<weights.size(); ++i)
    s += weights[i];
  for (std::size_t i = 0; i<weights.size(); ++i)
    weights[i] /= s;
  for (std::size_t i = 0; i < vref.size(); ++i) {
    val_int += weights[i] * vref[i];
  }
  return val_int;
}

std::vector<double> Polynomial3D::computeInterpolation(const NPoint &xpos,
                                                       const std::vector<NPoint > &xref) {
  std::vector<double> weights;
  double s = 0.0;

  if (xref.size() == 4)
    return computeInterpolation4P(xpos, xref);
  else if (xref.size() == 7)
    return computeInterpolation7P(xpos, xref);
  else if (xref.size() == 8)
    return computeInterpolation8P(xpos, xref);
  else if (xref.size() == 10)
    return computeInterpolation10P(xpos, xref);
  else if (xref.size() == 11)
    return computeInterpolation11P(xpos, xref);
  else {
    std::cout << "Interpolation Polynomial3D not implemented for size " << xref.size() << std::endl;
    exit(-1);
  }

  for (std::size_t i = 0; i<weights.size(); ++i)
    s += weights[i];
  for (std::size_t i = 0; i < weights.size(); ++i)
    weights[i] /= s;
  return weights;
}
}
