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
 * @file Polynomial2D.cpp
 * @brief This file contains the 2D polynomial Interpolator class
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
#include "Polynomial2D.hpp"

namespace neos {

Polynomial2D::Polynomial2D() {
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


double Polynomial2D::computeInterpolation3P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref) {
  int size = xref.size();
  double A[xref.size() * xref.size()];
  double b[xref.size()];

  std::copy(val_ref.begin(),val_ref.end(), b);
  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()]     = xref[i][NPX];
    A[i * xref.size() + 1] = xref[i][NPY];
    A[i * xref.size() + 2] = 1.0;
  }
  _solve(size, A);
  _matVecProduct(xref.size(), A, b);
  return b[0] * xpos[NPX] + b[1] * xpos[NPY] + b[2];
}

std::vector<double> Polynomial2D::computeInterpolation3P(const NPoint &xpos, const std::vector<NPoint > &xref) {
  double A[xref.size() * xref.size()];
  std::vector<double> res;
  int size = xref.size();
  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()]     = xref[i][NPX];
    A[i * xref.size() + 1] = xref[i][NPY];
    A[i * xref.size() + 2] = 1.0;
  }
  _solve(size, A);

  size_t nPoint = xref.size();
  for (size_t i = 0; i < nPoint; i++)
    res.push_back(xpos[NPX] * A[i] + xpos[NPY] * A[i + nPoint] + A[i + 2*nPoint]);
  return res;
}

double Polynomial2D::computeInterpolation4P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref) {
  double A[xref.size() * xref.size()];
  double b[xref.size()];

  std::copy(val_ref.begin(),val_ref.end(), b);

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()]     = xref[i][NPX] * xref[i][NPY];
    A[i * xref.size() + 1] = xref[i][NPX];
    A[i * xref.size() + 2] = xref[i][NPY];
    A[i * xref.size() + 3] = 1.0;
  }
  _solve(xref.size(), A);
  _matVecProduct(xref.size(), A, b);
  return b[0] * xpos[NPX] * xpos[NPY] + b[1] * xpos[NPX] + b[2] * xpos[NPY] + b[3];
}

std::vector<double> Polynomial2D::computeInterpolation4P(const NPoint &xpos, const std::vector<NPoint > &xref) {
  double A[xref.size() * xref.size()];
  std::vector<double> res;

  for (std::size_t i = 0; i < xref.size(); ++i) {
    A[i * xref.size()]     = xref[i][NPX] * xref[i][NPY];
    A[i * xref.size() + 1] = xref[i][NPX];
    A[i * xref.size() + 2] = xref[i][NPY];
    A[i * xref.size() + 3] = 1.0;
  }
  _solve(xref.size(), A);

  size_t nPoint = xref.size();
  for (size_t i = 0; i < nPoint; i++)
    res.push_back(xpos[NPX] * xpos[NPY] * A[i] + xpos[NPX] * A[i + nPoint] + xpos[NPY] * A[i + 2*nPoint]  + A[i + 3*nPoint]);
  return res;
}


double Polynomial2D::computeInterpolation(const NPoint &xpos,
                                          const std::vector<NPoint > &xref,
                                          const std::vector<double> &vref) {
  std::vector<double> weights;
  double val_int = 0.0;
  double s       = 0.0;

  if (xref.size() == 3)
    return computeInterpolation3P(xpos, xref, vref);
  else if (xref.size() == 4)
    return computeInterpolation4P(xpos, xref, vref);
  else {
    std::cout << "Interpolation Polynomial2D not implemented for size " << xref.size() << std::endl;
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

std::vector<double> Polynomial2D::computeInterpolation(const NPoint &xpos,
                                                       const std::vector<NPoint > &xref) {
  std::vector<double> weights;
  double s = 0.0;

  if (xref.size() == 3)
    return computeInterpolation3P(xpos, xref);
  else if (xref.size() == 4)
    return computeInterpolation4P(xpos, xref);
  else {
    std::cout << "Interpolation Polynomial2D not implemented for size " << xref.size() << std::endl;
    exit(-1);
  }

  for (std::size_t i = 0; i<weights.size(); ++i)
    s += weights[i];
  for (std::size_t i = 0; i < weights.size(); ++i)
    weights[i] /= s;
  return weights;
}
}
