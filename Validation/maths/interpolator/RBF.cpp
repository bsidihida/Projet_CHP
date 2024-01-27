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
 * @file Interpolator.cpp
 * @brief This file contains the Interpolator class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-05-14
 * @copyright Inria
 */

#include <math.h>
#include "RBF.hpp"

namespace neos {

Rbf::Rbf() {
  _epsilon = .01;
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

void Rbf::buildMatrix(const std::vector<NPoint > &xNeighs) {
  int m_nNeighs = xNeighs.size();

  m_matrix.resize(m_nNeighs * m_nNeighs);
  for (int i = 0; i < m_nNeighs; ++i) {
    for (int j = 0; j <m_nNeighs; ++j) {
      m_matrix[i * m_nNeighs + j] = computeRBFValue(xNeighs[j], xNeighs[i]);
    }
  }
}

double Rbf::computeRBFValue(const NPoint &x0, const NPoint &x1) {
  NPoint x;
  double phi;

  for (int i = 0; i < 3; ++i)
    x[i] = x0[i] - x1[i];
  phi = exp(-_epsilon * _epsilon * norm2(x) * norm2(x));
  return phi;
}

std::vector<double> Rbf::computeInterpolationRBF(const NPoint &x0,
                                                 const std::vector<NPoint > &xNeighs) {
  std::vector<double> m_weights;
  int m_nNeighs = xNeighs.size();
  double m[m_nNeighs * m_nNeighs];

  m_matrix.clear();
  m_weights.resize(m_nNeighs,0.0);
  buildMatrix(xNeighs);
  for (int i = 0; i < m_nNeighs * m_nNeighs; ++i)
    m[i] = m_matrix[i];
  _solve(m_nNeighs, m);
  for (int i = 0; i < m_nNeighs * m_nNeighs; ++i)
    m_matrix[i] = m[i];
  for (int i = 0; i < m_nNeighs; ++i) {
    for (int j = 0; j < m_nNeighs; ++j) {
      m_weights[i] += m_matrix[j * m_nNeighs + i] * computeRBFValue(xNeighs[j],x0);
    }
  }
  return m_weights;
}

double Rbf::computeInterpolation(const NPoint &xpos,
                                 const std::vector<NPoint > &xref,
                                 const std::vector<double> &vref) {
  std::vector<double> weights;
  double val_int = 0.0;
  double s       = 0.0;

  weights = computeInterpolationRBF(xpos, xref);
  for (std::size_t i = 0; i<weights.size(); ++i)
    s += weights[i];
  for (std::size_t i = 0; i<weights.size(); ++i)
    weights[i] /= s;
  for (std::size_t i = 0; i < vref.size(); ++i) {
    val_int += weights[i] * vref[i];
  }
  return val_int;
}

std::vector<double> Rbf::computeInterpolation(const NPoint &xpos,
                                              const std::vector<NPoint > &xref) {
  std::vector<double> weights;
  double s = 0.0;

  weights = computeInterpolationRBF(xpos, xref);
  for (std::size_t i = 0; i<weights.size(); ++i)
    s += weights[i];
  for (std::size_t i = 0; i < weights.size(); ++i)
    weights[i] /= s;
  return weights;
}
}
