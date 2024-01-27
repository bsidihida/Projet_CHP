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
 * @file Distance.cpp
 * @brief This file contains the Interpolator class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-05-14
 * @copyright Inria
 */

#include <math.h>
#include "Distance.hpp"

namespace neos {

Distance::Distance() {
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

std::vector<double> Distance::computeDistanceWeightedInterpolation(const NPoint &xpos,
                                                                   const std::vector<NPoint > &xref,
                                                                   int power)
{
  std::vector<double> weights;
  double s = 0;

  for (std::size_t i = 0; i < xref.size(); ++i) {
    double d = sqrt((xpos[NPX]-xref[i][NPX]) * (xpos[NPX]-xref[i][NPX]) +
                    (xpos[NPY]-xref[i][NPY]) * (xpos[NPY]-xref[i][NPY]) +
                    (xpos[NPZ]-xref[i][NPZ]) * (xpos[NPZ]-xref[i][NPZ]));
    weights.push_back(pow(1.0 / (d + 1e-16), power));
    s += weights[i];
  }
  s = 1./s;
  for (std::size_t i = 0; i < xref.size(); ++i)
    weights[i] *= s;
  return weights;
}

double Distance::computeInterpolation(const NPoint &xpos,
                                      const std::vector<NPoint > &xref,
                                      const std::vector<double> &vref) {
  std::vector<double> weights;
  double val_int = 0.0;
  double s       = 0.0;

  weights = computeDistanceWeightedInterpolation(xpos, xref, 3);
  for (std::size_t i = 0; i<weights.size(); ++i)
    s += weights[i];
  for (std::size_t i = 0; i<weights.size(); ++i)
    weights[i] /= s;
  for (std::size_t i = 0; i < vref.size(); ++i) {
    val_int += weights[i] * vref[i];
  }
  return val_int;
}

std::vector<double> Distance::computeInterpolation(const NPoint &xpos,
                                                   const std::vector<NPoint > &xref) {
  std::vector<double> weights;
  double s = 0.0;

  weights = computeDistanceWeightedInterpolation(xpos, xref, 3);
  for (std::size_t i = 0; i<weights.size(); ++i)
    s += weights[i];
  for (std::size_t i = 0; i < weights.size(); ++i)
    weights[i] /= s;
  return weights;
}
}
