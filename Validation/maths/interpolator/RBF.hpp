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
 * @file Interpolator.hpp
 * @brief This file contains the Interpolator class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-05-14
 * @copyright Inria
 */

#ifndef __RBF_HPP__
#define __RBF_HPP__

#include <vector>
#include "common.hpp"
#include "Function.hpp"
#include "IInterpolator.hpp"

namespace neos {

/**
 * @class Interpolator
 * @brief Interpolator class
 */
class Rbf : public IInterpolator {
private:
double _epsilon;                /**< Epsilon. For BRF interpolation */
std::vector<double> m_matrix;   /**< matrix. For BRF interpolation */
std::function<int(int, double*)> _solve;
std::function<int(int, double*, double*)> _matVecProduct;

void buildMatrix(const std::vector<NPoint> &xNeighs);
double computeRBFValue(const NPoint &x0, const NPoint &x1);

/**
 * @brief Return Rbf interpolated weights
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 */
std::vector<double> computeInterpolationRBF(const NPoint &x0, const std::vector<NPoint > &xref);

public:
/**
 * @brief Interpolator class
 */
Rbf();
~Rbf() {
}

/**
 * @brief Set epsilon value
 *
 * @param[i] epsilon Epsilon value
 */
void setEpsilon(double epsilon) {
  _epsilon = epsilon;
}

/**
 * @brief Get epsilon value
 *
 * @return Epsilon value
 */
double getEpsilon() {
  return _epsilon;
}

/**
 * @brief Return interpolated value
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] vref Neighbours values
 * @param[i] type Interpolation type
 */
double computeInterpolation(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &vref);

/**
 * @brief Return interpolated weight
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] type Interpolation type
 */
std::vector<double> computeInterpolation(const NPoint &xpos, const std::vector<NPoint > &xref);
};
}
#endif /* __Rbf_HPP__ */
