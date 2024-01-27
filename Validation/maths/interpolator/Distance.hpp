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

#ifndef __DISTANCE_HPP__
#define __DISTANCE_HPP__

#include <iostream>
#include <vector>
#include <array>
#include <functional>
#include "common.hpp"
#include "Function.hpp"
#include "IInterpolator.hpp"

namespace neos {

/**
 * @class Distance
 * @brief DistanceWeigthed interpolation algorithm class
 */
class Distance : public IInterpolator {
private:
std::function<int(int, double*)> _solve;
std::function<int(int, double*, double*)> _matVecProduct;

/**
 * @brief Return interpolated weights
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] val_ref Neighbours values
 */
std::vector<double> computeDistanceWeightedInterpolation(const NPoint &xpos, const std::vector<NPoint > &xref, int power = 3);

public:
/**
 * @brief Interpolator class
 */
Distance();
~Distance() {
}

/**
 * @brief Return interpolated value
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] vref Neighbours values
 */
double computeInterpolation(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &vref);

/**
 * @brief Return interpolated weight
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 */
std::vector<double> computeInterpolation(const NPoint &xpos, const std::vector<NPoint > &xref);
};
}

#endif /* __DISTANCE_HPP__ */
