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

#ifndef __IINTERPOLATOR_HPP__
#define __IINTERPOLATOR_HPP__
/**
 * @file IInterpolator.hpp
 *
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-07-14
 * @copyright Inria
 */

#include <functional>
#include <vector>
#include "common.hpp"

namespace neos {

class IInterpolator {
public:
IInterpolator() {
}
virtual ~IInterpolator() {
}
/**
 * @brief Return interpolated value
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] vref Neighbours values
 */
virtual double computeInterpolation(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &vref) = 0;

/**
 * @brief Return interpolated weight
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 */
virtual std::vector<double> computeInterpolation(const NPoint &xpos, const std::vector<NPoint > &xref) = 0;
};
}
#endif /* __IINTERPOLATOR_HPP__ */
