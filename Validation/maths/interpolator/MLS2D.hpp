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
 * @file   MLS2D.hpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Tue Aug  6 11:22:17 2019
 *
 * @brief  This file contains Compact Interpolator Class
 *
 *
 */

#ifndef __MLS2D_HPP__
#define __MLS2D_HPP__

#include "IInterpolator.hpp"
#include "Function.hpp"

namespace neos {

class MLS2D : public IInterpolator
{

private:
std::function<int(int, double*)> _solve;
std::function<int(int, double*, double*)> _matVecProduct;

public:
/**
 * @brief constructor
 */
MLS2D();

/**
 * @brief destructor
 */
~MLS2D() {
}

/**
 * @brief Return interpolated weight
 *
 * @param[in] xpos Point position
 * @param[in] xref Neighbours positions
 * @param[in] type Interpolation type
 */
std::vector<double> computeInterpolation(const NPoint &xpos,
                                         const std::vector<NPoint> &xref);

double computeInterpolation(const NPoint& xpos,
                            const std::vector<NPoint> &xref,
                            const std::vector<double> &vref);

std::vector<double> computeFirstDerivatives(const NPoint &xpos,
                                            const std::vector<NPoint > &xref,
                                            const std::vector<double>
                                            &vref);
std::vector<double> computeSecondDerivatives(const NPoint &xpos,
                                             const std::vector<NPoint > &xref,
                                             const std::vector<double> &vref);

std::vector<double> computeInterpolationAtFaceCenter(const NPoint &xpos,
                                                     const std::vector<NPoint> &xref);


private:
std::vector<double> computeInterpolation2ndOrder(const NPoint &xpos,
                                                 const std::vector<NPoint>&
                                                 xref);

};
}


#endif
