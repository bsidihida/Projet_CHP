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
 * @file Polynomial2D.hpp
 * @brief This file contains the Polynomial Interpolator class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-05-14
 * @copyright Inria
 */

#ifndef __POLYNOMIAL2D_HPP__
#define __POLYNOMIAL2D_HPP__

#include "IInterpolator.hpp"
#include "Function.hpp"

namespace neos {

class Polynomial2D : public IInterpolator {
private:
std::function<int(int, double*)> _solve;
std::function<int(int, double*, double*)> _matVecProduct;

public:
/**
 * @brief constructor
 */
Polynomial2D();

/**
 * @brief destructor
 */
~Polynomial2D() {
}

/**
 * @brief Return interpolated value
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] vref Neighbours values
 * @param[i] type Interpolation type
 */
double computeInterpolation(const NPoint &xpos, const std::vector<std::array<double,3> > &xref, const std::vector<double> &vref);

/**
 * @brief Return interpolated weight
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] type Interpolation type
 */
std::vector<double> computeInterpolation(const NPoint &xpos, const std::vector<std::array<double,3> > &xref);

private:
/**
 * @brief Return polynomial interpolated value with 3points (2D)
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] val_ref Neighbours values
 */
double computeInterpolation3P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref);

/**
 * @brief Return polynomial interpolated value with 3points (2D)
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 */
std::vector<double> computeInterpolation3P(const NPoint &xpos, const std::vector<NPoint > &xref);

/**
 * @brief Return polynomial interpolated value with 4points (2D)
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] val_ref Neighbours values
 */
double computeInterpolation4P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref);

/**
 * @brief Return polynomial interpolated value with 4points (2D)
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 */
std::vector<double> computeInterpolation4P(const NPoint &xpos, const std::vector<NPoint > &xref);
};
}

#endif /* __POLYNOMIAL2D_HPP__ */
