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
 * @file Polynomial3D.hpp
 * @brief This file contains the Polynomial Interpolator class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-05-14
 * @copyright Inria
 */

#ifndef __POLYNOMIAL3D_HPP__
#define __POLYNOMIAL3D_HPP__

#include <vector>
#include "common.hpp"
#include "IInterpolator.hpp"
#include "Function.hpp"

namespace neos {

class Polynomial3D : public IInterpolator {
private:
std::function<int(int, double*)> _solve;
std::function<int(int, double*, double*)> _matVecProduct;

public:
/**
 * @brief constructor
 */
Polynomial3D();

/**
 * @brief destructor
 */
~Polynomial3D() {
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

private:
/**
 * @brief Return polynomial interpolated value with 3points (3D)
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] val_ref Neighbours values
 */
double computeInterpolation4P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref);
std::vector<double> computeInterpolation4P(const NPoint &xpos, const std::vector<NPoint > &xref);

/**
 * @brief Return polynomial interpolated value with 7points (3D)
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] val_ref Neighbours values
 */
double computeInterpolation7P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref);
std::vector<double> computeInterpolation7P(const NPoint &xpos, const std::vector<NPoint > &xref);

/**
 * @brief Return polynomial interpolated value with 8points (3D)
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] val_ref Neighbours values
 */
double computeInterpolation8P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref);
std::vector<double> computeInterpolation8P(const NPoint &xpos, const std::vector<NPoint > &xref);

/**
 * @brief Return polynomial interpolated value with 13points (3D)
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] val_ref Neighbours values
 */
double computeInterpolation10P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref);
std::vector<double> computeInterpolation10P(const NPoint &xpos, const std::vector<NPoint > &xref);

/**
 * @brief Return polynomial interpolated value with 13points (3D)
 *
 * @param[i] xpos Point position
 * @param[i] xref Neighbours positions
 * @param[i] val_ref Neighbours values
 */
double computeInterpolation11P(const NPoint &xpos, const std::vector<NPoint > &xref, const std::vector<double> &val_ref);
std::vector<double> computeInterpolation11P(const NPoint &xpos, const std::vector<NPoint > &xref);
};
}

#endif /* __POLYNOMIAL3D_HPP__ */
