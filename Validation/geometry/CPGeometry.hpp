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

#ifndef __CPGEOMETRY_HPP__
#define __CPGEOMETRY_HPP__

#include "NGeometry.hpp"

/**
 * @file CPGeometry.hpp
 * @brief This file contains Cloud Point Geometry class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

namespace neos {

class CPGeometry : public NGeometry {
protected:
std::vector<std::array<double, 3> >   _bdr;     /**< Vector of boundary point */
int _precision;                                 /**< Number of boundary points */
int _slice;                                     /**< Number of slices (3D) */

public:
/**
 * @brief CPGeometry class constructor
 *
 * @param[in] x The X initial position
 * @param[in] y The Y initial position
 * @param[in] z The Z initial position
 * @param[in] radius Radius size
 * @param[in] precision Number of boundary points
 * @param[in] slice Number of slices
 */
CPGeometry(double x, double y, double z, double radius, int precision, int slice)
  : NGeometry(x, y, z, radius)
{
  _precision = precision;
  _slice     = slice;
}

/**
 * @brief CPGeometry class constructor
 *
 * @param[in] x The X initial position
 * @param[in] y The Y initial position
 * @param[in] radius Radius size
 * @param[in] precision Number of boundary points
 */
CPGeometry(double x, double y, double radius, int precision)
  : NGeometry(x, y, radius)
{
  _precision = precision;
  _slice     = 0;
}

/**
 * @brief Geometry class destructor
 */
virtual ~CPGeometry() {
}

/**
 * @brief precision getter
 *
 * @return _precision value
 */
int getPrecision() {
  return _precision;
}

/**
 * @brief precision getter
 *
 * @return _precision value
 */
int getSlice() {
  return _slice;
}

/**
 * @brief [] operator
 *
 * @return The element at the specified position in the vector _bdr.
 */
std::array<double, 3> operator[](const size_t i) {
  if (i >= _bdr.size()) {
    std::cout << "Index out of range" << std::endl;
    return { };
  }
  return _bdr[i];
}

/**
 * @brief callback for deformation
 */
void deformation(std::function<void(NGeometry*,std::vector<std::array<double, 3> >&)> Fct) {
  Fct(this, _bdr);
}

virtual std::array<double, 3> getPoint(const size_t i) = 0;
/**
 * @brief Compute boundary points
 *
 * Compute _precision boundary points.
 * Must be implemented in the childs classes
 */
virtual void compute() = 0;

/**
 * @brief Return the levelset value of the point
 *
 * @param[in] point coordinate array of the point { { x, y, z } }
 */
virtual double getLevelSet(const std::array<double, 3> &point) = 0;

/**
 * @brief Return the levelset
 *
 * @param[in] slice Number of slices
 */
// virtual std::vector<double> getLevelSet() = 0;
};
}
#endif /* __CPGEOMETRY_HPP__ */
