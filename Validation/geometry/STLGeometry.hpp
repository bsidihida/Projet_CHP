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

#ifndef __STLGEOMETRY_HPP__
#define __STLGEOMETRY_HPP__

/**
 * @file STLGeometry.hpp
 * @brief This file contains Geometry class declaration
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-03-14
 * @copyright Inria
 */

#include <iostream>
#include <vector>
#include <functional>

// #include "bitpit.hpp"
#include <levelSet.hpp>
#include <Geometry.hpp>
#include <Grid.hpp>


/**
 * @class STLGeometry
 * @brief STLGeometry class
 *
 * Class to load stl files
 */

namespace neos {

class STLGeometry : public Geometry {
private:
Grid        *_grid;
bitpit::LevelSet _lvst;
int _stlID = -1;
PiercedVector<double>  _phi;
void updatePhi();

public:
STLGeometry() {};
STLGeometry(const std::string &file, Grid *grid);
~STLGeometry();
void addSTL(const std::string &file, Grid *grid);
double getLevelSet(const std::array<double, 3> &point);
double getLevelSet(const long id) {
  return _phi[id];
}

PiercedVector<double> &getPhi() {
  return _phi;
}
void computeLevelSet(Grid * grid);
};
}
#endif /* __STLGEOMETRY_HPP__ */
