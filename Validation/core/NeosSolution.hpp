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

#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include "NeosPiercedVector.hpp"
#include "Elasticity.hpp"

namespace neos {

struct Solution
{
  ElasticStructures Elastic;
  PiercedVector<std::array<double,3> > VelocityCC;
  PiercedVector<double> Ux;
  PiercedVector<double> Uy;
  PiercedVector<double> Uz;
  PiercedVector<double> Pressure;
  PiercedVector<double> VelocityFC;
  PiercedVector<double> Vorticity;
  PiercedVector<double> FctInd;
  double t;
  double dt;
};
}
#endif
