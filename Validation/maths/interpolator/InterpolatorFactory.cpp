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
 * @file InterpolationFactory.cpp
 * @brief This file contain the source code for the InterpolationFactory class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-11-23
 * @copyright Inria
 */

#include "InterpolatorFactory.hpp"
#include "Polynomial2D.hpp"
#include "Polynomial3D.hpp"
#include "RBF.hpp"
#include "MLS2D.hpp"
#include "Distance.hpp"

namespace neos {

/**
 * The good interpolator class will be instancied and
 * a IInterpolator will be returned
 */
IInterpolator *InterpolatorFactory::get(interpoType itype=interpoType::DISTANCEWEIGHTED)
{
  IInterpolator *interpo = NULL;
  switch (itype)
  {
  case interpoType::POLY2D:
    interpo = new Polynomial2D();
    break;
  case interpoType::POLY3D:
    interpo = new Polynomial3D();
    break;
  case interpoType::RBF:
    interpo= new Rbf();
    break;
  case interpoType::MLS:
    interpo= new MLS2D();
    break;
  default:
    interpo = new Distance();
    break;
  }
  return interpo;
}
}
