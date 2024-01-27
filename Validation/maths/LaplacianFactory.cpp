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
 * @file LaplacianFactory.cpp
 * @brief This file contain the source code for the LaplacianFactory class
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-11-23
 * @copyright Inria
 */

#include "solvers/LaplacianPetsc.hpp"
#include "LaplacianFactory.hpp"
//#include "AliceLaplacian.hpp"

namespace neos {

/**
 * The good laplacian class will be instancied and
 * a Laplacian will be returned
 */
Laplacian *LaplacianFactory::get(lapType ltype, solverType solver,
                                  Grid *grid,
                                  int argc, char **argv)
{
  Laplacian *lap = NULL;
  switch (ltype)
  {
  case lapType::FINITEVOLUME:
    lap =  new LaplacianPetsc(grid, true, argc,argv);
    break;
  case lapType::FINITEDIFF:
    lap =  new LaplacianPetsc(grid, true, argc,argv);
    break;
  default:
    lap =  new LaplacianPetsc(grid, true, argc,argv);
    break;
  }
  return lap;
}

///**
// * The good laplacian class will be instancied and
// * a Laplacian will be returned
// */
//Laplacian *LaplacianFactory::get(lapType ltype, solverType solver, Grid *grid, int argc, char **argv)
//{
//  Laplacian *lap = NULL;
//  switch (ltype)
//  {
//  case lapType::FINITEVOLUME:
//    lap =  new Laplacian(grid,argc,argv);
//    break;
//  case lapType::FINITEDIFF:
//    lap =  new Laplacian(grid,argc,argv);
//    break;
//  default:
//    lap =  new Laplacian(grid,argc,argv);
//    break;
//  }
//  return lap;
//}

/**
 * The good laplacian class will be instancied and
 * a Laplacian will be returned
 */
Laplacian *LaplacianFactory::get(lapType ltype, solverType solver, Grid *grid)
{
  Laplacian *lap = NULL;
  switch (ltype)
  {
  case lapType::FINITEVOLUME:
    lap =  new LaplacianPetsc(grid,true, 0,NULL);
    break;
  case lapType::FINITEDIFF:
    lap =  new LaplacianPetsc(grid,true, 0, NULL);
    break;
  default:
    lap =  new LaplacianPetsc(grid,true, 0,NULL);
    break;
  }
  return lap;
}
}
