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

#include "Neos.hpp"
#include "gtest/gtest.h"

#ifdef USE_LAPLACIAN
#include "petscksp.h"
#include "petscpc.h"
#include "petscvec.h"

using namespace neos;

static char help[] = "";

#endif /* USE_LAPLACIAN */

int main(int ac, char **av) {
#if defined (ENABLE_MPI)
  MPI_Init(&ac,&av);
#endif
#ifdef BITPIT_FOUND
  bitpit::log::manager().initialize(bitpit::log::COMBINED);
  bitpit::log::cout() << consoleVerbosity(bitpit::log::QUIET);
#endif
  int result = 0;
#ifdef USE_LAPLACIAN
  PetscInitialize(&ac, &av, (char *)0, help);
#endif /* USE_LAPLACIAN */
  ::testing::InitGoogleTest(&ac, av);
  result = RUN_ALL_TESTS();

#ifdef USE_LAPLACIAN
  PetscFinalize();
#endif /* USE_LAPLACIAN */
#if defined (ENABLE_MPI)
  MPI_Finalize();
#endif
  return result;
}
