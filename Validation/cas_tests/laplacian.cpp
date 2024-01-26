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

#include "bitpit.hpp"
#include "LaplacianFactory.hpp"
#include "Grid.hpp"
#include "gtest/gtest.h"
#include <math.h>

using namespace neos;

TEST(Laplacian, lpn) {
  int level = 4;
  Grid *g = new Grid(0, 0, 0, 1.0, 1.0/pow(2,level), 2);
  int test_succeed = 1;
  std::vector<double> error(3,0.0);
  std::vector<std::vector<double> > vec_err;
  Laplacian *lap = LaplacianFactory::get(lapType::FINITEVOLUME, solverType::PETSC, g);

  lap->setBCPositionCoef(BCPosition::OnInterface);
  lap->setRHS();
  PiercedVector<double> sol = lap->solve();
  std::vector<double> error_rank = lap->error();

  std::vector<double> T(g->getCellCount());
  std::vector<double> err(g->getCellCount());
  std::vector<double>::iterator iT=T.begin();
  std::vector<double>::iterator ierr=err.begin();

  for ( auto cell=g->internalCellBegin(); cell!=g->internalCellEnd(); ++cell) {
    const long &id=cell.getId();
    *iT=sol[id];
    ++iT;
    std::array<double,3> c=g->evalCellCentroid(id);
    *ierr=abs(sol[id]-sinh(M_PI*(1-c[0]))*sin(M_PI*c[1])/sinh(M_PI));
    ++ierr;
  }

  MPI_Allreduce(&error_rank[0],&error[0],1,MPI_DOUBLE,MPI_SUM,g->getCommunicator());
  MPI_Allreduce(&error_rank[1],&error[1],1,MPI_DOUBLE,MPI_SUM,g->getCommunicator());
  MPI_Allreduce(&error_rank[2],&error[2],1,MPI_DOUBLE,MPI_MAX,g->getCommunicator());
  error[1]=sqrt(error[1]);
  vec_err.push_back(error);

  if (g->getRank()==0)
    std::cout<<error<<std::endl;

  MPI_Barrier(g->getCommunicator());
  std::vector<std::vector<double> > ref_err(5,std::vector<double>(3,0.0));

  ref_err[0]={{ 0.000498261, 0.000990569, 0.00402514 }};
  double s=0;
  for (std::size_t i=0; i<1; ++i) {
    s+=abs(ref_err[i][0]-vec_err[i][0]);
    s+=abs(ref_err[i][1]-vec_err[i][1]);
    s+=abs(ref_err[i][2]-vec_err[i][2]);
  }

  if (s < 1e-8)
    test_succeed = 0;

  if (g->getRank()==0)
    std::cout<<"Error with respect to reference solution (sum of L1, L2, Linf errors on all grids) = "<<s<<std::endl;

  EXPECT_EQ(0, test_succeed);
  delete g;
}
