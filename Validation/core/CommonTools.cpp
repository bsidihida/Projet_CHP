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
 * @file CommonTools.cpp
 *
 * @brief This file contains CommonTools namespace
 *
 * @author Antoine Gerard
 * @date 2019-08-02
 *
 */

#include "CommonTools.hpp"

namespace neos
{

namespace CommonTools
{
/**
 * @brief flatten function for vector 2D
 *
 * @param v vector to flatten
 *
 * @return 1D representation of v
 */
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T> >& v)
{
  std::size_t total_size = 0;
  for (const auto& sub : v)
    total_size += sub.size();
  std::vector<T> result;
  result.reserve(total_size);
  for (const auto& sub : v)
    result.insert(result.end(), sub.begin(), sub.end());
  return result;
}

template <typename T>
T getGlobalOp(Grid* grid, const std::vector<T>& v, MPI_Op op)
{
  T globOpLoc;

  if (op == MPI_MAX)
  {
    globOpLoc = *std::max_element(v.begin(), v.end());
  }
  else if (op == MPI_MIN)
  {
    globOpLoc = *std::min_element(v.begin(), v.end());
  }

  T globOp(globOpLoc);
  if (grid->getProcessorCount()>1)
  {
    MPI_Barrier(grid->getCommunicator());
    MPI_Allreduce(&globOpLoc,
                  &globOp,
                  1,
                  MPI_DOUBLE,
                  op,
                  grid->getCommunicator());
  }
  return globOp;
}

/**
 * @brief Function to compute variable time step
 *
 * @param grid: grid reference to the mesh
 * @param param: pointers to a SimulationParameters object
 * @param Velocity: vector of velocities
 *
 * @return computed time step for next step
 */
double getTimeStep(Grid* grid,
                   SimulationParameters* param,
                   const double &t,
                   const std::vector<std::vector<double> > Velocity)
{
  double VelMax(0.), dx(0);
  auto aggregVel = CommonTools::flatten(Velocity);

  // Take absolute value of all value in vector
  for (double &v: aggregVel)
  {
    v = std::fabs(v);
  }

  VelMax = CommonTools::getGlobalOp(grid, aggregVel, MPI_MAX);
  VelMax = std::max(VelMax, 1.);

  dx = grid->getMinSize();
  double c = std::sqrt(param->getShearModulus() / param->getSolidDensity());

  double minDt = std::min(param->getCFL() * dx / c, param->getCFL() * dx / VelMax);

  return std::min(minDt, param->getTmax() - t + 1e-12);
}

void update(Solution& SolPrev, Solution& Sol)
{
  //SolPrev.t = Sol.t;
  //SolPrev.dt = Sol.dt;
  //SolPrev.Vorticity = Sol.Vorticity;
  //SolPrev.Pressure = Sol.Pressure;
  //SolPrev.Ux = Sol.Ux;
  //SolPrev.Uy = Sol.Uy;
  //SolPrev.Uz = Sol.Uz;
  //SolPrev.VelocityCC = Sol.VelocityCC;
  //SolPrev.Elastic = Sol.Elastic;
  //SolPrev.FctInd = Sol.FctInd;
  SolPrev = Sol;
}

void update(Solution& SolPrev, Solution& Sol, Solution& SolNext)
{
  update(SolPrev, Sol);
  update(Sol, SolNext);
}
}
}
