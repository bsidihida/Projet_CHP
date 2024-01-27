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
 * @file Laplacian.cpp
 * @brief This file contains florian laplacian algo
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-12-05
 * @copyright Inria
 */

#include "Laplacian.hpp"
#include "Utils.hpp"
#include <VTK.hpp>
#include <numeric>
#include "InterpolatorFactory.hpp"

namespace neos {

UserDataComm<double>* Laplacian::_userComm;
int Laplacian::_userComm_only_one_init_dirty=0;

Laplacian::Laplacian(Grid *grid_, bool buildMatrices_)
  : _grid(grid_),
  _mapGlobal(_grid->getCellGlobalMap()),
  _positionCoef(1),
  _gradient(NULL),
  _hasStencil(false),
  _neumann(true),
  _buildMatrices(buildMatrices_),
  _lambda({}),
  _lambdaf(NULL)
{}

Laplacian::Laplacian(Grid* grid, StencilBuilder* stencils)
  : _grid(grid),
  _mapGlobal(_grid->getCellGlobalMap()),
  _positionCoef(1),
  _gradient(NULL),
  _stencils(stencils),
  _hasStencil(stencils != NULL),
  _neumann(true),
  _buildMatrices(true),
  _lambda({}),
  _lambdaf(NULL)
{}

void Laplacian::setStencils(StencilBuilder* stencils)
{
  _stencils = stencils;
  _hasStencil = stencils != NULL;
}

void Laplacian::setLambda(std::vector<double>* lambda)
{
  _lambda = VtoPV(*lambda, _grid);

  if (Laplacian::_userComm_only_one_init_dirty!=10101)
  {
    Laplacian::_userComm = new UserDataComm<double>(*_grid, 1);
    Laplacian::_userComm_only_one_init_dirty=10101;
  }

  Laplacian::_userComm->update(_lambda.size());
  Laplacian::_userComm->communicate(_lambda);

  _lambdaf = NULL;
}

void Laplacian::setLambda(double (*lambdaf)(double,double,double))
{
  _lambdaf = lambdaf;
  if (_lambda.size()) _lambda.clear();
}

void Laplacian::toggleNeumann(bool neumann) {
  _neumann = neumann;
}

void Laplacian::putRHS(const PiercedVector<double> &rhs, RHSOpType operation)
{
  for (auto &cell : _grid->getCells())
  {
    const long &id = cell.getId();
    int gl = _mapGlobal.at(id);
    if (cell.isInterior())
    {
      if (operation == RHSOpType::ADD) {
        addRHSValue(gl, rhs[id]);
      } else {
        setRHSValue(gl, rhs[id]);
      }
    }
  }
}

void Laplacian::setRHS(const PiercedVector<double> &rhs)
{
  putRHS(rhs, RHSOpType::SET);
}

void Laplacian::setRHS(const std::vector<double> &rhs)
{
  PiercedVector<double> val = VtoPV(rhs, _grid);

  val.checkIntegrity();

  putRHS(val, RHSOpType::SET);
}

void Laplacian::addRHS(const PiercedVector<double> &rhs)
{
  putRHS(rhs, RHSOpType::ADD);
}

void Laplacian::addRHS(const std::vector<double> &rhs)
{
  PiercedVector<double> val = VtoPV(rhs, _grid);

  val.checkIntegrity();

  putRHS(val, RHSOpType::ADD);
}

void Laplacian::setRHS(double (*callback)(int, double))
{
  PiercedVector<double> sol_pv = getSolutionPV();
  for (auto &cell : _grid->getCells())
  {
    const long &id = cell.getId();
    int gl = _mapGlobal.at(id);
    if (cell.isInterior())
    {
      double sol[1];
      sol[0] = callback(id, sol_pv[gl]);
      setRHSValue(gl, sol[0]);
    }
  }
}

void Laplacian::penalize(const std::vector<double> &phi)
{
  PiercedVector<double> val = VtoPV(phi, _grid);

  for (auto &cell : _grid->getCells())
  {
    const long &id = cell.getId();
    int gl = _mapGlobal.at(id);
    if (cell.isInterior() && val[id]!=0.0)
    {
      double t_val = val[id];
      addMatrixValue(gl, gl, t_val);
    }
  }
}

void Laplacian::penalizeAtOrder1(const PiercedVector<double> &fctInd,
                                 double coeffPen)
{
  for (auto &cell : _grid->getCells())
  {
    const long &id = cell.getId();
    int gl = _mapGlobal.at(id);
    if (cell.isInterior() && fctInd[id] < 0.0)
    {
      addMatrixValue(gl, gl, coeffPen);
    }
  }
}

void Laplacian::computeRHSHeat2D()
{
  std::array<double, 3> ic;
  double vp;

  for (auto &inter : _grid->getInterfaces())
  {
    const long &id=inter.getId();
    if (inter.isBorder())
    {
      int gl=_mapGlobal.at(inter.getOwner());
      ic=_grid->evalInterfaceCentroid(id);
      vp=-sinh(M_PI*(1-ic[0]))*sin(M_PI*ic[1])/sinh(M_PI)*_grid->evalInterfaceArea(id)/
          (_grid->evalCellSize(inter.getOwner())/2.0)/_grid->evalCellVolume(inter.getOwner());
      addRHSValue(gl, vp);
    }
  }
}

PiercedVector<double> Laplacian::getSolutionPV() {
  return VtoPV(getSolution(),_grid);
}

void Laplacian::buildMatrix()
{
  double vp, vm;

  if (_gradient == NULL) {
    _gradient = new Gradient(_grid->getDimension(), _grid);
  }

  for (auto &inter : _grid->getInterfaces())
  {
    const long &interfaceId = inter.getId();
    if (!inter.isBorder())
    {
      std::array<long, 2> owners = inter.getOwnerNeigh();

      int gr = _mapGlobal.at(owners[0]);
      int gl = _mapGlobal.at(owners[1]);

      std::vector<std::pair<long,double> > weights;

      // old version before Florian      weights = stencil.buildInterface(interfaceId);

      weights = _gradient->buildInterfaceNormalGradientStencil(interfaceId);

      double lambda_int = computeLambdaOnInterface(inter);

      for (std::size_t i = 0; i<weights.size(); ++i)
      {
        const long id = weights[i].first;
        int neigh = _mapGlobal.at(id);

        vp =  _grid->evalInterfaceArea(interfaceId)*weights[i].second/_grid->evalCellVolume(owners[0])*lambda_int;
        vm = -_grid->evalInterfaceArea(interfaceId)*weights[i].second/_grid->evalCellVolume(owners[1])*lambda_int;

        if (_grid->getCell(owners[0]).isInterior())
        {
          addMatrixValue(gr, neigh, vp);
        }
        if (_grid->getCell(owners[1]).isInterior())
        {
          addMatrixValue(gl, neigh, vm);
        }
      }
    } else {
      int gl = _mapGlobal.at(inter.getOwner());
      double vp = 0;
      if (_grid->getCell(inter.getOwner()).isInterior())
      {
        addMatrixValue(gl, gl, vp);
      }
    }
  }
}

void Laplacian::buildFVMatrix(PiercedVector<double>& kappaCC,
                              PiercedVector<double>& kappaFC,
                              double t,
                              Var type)
{
  Gradient grad(_grid->getDimension(), _grid, _stencils);
  BoundaryConditions* BC(_grid->getBoundaryConditions());

  auto& stencils = (_stencils->getInterfaceGradientStencil());
  Stencil stencil;

  for (auto &inter: _grid->getInterfaces())
  {
    const long &interId = inter.getId();
    const long& ownerId = inter.getOwner();
    int g0 = _mapGlobal.at(ownerId);

    if (_hasStencil)
    {
      stencil = stencils[interId];
    }
    else
    {
      stencil = grad.buildFaceNormalGradientStencil(interId);
    }

    const double& area = _grid->evalInterfaceArea(interId);

    double volume = _grid->evalCellVolume(ownerId);

    if (inter.isBorder())
    {
      if (_grid->getCell(ownerId).isInterior())
      {
        if ( BC->getBCType(inter.getOwnerFace(), type)
        == "Dirichlet"  )        //Ghost boundary cell
        {
          for (size_t i = 0; i<stencil.neighs.size(); i++)
          {
            const long cellId = stencil.neighs[i];

            const int isOutsideCell = stencil.isOutsideCell[i];
            double value;
            const double& weight = stencil.weights.weights[i];
            if( isOutsideCell < 0 )
            {
              int boundaryFace = std::abs(isOutsideCell + 1 );
              NPoint ghostOctCenter =  _grid->evalCellCentroid(cellId);
              ghostOctCenter[boundaryFace/2] += pow(-1, boundaryFace + 1)
              * _grid->evalCellSize(cellId);

              value = -kappaCC[cellId] * kappaFC[interId] * weight
              * BC->getBCFunction(boundaryFace, type)(
                ghostOctCenter,
                t
              ) * area/volume;
              addRHSValue(g0, value);
            }
            else
            {
              value = kappaCC[cellId] * kappaFC[interId] * weight
              * area/volume;
              addMatrixValue(g0, g0, value);
            }
          }
        }
      }
    }    //ELSE HOMOGENEOUS NEUMANN: Nothing to add
    else
    {
      const long& neighId = inter.getNeigh();

      //int g1 = _mapGlobal.at(neighId);
      for (size_t i = 0; i<stencil.neighs.size(); i++)
      {
        const long cellId = stencil.neighs[i];
        volume = _grid->evalCellVolume(ownerId);

        const int isOutsideCell = stencil.isOutsideCell[i];
        double value;
        const double& weight = stencil.weights.weights[i];

        if( isOutsideCell < 0 )
        {
          if(BC->getBCType(std::abs(isOutsideCell + 1 ), type)
             == "Dirichlet")            //Ghost boundary cell
          {
            int boundaryFace = std::abs(isOutsideCell + 1 );
            NPoint ghostOctCenter =  _grid->evalCellCentroid(cellId);
            ghostOctCenter[boundaryFace/2] += pow(-1, boundaryFace + 1)
                                              * _grid->evalCellSize(cellId);

            value = -kappaCC[cellId] * kappaFC[interId] * weight
                    * BC->getBCFunction(boundaryFace, type)(ghostOctCenter,
                                                            t) * area/volume;

            if (_grid->getCell(ownerId).isInterior())
            {
              addRHSValue(g0,value);
            }

            if (_grid->getCell(neighId).isInterior())
            {
              int g1 = _mapGlobal.at(neighId);
              value *= volume;

              volume = _grid->evalCellVolume(neighId);
              value  = -value/volume;

              addRHSValue(g1,value);
            }
          }
          //ELSE HOMOGENEOUS NEUMANN
        }
        else
        {
          int gn = _mapGlobal.at(cellId);
          value = kappaCC[cellId] * kappaFC[interId] * weight
                  * area/volume;

          if (_grid->getCell(ownerId).isInterior())
          {
            addMatrixValue(g0,gn,value);
          }

          if (_grid->getCell(neighId).isInterior())
          {
            int g1 = _mapGlobal.at(neighId);
            value *= volume;

            volume = _grid->evalCellVolume(neighId);
            value  = -value/volume;

            addMatrixValue(g1,gn,value);
          }
        }
      }
    }
  }
}

std::vector<std::array<double,3> > Laplacian::getGradientFromSolution()
{
  Gradient stencil(_grid->getDimension(), _grid);
  PiercedVector<double> sol = getSolutionPV();

  return stencil.computeLSGradient(sol);
}

void Laplacian::setDirichletCondition(double pos)
{
  double vp;

  _neumann = false;
  for (auto &inter : _grid->getInterfaces())
    {
    const long &id = inter.getId();
    if (inter.isBorder())
    {
      int cellId= inter.getOwner();
      int gl = _mapGlobal.at(cellId);

      double lambda_int = computeLambdaOnInterface(inter);

      vp = -_grid->evalInterfaceArea(id) / (_grid->evalCellSize(inter.getOwner()) * pos) / _grid->evalCellVolume(inter.getOwner())*lambda_int;
      if (_grid->getCell(inter.getOwner()).isInterior())
	    {
        addMatrixValue(gl,gl,vp);
      }
    }
  }
}

void Laplacian::setDirichletConditionRHS(double t)
{
  _neumann = false;
  for (auto &inter : _grid->getInterfaces())
  {
    if (inter.isBorder())
    {
      int gl = _mapGlobal.at(inter.getOwner());
      if (_grid->getCell(inter.getOwner()).isInterior())
      {
        double val, temp;
        temp = getRHSValue(gl);
        val = temp-2*t;
        setRHSValue(gl,val);
      }
    }
  }
}

void Laplacian::addSolutionToVTK(std::string tag) {

  std::vector<double> v = PVtoV(getSolutionPV(),_grid);
  _grid->getVTK().addData(tag, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, _sol_v);
}

double Laplacian::computeLambdaOnInterface(bitpit::Interface &inter) {
  double lambda_int = 1.0;

  if (_lambdaf) {

    NPoint xpos = _grid->evalInterfaceCentroid(inter.getId());
    lambda_int = (*_lambdaf)(xpos[0], xpos[1], xpos[2]);

  } else if (_lambda.size()){

     if (_gradient == NULL) {
       _gradient = new Gradient(_grid->getDimension(), _grid);
     }


    IInterpolator *interpo = InterpolatorFactory::get(interpoType::DISTANCEWEIGHTED);

    if (!inter.isBorder())
    {
      const long &interfaceId = inter.getId();

      std::array<long, 2> owners = inter.getOwnerNeigh();

      std::vector<std::pair<long,double> > weights;
      weights = _gradient->buildInterfaceNormalGradientStencil(interfaceId);

      NPoint xpos = _grid->evalInterfaceCentroid(interfaceId);

      lambda_int = (_lambda[owners[0]]+_lambda[owners[1]])/2.0;
      std::vector<NPoint> xref;
      std::vector<double> val_ref;

      for (std::size_t i = 0; i<weights.size(); ++i)
      {
        const long cellId = weights[i].first;

        xref.push_back(_grid->evalCellCentroid(cellId));
        val_ref.push_back(_lambda[cellId]);
      }

      if (weights.size()>0) {
        lambda_int = interpo->computeInterpolation(xpos, xref, val_ref);
      }
    } else {

      const long &interfaceId = inter.getId();
      const long &cellId = inter.getOwner();
      const bitpit::Cell cell = _grid->getCell(cellId);
      int norm = -1;
      double oppositeLambda = 0.0;
      int nbfaces = 1;

      NPoint interCentroid = _grid->evalInterfaceCentroid(interfaceId);
      NPoint cellCentroid = _grid->evalCellCentroid(cellId);

      NPoint oppositeCentroid;
      oppositeCentroid[NPX] = (cellCentroid[NPX]-interCentroid[NPX])+cellCentroid[NPX];
      oppositeCentroid[NPY] = (cellCentroid[NPY]-interCentroid[NPY])+cellCentroid[NPY];
      oppositeCentroid[NPZ] = (cellCentroid[NPZ]-interCentroid[NPZ])+cellCentroid[NPZ];

      if ((cellCentroid[NPY] == interCentroid[NPY]) && (cellCentroid[NPZ] == interCentroid[NPZ])) {
        norm = NPX;
      } else if ((cellCentroid[NPX] == interCentroid[NPX]) && (cellCentroid[NPZ] == interCentroid[NPZ])) {
        norm = NPY;
      } else {
        norm = NPZ;
      }

      for (auto &fid : _grid->findCellFaceNeighs(cellId)) {
        NPoint centroid = _grid->evalCellCentroid(fid);
        if ((
          (centroid[norm] < interCentroid[norm]) &&
          (centroid[norm] < oppositeCentroid[norm])
        ) || (
          (centroid[norm] > interCentroid[norm]) &&
          (centroid[norm] > oppositeCentroid[norm])
        ))
        oppositeLambda += _lambda[fid];
        nbfaces++;
      }

      oppositeLambda = oppositeLambda / nbfaces;

      if (nbfaces == 1) {
        lambda_int = _lambda[cellId] + (_lambda[cellId] - oppositeLambda)/2.0;
      } else {
        lambda_int =  _lambda[cellId] + (_lambda[cellId] - oppositeLambda)/1.5;
      }
    }
  }
  return lambda_int;
}

}
