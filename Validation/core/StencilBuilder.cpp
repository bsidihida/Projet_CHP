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

#include "StencilBuilder.hpp"
#include "Gradient.hpp"
#include "RBF.hpp"
//#include "Projection.hpp"

namespace neos {

StencilBuilder::StencilBuilder(Grid* grid) : _grid(grid),
  _stencilCC({}),
  _stencilFC({}),
  _stencilCCToFC({})
{
}

/* ----------------- Cell Centered gradient stencil ----------------- */

void StencilBuilder::buildCellGradientStencil()
{
  Gradient grad(_grid->getDimension(), _grid);
  for (auto &cell: _grid->getCells())
  {
    const long &cellId = cell.getId();
    _stencilCC.emplace(cellId);

    if (cell.isInterior())
      _stencilCC[cellId] = grad.buildCellGradientStencil(cellId);
  }
}

const PiercedVector<Stencil>& StencilBuilder::getCellGradientStencil() const
{
  return _stencilCC;

}

/* ----------------- face Centered gradient stencil ----------------- */

void StencilBuilder::buildInterfaceGradientStencil()
{
  Gradient grad(_grid->getDimension(), _grid);
  for (auto &inter: _grid->getInterfaces())
  {
    const long &interId = inter.getId();
    _stencilFC.emplace(interId);
    _stencilFC[interId] = grad.buildFaceNormalGradientStencil(interId);
  }
}

const PiercedVector<Stencil>& StencilBuilder::getInterfaceGradientStencil() const
{
  return _stencilFC;

}

/* --------- Cell Centered to face centered interpolation stencil ------- */
void StencilBuilder::buildCCToFCStencils(interpoType interpType)
{
//  Projection proj(_grid);
  Gradient grad(_grid->getDimension(), _grid);
  for (auto &inter: _grid->getInterfaces())
  {
    const long &interId = inter.getId();
    if (!inter.isBorder())
    {
      _stencilCCToFC.emplace(interId);
      _stencilCCToFC[interId] = buildCCToFCStencil(interId, interpType);
    }
  }
}

const Stencil StencilBuilder::buildCCToFCStencil(const long& interId,
                                       interpoType interpType)
{

  // Careful: stencil only built for non border interfaces

  Stencil stencil;
  std::vector<NPoint> posNeighs;

  std::vector<int> vertices;
  std::array<long, 2> owners = _grid->getInterface(interId).getOwnerNeigh();

  int face = _grid->getInterface(interId).getOwnerFace();
  if ( _grid->getCellLevel(owners[0]) !=
       _grid->getCellLevel(owners[1]) )
  {
    stencil.neighs.push_back(owners[0]);
    posNeighs.push_back(_grid->evalCellCentroid(owners[0]));


    stencil.neighs.push_back(owners[1]);
    posNeighs.push_back(_grid->evalCellCentroid(owners[1]));

    // FIXME: Owners[0] is always the finest octant ?
    vertices = _grid->getFaceVertexLocalIds(face);
    for (auto &v: vertices)
    {
      auto neighBuff = _grid->findCellVertexNeighs(owners[0], v);
      // FIXME: Pass vector neighs as pointer for bitpit function
      // instead of testing with std::find
      for (auto &n: neighBuff)
      {
        if (std::find(stencil.neighs.begin(), stencil.neighs.end(), n)
            == stencil.neighs.end())
        {
          stencil.neighs.push_back(n);
          posNeighs.push_back(_grid->evalCellCentroid(n));
        }
      }
    }

    IInterpolator* interpo = InterpolatorFactory::get(interpType);
    if (interpType == interpoType::RBF)
      ((Rbf*)interpo)->setEpsilon(0.01/_grid->evalInterfaceArea(interId));

    stencil.weights.weights_CCToFC =
      interpo->computeInterpolation(_grid->evalInterfaceCentroid(interId),
                                    posNeighs);
  }
  else
  {
    stencil.neighs.push_back(owners[0]);
    stencil.weights.weights_CCToFC.push_back(0.5);

    stencil.neighs.push_back(owners[1]);
    stencil.weights.weights_CCToFC.push_back(0.5);
  }

  return stencil;
}

const PiercedVector<Stencil>& StencilBuilder::getCCToFCStencil() const
{
  return _stencilCCToFC;
}
}
