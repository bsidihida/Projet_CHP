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

#ifndef __STENCIL_HPP__
#define __STENCIL_HPP__

#include <Grid.hpp>
#include <common.hpp>
#include <NeosPiercedVector.hpp>

namespace neos {

struct Stencil_weights
{
  std::vector<double> weights;
  std::vector<double> weights_dx;
  std::vector<double> weights_dy;
  std::vector<double> weights_CCToFC;
};

struct Stencil
{
  Stencil_weights weights;
  std::vector<long> neighs;
  std::vector<int> isOutsideCell;
};

class StencilBuilder
{
public:
StencilBuilder(Grid* grid);
~StencilBuilder(){
  ;
}

void buildCellGradientStencil();
const PiercedVector<Stencil>& getCellGradientStencil() const;

void buildInterfaceGradientStencil();
const PiercedVector<Stencil>& getInterfaceGradientStencil() const;

void buildCCToFCStencils(interpoType interpType = interpoType::RBF);
const Stencil buildCCToFCStencil(const long& interId,
                            interpoType interpType = interpoType::RBF);
const PiercedVector<Stencil>& getCCToFCStencil() const;

private:
Grid* _grid;
PiercedVector<Stencil> _stencilCC;
PiercedVector<Stencil> _stencilFC;
PiercedVector<Stencil> _stencilCCToFC;
};

}

#endif
