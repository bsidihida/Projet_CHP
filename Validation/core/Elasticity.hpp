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

#ifndef __ELASTICITY_HPP__
#define __ELASTICITY_HPP__

/**
 * @file Elasticity
 * @brief This file contains the Elasticity class
 * @author Antoine Gerard
 * @date 2019-08-01
 * @copyright Inria
 */

/**
 *@class Elasticity
 *@brief Class to handle elastics solids
 *
 *
 **/

#include "Grid.hpp"
#include "Transport.hpp"
#include "NeosPiercedVector.hpp"

namespace neos {


struct ElasticStructures
{
  std::vector<PiercedVector<double> > InvDef;
  PiercedVector<double> LevSet;
};

class Elasticity
{
private:
Grid *m_grid;                   /**< Pointer to the current domain of computation*/
public:
/**
 * @brief Default constructor
 *
 * @param grid: pointer to current domain of computation
 */
Elasticity(Grid* grid);

/**
 * Destructor
 *
 */
~Elasticity(){
  ;
}


ElasticStructures transportElasticStuctures(Transport& transport,
                                            ElasticStructures& elastic,
                                            ElasticStructures& elasticPrev,
                                            PiercedVector<double>& velFC,
                                            double t,
                                            std::vector<double> timeSteps,
                                            Var type = Var::LS);

};
}

#endif
