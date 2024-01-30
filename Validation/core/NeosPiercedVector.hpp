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

#ifndef __PIERCEDVECTOR_HPP__
#define __PIERCEDVECTOR_HPP__

/**
 * @file   PiercedVector.hpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Fri Aug 30 12:00:10 2019
 *
 * @brief  This file contains PiercedVector class
 *
 *
 */
#include <piercedVector.hpp>
#include <NeosAssert.hpp>

namespace neos {
/**
 * @class PiercedVector
 * @brief Class derived from bitpit::PiercedVector
 *
 * This class allows to overload some operators for PiercedVector.
 * People which uses this operators has to be careful whith manipulated
 * PiercedVector. Indeed, to add or substract two PiercedVector, they have
 * to have the same Ids or a "segmentation fault" will occured.
 */
template<typename value_t, typename id_t = long>
class PiercedVector : public bitpit::PiercedVector<value_t, id_t>
{
private:

public:
//PiercedVector();
using bitpit::PiercedVector<value_t, id_t>::PiercedVector;
PiercedVector(bitpit::PiercedVector<value_t, id_t>& v);

PiercedVector<value_t, id_t> operator+(const PiercedVector<value_t, id_t>&x);
PiercedVector<value_t, id_t> operator-(const PiercedVector<value_t, id_t>&x);
PiercedVector<value_t, id_t>& operator+=(const PiercedVector<value_t, id_t>&x);
PiercedVector<value_t, id_t>& operator-=(const PiercedVector<value_t, id_t>&x);

template<typename T>
PiercedVector<value_t, id_t> operator*(const T& alpha);
};

template<typename value_t, typename T, typename id_t = long>
PiercedVector<value_t, id_t> operator*(const T& alpha,
                                       const PiercedVector<value_t,id_t>& v);

}
#include "NeosPiercedVector.tpp"

#endif
