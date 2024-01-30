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

//template<typename value_t, typename id_t>
//PiercedVector<value_t, id_t>::PiercedVector():bitpit::PiercedVector<value_t, id_t>()
//{}

namespace neos {

  template<typename value_t, typename id_t>
  PiercedVector<value_t, id_t>::PiercedVector(bitpit::PiercedVector<value_t, id_t>&v) : bitpit::PiercedVector<value_t, id_t>()
  {
    const std::vector<id_t> ids = v.getIds();
    for (auto &id: ids)
    //for (auto &id: v.getIds)
    {
      this->emplace(id);
      this->at(id) = v[id];
    }
  }

  template<typename value_t, typename id_t>
  PiercedVector<value_t, id_t> PiercedVector<value_t, id_t>::operator+(const PiercedVector<value_t, id_t>&x)
  {
    PiercedVector<value_t, id_t> res;
    bool hasId(true);
    ASSERT(this->size() == x.size(),
           "You can't add PiercedVector which have not the same size");
    for (auto& id: this->getIds())
    {
      if (x.find(id) == x.end())
        hasId = false;

      ASSERT(hasId,
             "You can't add PiercedVector which have not the same Ids");
      res.emplace(id, this->at(id) + x.at(id));
    }
    return res;
  }

  template<typename value_t, typename id_t>
  PiercedVector<value_t, id_t> PiercedVector<value_t, id_t>::operator-(const PiercedVector<value_t, id_t>&x)
  {
    PiercedVector<value_t, id_t> res;
    bool hasId(true);
    ASSERT(this->size() == x.size(),
           "You can't substract PiercedVector which have not the same size");
    for (auto& id: this->getIds())
    {
      if (x.find(id) == x.end())
        hasId = false;
      ASSERT(hasId,
             "You can't substract PiercedVector which have not the same Ids");
      res.emplace(id, this->at(id) - x.at(id));
    }
    return res;
  }

  template<typename value_t, typename id_t>
  PiercedVector<value_t, id_t>&PiercedVector<value_t, id_t>::operator+=(const PiercedVector<value_t, id_t>&x)
  {
    bool hasId(true);
    ASSERT(this->size() == x.size(),
           "You can't substract PiercedVector which have not the same size");
    for (auto& id: this->getIds())
    {
      if (x.find(id) == x.end())
        hasId = false;
      ASSERT(hasId,
             "You can't substract PiercedVector which have not the same Ids");
      this->at(id) += x.at(id);
    }
    return *this;
  }

  template<typename value_t, typename id_t>
  PiercedVector<value_t, id_t>&PiercedVector<value_t, id_t>::operator-=(const PiercedVector<value_t, id_t>&x)
  {
    bool hasId(true);
    ASSERT(this->size() == x.size(),
           "You can't substract PiercedVector which have not the same size");
    for (auto& id: this->getIds())
    {
      if (x.find(id) == x.end())
        hasId = false;
      ASSERT(hasId,
             "You can't substract PiercedVector which have not the same Ids");
      this->at(id) -= x.at(id);
    }
    return *this;
  }

  template<typename value_t, typename id_t>
  template<typename T>
  PiercedVector<value_t, id_t> PiercedVector<value_t,id_t>::operator*(const T& alpha)
  {
    PiercedVector<value_t, id_t> res;
    for (auto& id: this->getIds())
    {
      res.emplace(id, alpha * this->at(id));
    }
    return res;
  }

  template<typename value_t, typename T, typename id_t>
  PiercedVector<value_t, id_t> operator*(const T&
                                         alpha,  const PiercedVector<value_t,id_t>&v)
  {
    PiercedVector<value_t, id_t> res;
    for (auto& id: v.getIds())
    {
      res.emplace(id, alpha * v.at(id));
    }
    return res;
  }
}
