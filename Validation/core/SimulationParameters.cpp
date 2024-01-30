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
 * @file   SimulationParameters.cpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Fri Aug  2 12:40:24 2019
 *
 * @brief
 *
 *
 */

#include "SimulationParameters.hpp"

SimulationParameters::SimulationParameters() : m_domainOrigin(3, 0.),
  m_length(1.),
  m_CFL(0.5),
  m_rhoS(1.),
  m_Tmax(0.),
  m_rhoF(1.),
  m_muF(1.)
{
  //FIXME: Default value ?
}

void SimulationParameters::readFromInputFile(const string& inputFile)
{
  InputFileReader reader(inputFile);

  if (reader.hasKey("Shear Modulus"))
    m_shearModulus = reader.getReal("Shear Modulus");

  if (reader.hasKey("Domain Origin"))
  {
    this->setDomainOriginFromString(reader("Domain Origin"));
  }

  if (reader.hasKey("Solid Density"))
  {
    m_rhoS = reader.getReal("Solid Density");

  }
  if (reader.hasKey("Final time"))
  {
    m_Tmax = reader.getReal("Final time");
  }
  if (reader.hasKey("CFL"))
  {
    m_CFL = reader.getReal("CFL");
  }
  if (reader.hasKey("Fluid Density"))
  {
    m_rhoF = reader.getReal("Fluid Density");
  }
  if (reader.hasKey("Fluid Viscosity"))
  {
    m_muF = reader.getReal("Fluid Viscosity");
  }
}

void SimulationParameters::setDomainOriginFromString(const string& s)
{
  string s2(s);
  std::istringstream iss(s2);

  //split string into a vector
  std::vector<std::string> originStr(std::istream_iterator<std::string>(iss),{}
                                     );

  std::vector<double> origin(originStr.size());
  std::transform(originStr.begin(),
                 originStr.end(),
                 origin.begin(),
                 [](const std::string& val){
    return std::stod(val);
  });

  m_domainOrigin = origin;
}

const std::vector<double>& SimulationParameters::getDomainOrigin() const
{
  return m_domainOrigin;
}



const double& SimulationParameters::getShearModulus() const
{
  return m_shearModulus;
}

const double& SimulationParameters::getSolidDensity() const
{
  return m_rhoS;
}

const double& SimulationParameters::getCFL() const
{
  return m_CFL;
}

const double& SimulationParameters::getTmax() const
{
  return m_Tmax;
}

const double& SimulationParameters::getFluidDensity() const
{
  return m_rhoF;
}

double SimulationParameters::getFluidKinematicViscosity() const
{
  return m_muF/m_rhoF;
}
