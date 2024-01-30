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
 * @file   SimulationParameters.hpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Fri Aug  2 11:31:01 2019
 *
 * @brief  This file contains the SimulationParameters class
 *
 *
 */

#ifndef __SIMULATIONPARAMETERS_HPP__
#define __SIMULATIONPARAMETERS_HPP__

#include <iostream>
#include <vector>

#include "InputFileReader.hpp"

using std::string;

/**
 * @class SimulationParameters
 * @brief Class to handle useful parameters for simulations
 */
class SimulationParameters
{
public:
/**
 * @brief Default constructor
 *
 */
SimulationParameters();

/**
 * @brief Destructor
 *
 */
~SimulationParameters(){
  ;
}

/**
 * @brief Read parameters from an input file
 *
 * @param inputFile: name of file to read
 */
void readFromInputFile(const string& inputFile);

void setDefaultParameters();

/**
 * @brief Set origin of the domain from a string
 *
 * @param s: string containing X Y Z coordinates of the origin
 */
void setDomainOriginFromString(const string& s);

/**
 *
 * @brief Get origin of the domain
 *
 * @return origin of the domain as a vector (X,Y,Z)
 */
const std::vector<double>& getDomainOrigin() const;

/**
 * @brief Get Shear Modulus
 *
 *
 * @return Shear Modulus as a double
 */
const double& getShearModulus() const;

/**
 * @brief Get CFL condition
 *
 *
 * @return cfl condition
 */
const double& getCFL() const;

/**
 * @brief get density of the solid
 *
 *
 * @return density of the solid
 */
const double& getSolidDensity() const;

/**
 * @brief get final time of simulation
 *
 *
 * @return final time of simulation
 */
const double& getTmax() const;

/**
 * @brief get density of the fluid
 *
 *
 * @return density of the fluid
 */
const double& getFluidDensity() const;

/**
 * @brief get kinematic viscosity of the fluid
 *
 *
 * @return kinematic viscosity of the fluid
 */
double getFluidKinematicViscosity() const;

protected:
std::vector<double> m_domainOrigin;     /**< Origin of the domain */
double m_length;                /**< Length of domain */
double m_shearModulus;          /**< Shear Modulus for elastic solids */
double m_CFL;                   /**< CFL condition */
double m_rhoS;                  /**< Solid density */
double m_Tmax;                  /**< Final time */
double m_rhoF;                  /**< Fluid density */
double m_muF;                   /**< Fluid viscosity */
};

#endif
