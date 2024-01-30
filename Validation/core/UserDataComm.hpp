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
 * @file   UserDataComm.hpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Thu Oct 24 16:23:41 2019
 *
 * @brief  This file contains UserDataComm class
 *
 *
 */


#ifndef __USERDATACOMM_HPP
#define __USERDATACOMM_HPP

#if defined (BITPIT_ENABLE_MPI)

// was #if ENABLE_BITPIT==1

//#include <bitpit.hpp>
#include <communications.hpp>
#include <NeosPiercedVector.hpp>
#include <Grid.hpp>

using bitpit::DataCommunicator;

namespace neos {

/* @brief This class handles mpi communications*/
template<typename T>
class UserDataComm : public DataCommunicator
{
public:

/**
 * @brief Constructor
 *
 * @param grid: Reference to the current grid
 * @param dataSize: number of value which will be send
 *
 * @return
 */
UserDataComm(const Grid &grid, int dataSize);
UserDataComm(const Grid &grid);

/**
 * @brief Destructor
 */
~UserDataComm();

/**
 * @brief Get the data size which will be communicated in bytes.
 *
 *
 * @return Data size in bytes
 */
size_t getDataSize();

/**
 * @brief Prepare data which will be send to other processes
 *
 * @param[in] cellData: Data to send
 */
void prepareSend(PiercedVector<T> &cellData);

/**
 * @brief: Receive data send by other processes
 *
 * @param[in,out] cellData: Where to put received data
 */
void receiveDataset(PiercedVector<T> &cellData);

/**
 * @brief Communicate data => Send and receive are done here
 *
 * @param[in,out] cellData: PiercedVector to send and also receiver.
 */
void communicate(PiercedVector<T> &cellData);

/**
 * @brief Prepare data which will be send to other processes
 *
 * @param[in] cellData: Data to send
 */
void prepareSend(PiercedVector<std::vector<T> > &cellData);

/**
 * @brief: Receive data send by other processes
 *
 * @param[in,out] cellData: Where to put received data
 */
void receiveDataset(PiercedVector<std::vector<T> > &cellData);


void communicate(PiercedVector<std::vector<T> > &cellData);

/**
 * @brief Prepare data which will be send to other processes
 *
 * @param[in] cellData: Data to send
 */
void prepareSend(PiercedVector<std::array<T,3> > &cellData);

/**
 * @brief Prepare data which will be send to other processes
 *
 * @param[in] cellData: Data to send
 */
void receiveDataset(PiercedVector<std::array<T,3> > &cellData);

/**
 * @brief Communicate data => Send and receive are done here
 *
 * @param[in,out] cellData: PiercedVector to send and also receiver.
 */
void communicate(PiercedVector<std::array<T,3> > &cellData);

/**
 * @brief Update information needed for ghost data exchange
 *
 */
void update();
void update(int dataSize);

private:
const Grid &m_grid;             /**< Reference to the grid */
int m_dataSize;                 /**< Size of data in bytes */
int m_rank;                     /**< Rank of process */
};
}

#include "UserDataComm.tpp"

#endif

#endif
