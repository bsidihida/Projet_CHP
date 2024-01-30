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

namespace neos {

  template <typename T>
  UserDataComm<T>::UserDataComm(const Grid &grid, int dataSize)
    : DataCommunicator(grid.getCommunicator()),
    m_grid(grid),
    m_dataSize(dataSize*sizeof(T))
  {
    // Get MPI information
    MPI_Comm_rank(grid.getCommunicator(), &m_rank);
    update(dataSize*sizeof(T));
  }

  template <typename T>
  UserDataComm<T>::UserDataComm(const Grid &grid)
    : DataCommunicator(grid.getCommunicator()),
    m_grid(grid)
  {
    update((grid.getDimension()+2)*sizeof(T));
  }

  template <typename T>
  UserDataComm<T>::~UserDataComm()
  {
    std::cout<<m_rank<<" call destructor"<<std::endl;
  }

/*!
        Gets the size in bytes of the data that will be communicated.

        \result The size in bytes of the data that will be communicated
 */
  template <typename T>
  size_t UserDataComm<T>::getDataSize()
  {
    return m_dataSize;
  }

/*!
        Updates the information needed for ghost data exchange

        \param dataSize is the size in bytes of the data that will be communicated
 */
  template <typename T>
  void UserDataComm<T>::update(int dataSize)
  {
    m_dataSize = dataSize * sizeof(T);
    update();
  }

/*!
        Updates the information needed for ghost data exchange
 */
  template <typename T>
  void UserDataComm<T>::update()
  {
    for (const auto entry : m_grid.getGhostCellExchangeSources())
    {
      int rank = entry.first;
      const std::vector<long> &sources = entry.second;
      long bufferSize = sources.size() * m_dataSize;
      this->setSend(rank, bufferSize);
    }
  }


/*!
        Send ghosts data using non-blocking communications
 */
  template <typename T>
  void UserDataComm<T>::prepareSend(PiercedVector<T> &cellData)
  {
    // Fill the buffer with the given field and start sending the data
    for (const auto entry : m_grid.getGhostCellExchangeSources()) {
      // Get the send buffer
      short rank = entry.first;
      bitpit::SendBuffer buffer = this->getSendBuffer(rank);

      // Store data in the buffer
      const std::vector<long> &sources = entry.second;
      for (const long id : sources)
      {
        buffer<<cellData[id];
      }
    }
  }

/*!
   Receive ghosts data using non-blocking communications
 */
  template <typename T>
  void UserDataComm<T>::receiveDataset(PiercedVector<T> &cellData)
  {
    int nCompletedRecvs = 0;
    std::vector<int> recvRanks = this->getRecvRanks();
    std::sort(recvRanks.begin(),recvRanks.end());
    for(auto rank: recvRanks)
    {
      this->waitRecv(rank);
      bitpit::RecvBuffer buffer = this->getRecvBuffer(rank);
      for (long id : m_grid.getGhostCellExchangeTargets(rank))
      {
        if (!cellData.exists(id))
          cellData.emplace(id);
        buffer>>cellData[id];
      }

      ++nCompletedRecvs;
    }
  }

/*!
        Send and receive ghosts in the same PiercedVector
 */
  template <typename T>
  void UserDataComm<T>::communicate(PiercedVector<T> &cellData)
  {
    // Prepare buffers to send ans start sends
    prepareSend(cellData);

    // Let's each process discover received data
    // and start receives.
    this->discoverRecvs();
    this->startAllRecvs();

    this->startAllSends();
    // Write received data in cellData
    receiveDataset(cellData);

    this->waitAllSends();

  }

/*!
        Send ghosts data using non-blocking communications
 */
  template <typename T>
  void UserDataComm<T>::prepareSend(PiercedVector<std::vector<T>> &cellData)
  {
    // Wait previous sends
    //waitAllSends();

    // Fill the buffer with the given field and start sending the data
    for (const auto entry : m_grid.getGhostCellExchangeSources())
    {
      // Get the send buffer
      short rank = entry.first;
      bitpit::SendBuffer buffer = this->getSendBuffer(rank);

      // Store data in the buffer
      const std::vector<long> &sources = entry.second;
      for (const long id : sources)
      {
        for (uint64_t i=0; i<cellData[id].size(); ++i)
          buffer<<cellData[id][i];
      }
    }
  }

/*!
        Receive ghosts data using non-blocking communications
 */
  template <typename T>
  void UserDataComm<T>::receiveDataset(PiercedVector<std::vector<T>> &cellData)
  {
    int nCompletedRecvs = 0;
    std::vector<int> recvRanks = this->getRecvRanks();
    std::sort(recvRanks.begin(),recvRanks.end());
    for(auto rank: recvRanks) {
      this->waitRecv(rank);
      bitpit::RecvBuffer buffer = this->getRecvBuffer(rank);
      for (long id : m_grid.getGhostExchangeTargets(rank))
      {
        if (!cellData.exists(id))
          cellData.emplace(id);
        for (uint64_t i=0; i<cellData[id].size(); ++i)
          buffer>>cellData[id][i];
      }

      ++nCompletedRecvs;
    }
  }

/*!
        Send and receive ghosts in the same PiercedVector
 */
  template <typename T>
  void UserDataComm<T>::communicate(PiercedVector<std::vector<T>> &cellData)
  {
    // Prepare buffers to send ans start sends
    prepareSend(cellData);

    // Let's each process discover received data
    // and start receives.
    this->discoverRecvs();
    this->startAllRecvs();

    this->startAllSends();

    // Write received data in cellData
    receiveDataset(cellData);

    this->waitAllSends();
  }

/*!
        Send ghosts data using non-blocking communications
 */
  template <typename T>
  void UserDataComm<T>::prepareSend(PiercedVector<std::array<T,3>> &cellData)
  {
    // Fill the buffer with the given field and start sending the data
    for (const auto entry : m_grid.getGhostCellExchangeSources())
    {
      // Get the send buffer
      short rank = entry.first;
      bitpit::SendBuffer buffer = this->getSendBuffer(rank);
      // Store data in the buffer
      const std::vector<long> &sources = entry.second;
      for (const long id : sources)
      {
        for (int i=0; i<3; ++i)
        {
          buffer<<cellData[id][i];
        }
      }
    }
  }

/*!
        Receive ghosts data using non-blocking communications
 */
  template <typename T>
  void UserDataComm<T>::receiveDataset(PiercedVector<std::array<T,3>> &cellData)
  {
    int nCompletedRecvs = 0;
    std::vector<int> recvRanks = this->getRecvRanks();
    std::sort(recvRanks.begin(),recvRanks.end());
    for(auto rank: recvRanks)
    {
      this->waitRecv(rank);
      bitpit::RecvBuffer buffer = this->getRecvBuffer(rank);
      for (long id : m_grid.getGhostCellExchangeTargets(rank))
      {
        if (!cellData.exists(id))
        {
          cellData.emplace(id);
        }

        for (uint64_t i=0; i<3; ++i)
        {
          buffer>>cellData[id][i];
        }
      }

      ++nCompletedRecvs;
    }
  }

/*!
        Send and receive ghosts in the same PiercedVector
 */
  template <typename T>
  void UserDataComm<T>::communicate(PiercedVector<std::array<T,3>> &cellData)
  {
    // Prepare buffers to send ans start sends
    prepareSend(cellData);

    // Let's each process discover received data
    // and start receives.
    this->discoverRecvs();
    this->startAllRecvs();

    this->startAllSends();

    // Write received data in cellData
    receiveDataset(cellData);

    this->waitAllSends();
  }
}
