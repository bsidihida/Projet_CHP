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
 * @file   InputFileReader.cpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Fri Aug  2 11:54:58 2019
 *
 * @brief
 *
 *
 */

#include "InputFileReader.hpp"
//#include "bitpit.hpp"
#include <mpi.h>
#include "NeosAssert.hpp"

InputFileReader::InputFileReader(const string& filename) :
  m_filename(filename)
{
  this->read();
  this->showAllParameters();
}

void InputFileReader::read()
{
  //FIXME: Add a test to know if file is open
  //TODO: Define an ASSERT MACRO
  this->open();

  std::vector<std::string> lines;
  std::string line, key, value;

  while(std::getline(m_file, line))
  {
    lines.push_back(line);
  }
  for (auto i = 0u; i<lines.size(); i++)
  {
    line = lines.at(i);
    if (line.find(':') != line.npos && line.at(0) != '#')
    {
      std::stringstream ss(line);
      std::getline(ss, key, ':');
      std::getline(ss,value);

      // Remove spaces from keys
      key.erase(std::remove(key.begin(), key.end(), ' '), key.end());

      //arrange key and values
      trim(key);
      toUpper(key);
      trim(value);

      if (key.length()==0 and this->isMaster())
        std::cerr << std::endl << "\033[1m\033[34mWarning:\033[0m" << std::endl
                  << "InputFileReader: key is empty at line "
                  << i <<std::endl;
      if (value.length()==0 and this->isMaster())
        std::cerr << std::endl << "\033[1m\033[34mWarning:\033[0m" << std::endl
                  << "InputFileReader: value is empty at line "
                  << i << std::endl;


      m_keyVals.insert(std::make_pair(key, value));

    }
    else
    {
      continue;
    }
  }
  m_file.close();
}

double InputFileReader::getReal(const string& key)
{
  string KEY(key);


  //arrange key
  trim(KEY);
  toUpper(KEY);

  //Remove white spaces
  KEY.erase(std::remove(KEY.begin(),KEY.end(),' '),KEY.end());

  if(m_keyVals.find(KEY) == m_keyVals.end())
  {
    //FIXME: Definir des macros ici
    std::cerr<<"No key: "<<KEY<<" find in configuration file"<<std::endl;
    std::abort();
  }
  else
  {
    try
    {
      return (double)std::stod(m_keyVals.at(KEY));
    }
    catch(...)
    {
      //FIXME: Definir des macros ici
      std::cerr<<"InputFileReader::getReal cannot convert string"
               <<m_keyVals.at(KEY)<<" to double"<<std::endl;
      std::cerr<<"Value for key: "<<key
               <<" is set to default value"<<std::endl;
    }
  }
  ASSERT(false, "Unreachable statement");
  return -1.0;
}

bool InputFileReader::hasKey(const string& key)
{
  string KEY(key);


  //arrange key
  trim(KEY);
  toUpper(KEY);

  //Remove white spaces
  KEY.erase(std::remove(KEY.begin(),KEY.end(),' '),KEY.end());

  return m_keyVals.find(KEY) != m_keyVals.end();
}

string InputFileReader::getString(const string& key)
{
  string KEY(key);

  //arrange key
  trim(KEY);
  toUpper(KEY);

  //remove white spaces
  KEY.erase(std::remove(KEY.begin(),KEY.end(),' '),KEY.end());
  if (m_keyVals.find(KEY) == m_keyVals.end())
  {
    //FIXME: Definir des macros ici
    std::cerr<<"No key: "<<KEY<<" find in configuration file"<<std::endl;
    std::abort();
  }

  return m_keyVals.at(KEY);
}

bool InputFileReader::open()
{
  m_file.open(m_filename.c_str(), std::ifstream::in);
  if (!m_file.good())
  {
    std::cerr << "Could not open the file: "<<m_filename<<std::endl;
    return false;
  }
  return true;
}

void InputFileReader::toUpper(string& s)
{
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);
}

void InputFileReader::trim(string& s)
{
  // trim string by removing tab
  size_t posBegin = s.find_first_not_of(" \t");
  size_t posEnd   = s.find_last_not_of (" \t")+1;

  s = (posEnd > posBegin ? s.substr(posBegin, posEnd-posBegin) : string(""));

}

void InputFileReader::showAllParameters()
{
  for (auto& pair: m_keyVals)
  {
    std::cout<< "key: "<<pair.first
             << " value: "<<pair.second
             <<std::endl;
  }
}

bool InputFileReader::isMaster()
{
  #if defined (ENABLE_MPI)
  int rank(0.);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    return true;
  }
  return false;
  #else
  return true;
  #endif
}

string InputFileReader::operator()(const string& key)
{
  return this->getString(key);
}
