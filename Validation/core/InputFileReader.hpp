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
 * @file   InputFileReader.hpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Fri Aug  2 11:45:53 2019
 *
 * @brief  This file contains InputFileReader class
 *
 *
 */

#ifndef __INPUTFILEREADER_HPP__
#define __INPUTFILEREADER_HPP__

#include <cerrno>
#include <utility>
#include <map>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <stdlib.h>

using std::string;

class InputFileReader
{
public:

/** @brief Constructor with file name
 *  @param inputFile   Can be absolute or relative file name
 */
InputFileReader(const string& inputFile);


/**
 * @brief Destructor
 *
 */
~InputFileReader(){
  ;
}

/**
 * @brief function to read an input file
 *
 */
void read();

/**
 * @brief Transform a key string into a double
 *
 * @param key: string to transform into a double
 *
 * @return: key as a double
 */
double getReal (const string& key);

/**
 * @brief Return string key in a formatted way
 *
 * @param key: string to arrange
 *
 * @return formatted key
 */
string  getString (const string& key);

/**
 * @brief Check if key is available
 *
 * @param key: string to analyze (not case sensitive)
 *
 * @return true if key exists
 */
bool hasKey(const string& key);

/**
 * Check if the file is well opened
 *
 *
 * @return True if file is well opened
 */
bool open();

/**
 * Trim left and right spaces
 *
 * @param s: string to trim
 */
void trim(string& s);

/**
 * Put all letters of s in capital.
 *
 * @param s: string to capitalize
 */
void toUpper(string& s);

/**
 * Function to display all parameters
 *
 */
void showAllParameters();

bool isMaster();

/**
 * Overload operator to get formatted key
 *
 *
 * @return formatted key
 */
string operator()(const string& key);



private:
std::map<string, string> m_keyVals;     /**< The parameters stored as strings*/
bool m_warnDefault;     /**< Print missing values warnings if true */
std::string m_filename;       /**< file to open */
std::ifstream m_file;       /**< stream */
};

#endif
