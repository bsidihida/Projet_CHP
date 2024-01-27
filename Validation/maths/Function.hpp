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
 * @file Function.hpp
 * @brief List of math function
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-12-05
 * @copyright Inria
 */


#ifndef __FUNCTION_HPP__
#define __FUNCTION_HPP__

#include <array>
#include <vector>
#include <common.hpp>

double norm2(const NPoint &vec);
double norm2(const std::vector<double> &vec);
int inverse(int n, double* A);
int fullPiv(int n, double *A);
int partialPiv(int n, double *A);
int matVecProduct(int n,double *A, double *b);
int matVecProductNS(int lda, int cda,
                    double *A, double *X, double *Y,
                    bool transposeA=false,
                    double alpha=1., double beta=0.,
                    int incX=1, int incY=1);
int eiMatVecProduct(int n,double *A, double *b);
int matMatProduct(int lda, int cda, int ldb, int cdb,
                  double *A, double *B, double *C,
                  bool transposeA= false, bool transposeB= false,
                  double alpha = 1.0, double beta = 0.0);
int multTransposeLeft(int lda, int cda,
                      double *A, double *C);

#endif /* __FUNCTION_HPP__ */
