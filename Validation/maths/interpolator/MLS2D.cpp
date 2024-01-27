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
 * @file   MLS2D.cpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Tue Aug  6 11:26:50 2019
 *
 * @brief  This file contains Compact Interpolator Class
 *
 *
 */

#include <iostream>
#include "MLS2D.hpp"

namespace neos {

MLS2D::MLS2D() {
#if USE_EIGEN3
  _solve = fullPiv;
  _matVecProduct = eiMatVecProduct;
#else
  _matVecProduct = matVecProduct;
#endif /* USE_EIGEN3 */
#ifdef USE_LAPACKE
  _solve = inverse;
#endif /* USE_LAPACKE */
}

std::vector<double> MLS2D::computeInterpolation(const NPoint &xpos,
                                                const std::vector<NPoint> &xref)
{
  std::vector<double> weights;
  weights = this->computeInterpolation2ndOrder(xpos, xref);
  return weights;
}

double MLS2D::computeInterpolation(const NPoint &xpos,
                                   const std::vector<NPoint > &xref,
                                   const std::vector<double> &vref)
{
  //Get weights
  std::vector<double> weightsMat = this->computeInterpolation(xpos, xref);
  std::vector<double> res(6);
  std::vector<double> b(vref);

  //Multiply by vref to obtain weights
  matVecProductNS(6, xref.size(),
                  weightsMat.data(), b.data(), res.data());

  return res[0];
}

std::vector<double> MLS2D::computeFirstDerivatives(const NPoint &xpos,
                                                   const std::vector<NPoint > &xref,
                                                   const std::vector<double> &vref)
{
  //Get weights
  std::vector<double> weightsMat = this->computeInterpolation(xpos, xref);
  std::vector<double> res(6);
  std::vector<double> b(vref);

  //Multiply by vref to obtain weights
  matVecProductNS(6, xref.size(),
                  weightsMat.data(), b.data(), res.data());

  return {res[1], res[2]};
}

std::vector<double> MLS2D::computeSecondDerivatives(const NPoint &xpos,
                                                    const std::vector<NPoint > &xref,
                                                    const std::vector<double> &vref)
{
  //Get weights
  std::vector<double> weightsMat = this->computeInterpolation(xpos, xref);
  std::vector<double> res(6);
  std::vector<double> b(vref);

  //Multiply by vref to obtain weights
  matVecProductNS(6, xref.size(),
                  weightsMat.data(), b.data(), res.data());

  return {res[3], res[4], res[5]};
}

std::vector<double> MLS2D::computeInterpolation2ndOrder(const NPoint &xpos,
                                                        const std::vector<NPoint> &xref)
{
  int Nneighs = xref.size();
  auto A = new double[6*Nneighs];
  double *tAA = new double[6*6]();
  std::vector<double> weights(6*Nneighs);

  //FIXME: Put an assert macro here
  if (Nneighs < 6)
  {
    // You must have more than 6 points to do MLS interpolation
    std::cerr<<"Problem in MLS interpolation : interpolation with less than "
             << Nneighs << " points...( < 6 ) \n";
    std::abort();
  }
  // Construction of the Nn*6 matrix
  for (int i=0; i < Nneighs; i++)
  {
    A[i + 0*Nneighs] = 1;
    A[i + 1*Nneighs] = xref[i][0] - xpos[0];
    A[i + 2*Nneighs] = xref[i][1] - xpos[1];
    A[i + 3*Nneighs] = (xref[i][0] - xpos[0]) * (xref[i][1] - xpos[1]);
    A[i + 4*Nneighs] = 0.5 * (xref[i][0] - xpos[0])
                       * (xref[i][0] - xpos[0]);
    A[i + 5*Nneighs] = 0.5 * (xref[i][1] - xpos[1])
                       * (xref[i][1] - xpos[1]);;
  }

  // tAA = A^{T} * A
  multTransposeLeft(Nneighs, 6, A, tAA);


  // tAA = (A^{T} * A)^{-1}
  _solve(6, tAA);

  // weights = (A^{T} * A)^{-1} * A^{T}
  matMatProduct(6, 6,
                Nneighs, 6,
                tAA, A,
                weights.data(),
                false, true);

  //Delete pointers
  delete[] A;
  delete[] tAA;

  return weights;
}

std::vector<double> MLS2D::computeInterpolationAtFaceCenter(const NPoint &xpos,
                                                            const std::vector<NPoint> &xref)
{
  int Nneighs = xref.size();
  double *A = new double[3*Nneighs];
  double *tAA = new double[3*3]();
  std::vector<double> weights(3*Nneighs);

  //FIXME: Put an assert macro here
  if (Nneighs < 3)
  {
    // You must have more than 4 points to do MLS interpolation at Face
    // center
    std::cerr<<"Problem in MLS interpolation : interpolation with less than "
             << Nneighs << " points...( < 4 ) \n";
    std::abort();
  }

  // Construction of the Nn*4 matrix
  for (int i=0; i < Nneighs; i++)
  {
    A[i + 0*Nneighs] = 1;
    A[i + 1*Nneighs] = xref[i][0] - xpos[0];
    A[i + 2*Nneighs] = xref[i][1] - xpos[1];
  }

  // tAA = A^{T} * A
  multTransposeLeft(Nneighs, 3, A, tAA);


  // tAA = (A^{T} * A)^{-1}
  _solve(3, tAA);

  // weights = (A^{T} * A)^{-1} * A^{T}
  matMatProduct(3, 3,
                Nneighs, 3,
                tAA, A,
                weights.data(),
                false, true);

  //Delete pointers
  delete[] A;
  delete[] tAA;

  return weights;
}
}
