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
 * @file Function.cpp
 * @brief List of math function
 * @author Matias Hastaran
 * @version 0.1
 * @date 2017-12-05
 * @copyright Inria
 */

#include <math.h>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string.h>
#include "Function.hpp"

#ifdef USE_LAPACKE
#include "lapacke.h"
extern "C"
{
void dgemm_(char *transa, char *transb,
            int *m, int *n, int *k,
            double *alpha, double *a, int *lda,
            double *b, int *ldb, double *beta,
            double *c, int *ldc );

void dgemv_(char* TRANS, const int* M, const int* N,
            double* alpha, double* A, const int* LDA, double* X,
            const int* INCX, double* beta, double* C, const int* INCY);
}
#include "cblas.h"
#endif /* USE_LAPACKE */

#ifdef USE_EIGEN3
#include "Eigen/Dense"
#include "Eigen/LU"
#endif /* USE_EIGEN3 */

double norm2(const NPoint &vec)
{
  double s = 0.0;

  for (uint64_t i = 0; i < 3; ++i)
  {
    s += vec[i] * vec[i];
  }
  return sqrt(s);
}

double norm2(const std::vector<double> &vec)
{
  double s = 0.0;

  for (uint64_t i = 0; i < vec.size(); ++i)
  {
    s += vec[i] * vec[i];
  }
  return sqrt(s);
}

#ifdef USE_EIGEN3
int fullPiv(int n, double *A)
{
  Eigen::MatrixXd eA = Eigen::Map<Eigen::MatrixXd>(A, n, n);
/*  if (abs(eA.determinant()) < 1E-12) {
    std::cout << "Matrice non inversible. Determinant : " << eA.determinant() << std::endl;
    }*/
  eA.fullPivLu();
  Eigen::Map<Eigen::MatrixXd>(A, n, n) = eA.inverse();
  return 0;
}

int partialPiv(int n, double *A)
{
  Eigen::MatrixXd eA = Eigen::Map<Eigen::MatrixXd>(A, n, n);
  eA.partialPivLu();
  Eigen::Map<Eigen::MatrixXd>(A, n, n) = eA.inverse();
  return 0;
}

int eiMatVecProduct(int n,double *A,double *b)
{
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eA
    = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(A, n, n);
  Eigen::VectorXd eB = Eigen::Map<Eigen::VectorXd>(b, n);
  eB = eA * eB;
  Eigen::Map<Eigen::VectorXd>(b, n) = eB;
  return 0;
}

int inverse(int n, double *A)
{
  Eigen::MatrixXd eA = Eigen::Map<Eigen::MatrixXd>(A, n, n);
  Eigen::Map<Eigen::MatrixXd>(A, n, n) = eA.inverse();
  return 0;
}

int matMatProduct(int lda, int cda, int ldb, int cdb,
                  double *A, double *B, double *C,
                  bool transposeA, bool transposeB,
                  double alpha, double beta)
{
  Eigen::MatrixXd eA = Eigen::Map<Eigen::MatrixXd>(A, lda, cda);
  Eigen::MatrixXd eB = Eigen::Map<Eigen::MatrixXd>(B, ldb, cdb);
  Eigen::MatrixXd eC;
  int ldc(0), cdc(0);
  if ((transposeA) && (!transposeB) )
  {
    eC = alpha * eA.transpose() * eB;
    eC += beta*eC;
    ldc = cda;
    cdc = cdb;
  }
  else if ((!transposeA) && (transposeB))
  {
    eC = alpha * eA * eB.transpose();
    eC += beta * eC;
    ldc = lda;
    cdc = ldb;
  }
  else if ((transposeA) && (transposeB))
  {
    eC = alpha * eA.transpose() * eB.transpose();
    eC += beta*eC;
    ldc = cda;
    cdc = ldb;
  }
  else
  {
    eC = alpha * eA * eB;
    eC += beta*eC;
    ldc = lda;
    cdc = cdb;
  }
  Eigen::Map<Eigen::MatrixXd>(C, ldc, cdc) = eC;

  return 0;
}

int matVecProductNS(int lda, int cda,
                    double *A, double *X, double *Y,
                    bool transposeA,
                    double alpha, double beta,
                    int incX, int incY)
{
  Eigen::MatrixXd eA = Eigen::Map<Eigen::MatrixXd>(A, lda, cda);
  Eigen::VectorXd eX;
  Eigen::MatrixXd eY;
  int ldy(0);

  if (transposeA)
  {
    ldy = cda;
    eX  = Eigen::Map<Eigen::VectorXd>(X, lda);
    eY  = alpha * eA.transpose() * eX;
    eY += beta*eY;
  }
  else
  {
    ldy = lda;
    eX  = Eigen::Map<Eigen::VectorXd>(X, cda);
    eY  = alpha * eA * eX;
    eY += beta * eY;
  }

  Eigen::Map<Eigen::VectorXd>(Y, ldy) = eY;

  return 0;

}
#endif /* USE_EIGEN3 */

#ifdef USE_LAPACKE
int inverse(int n, double* A)
{
  int ipiv[n];
  int lwork = n;
  double work[lwork];
  int info;

  dgetrf_(&n, &n, A, &n, ipiv, &info);
  if (info != 0) {
    std::cout << "Matrix facto failed " << info <<" "<< n <<std::endl;
    for (int i = 0; i < n; ++i)
    {
      std::cout << "| ";
      for (int j = 0; j < n; ++j)
      {
        std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) <<
          A[i * n + j]  << "| ";
      }
      std::cout << std::endl;
    }
    exit(0);
  }

  dgetri_(&n, A, &n, ipiv, work, &lwork, &info);
  if (info != 0)
  {
    std::cout << "Matrix inversion failed " << info << " " << n << std::endl;
    for (int i = 0; i < n; ++i)
    {
      std::cout << "| ";
      for (int j = 0; j < n; ++j)
      {
        std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << A[i * n + j]  << "| ";
      }
      std::cout << std::endl;
    }
    exit(0);
  }
  return 0;
}

int matMatProduct(int lda, int cda, int ldb, int cdb,
                  double *A, double *B, double *C,
                  bool transposeA, bool transposeB,
                  double alpha, double beta)
{
  // FIXME: Add an assert macro here
  char transA('N');
  char transB('N');
  int M(lda), N(cdb), K(cda);

  if (transposeA)
  {
    transA = 'T';
    M = cda;
    K = lda;
  }
  if (transposeB)
  {
    transB = 'T';
    N = ldb;
  }

  dgemm_(&transA, &transB,
         &M, &N, &K,
         &alpha, A, &lda,
         B, &ldb, &beta,
         C, &M);

  return 0;

}

int matVecProductNS(int lda, int cda,
                    double *A, double *X, double *Y,
                    bool transposeA,
                    double alpha, double beta,
                    int incX, int incY)
{
  // FIXME: Add an assert macro here
  char transA('N');

  if (transposeA)
  {
    transA = 'T';
  }

  dgemv_(&transA,
         &lda, &cda,
         &alpha, A, &lda,
         X, &incX, &beta,
         Y, &incY);

  return 0;

}

#endif /* USE_LAPACKE */
int matVecProduct(int n,double *A,double *b)
{
  int i, j;
  double *x = (double*)calloc((size_t)(n), sizeof(double));

  for(i = 0; i < n; i++)
  {
    double sdum = 0;

    for(j = 0; j < n; j++)
    {
      sdum = sdum + A[i * n + j] * b[j];
    }
    x[i] = sdum;
  }
  for(i = 0; i < n; i++)
  {
    b[i] = x[i];
  }
  free(x);
  return 0;
}

int multTransposeLeft(int lda, int cda,
                      double *A, double *C)
{
  double* tA = new double[cda * lda]();
  std::copy(A, A + lda*cda, tA);

  matMatProduct(lda, cda, lda, cda,
                tA, A, C,
                true, false);

  delete [] tA;
  return 0;

}

