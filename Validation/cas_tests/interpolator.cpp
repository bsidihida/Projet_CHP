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

#include "InterpolatorFactory.hpp"
#include "Function.hpp"
#include "common.hpp"
#include "gtest/gtest.h"
#include <math.h>

using namespace neos;

TEST(interpolator2D, 3P1)
{
  std::array<double, 3> xpos { { 5.0/3.0, 1.25} };
  std::array<double, 3> q1 { { 1.0, 1.0} };
  std::array<double, 3> q2 { { 2.0, 1.0} };
  std::array<double, 3> q3 { { 2.0, 2.0} };
  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  std::vector<double> val_ref;
  val_ref.push_back(2.0);
  val_ref.push_back(2.0);
  val_ref.push_back(1.0);

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::POLY2D);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);
  EXPECT_NEAR(res, 1.75, 1e-5);

  std::vector<double> weight = interpo->computeInterpolation(xpos, xref);
  double res2 = weight[0] * val_ref[0] + weight[1] * val_ref[1] + val_ref[2] * weight[2];
  EXPECT_NEAR(res2, 1.75, 1e-5);
  delete interpo;
}

TEST(interpolator2D, 3P2)
{
  std::array<double, 3> xpos { { 5.0/3.0, 1.25} };
  std::array<double, 3> q1 { { 1.0, 1.0} };
  std::array<double, 3> q2 { { 2.0, 1.0} };
  std::array<double, 3> q3 { { 2.0, 2.0} };
  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  std::vector<double> val_ref;
  val_ref.push_back(2.25);
  val_ref.push_back(1.23);
  val_ref.push_back(0.45);

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::POLY2D);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);
  EXPECT_NEAR(res, 1.375, 1e-5);

  std::vector<double> weight = interpo->computeInterpolation(xpos, xref);
  double res2 = weight[0] * val_ref[0] + weight[1] * val_ref[1] + val_ref[2] * weight[2];
  EXPECT_NEAR(res2, 1.375, 1e-5);
  delete interpo;
}

TEST(interpolator2D, 3P3)
{
  std::array<double, 3> xpos { { 1.0, 1.0} };
  std::array<double, 3> q1 { { 0.5, 0.25} };
  std::array<double, 3> q2 { { 0.75, 1.75} };
  std::array<double, 3> q3 { { 2.0, 1.0} };
  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  std::vector<double> val_ref;
  val_ref.push_back(-2.25);
  val_ref.push_back(1.23);
  val_ref.push_back(0.45);

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::POLY2D);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);
  EXPECT_NEAR(res, -0.248178, 1e-5);

  std::vector<double> weight = interpo->computeInterpolation(xpos, xref);
  double res2 = weight[0] * val_ref[0] + weight[1] * val_ref[1] + val_ref[2] * weight[2];
  EXPECT_NEAR(res2, -0.248178, 1e-5);
  delete interpo;
}


TEST(interpolator2D, 4P1)
{
  std::array<double, 3> xpos { { 5.0/3.0, 1.25} };
  std::array<double, 3> q1 { { 1.0, 1.0} };
  std::array<double, 3> q2 { { 2.0, 1.0} };
  std::array<double, 3> q3 { { 2.0, 2.0} };
  std::array<double, 3> q4 { { 1.0, 2.0} };
  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  xref.push_back(q4);
  std::vector<double> val_ref;
  val_ref.push_back(2.0);
  val_ref.push_back(2.0);
  val_ref.push_back(1.0);
  val_ref.push_back(1.0);

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::POLY2D);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);
  EXPECT_NEAR(res, 1.75, 1e-8);

  std::vector<double> weight = interpo->computeInterpolation(xpos, xref);
  double res2 = weight[0] * val_ref[0] + weight[1] * val_ref[1] + val_ref[2] * weight[2] + val_ref[3] * weight[3];
  EXPECT_NEAR(res2, 1.75, 1e-8);
  delete interpo;
}

TEST(interpolator2D, 4P2)
{
  std::array<double, 3> xpos { { 5.0/3.0, 1.25} };
  std::array<double, 3> q1 { { 1.0, 1.0} };
  std::array<double, 3> q2 { { 2.0, 1.0} };
  std::array<double, 3> q3 { { 2.0, 2.0} };
  std::array<double, 3> q4 { { 1.0, 2.0} };
  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  xref.push_back(q4);
  std::vector<double> val_ref;
  val_ref.push_back(0.23);
  val_ref.push_back(3.46);
  val_ref.push_back(-2.67);
  val_ref.push_back(-0.1);

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::POLY2D);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);
  EXPECT_NEAR(res, 1.334166666, 1e-8);

  std::vector<double> weight = interpo->computeInterpolation(xpos, xref);
  double res2 = weight[0] * val_ref[0] + weight[1] * val_ref[1] + val_ref[2] * weight[2] + val_ref[3] * weight[3];
  EXPECT_NEAR(res2, 1.334166666, 1e-8);
  delete interpo;
}

TEST(interpolator2D, 4P3)
{
  std::array<double, 3> xpos { { 1.5, 1} };
  std::array<double, 3> q1 { { 0.5, 0.25} };
  std::array<double, 3> q2 { { 0.75, 1.75} };
  std::array<double, 3> q3 { { 2.0, 1.0} };
  std::array<double, 3> q4 { { 2, 1.5} };
  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  xref.push_back(q4);
  std::vector<double> val_ref;
  val_ref.push_back(0.23);
  val_ref.push_back(3.46);
  val_ref.push_back(-2.67);
  val_ref.push_back(-0.1);

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::POLY2D);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);
  EXPECT_NEAR(res, -1.08916666668, 1e-8);


  std::vector<double> weight = interpo->computeInterpolation(xpos, xref);
  double res2 = weight[0] * val_ref[0] + weight[1] * val_ref[1] + val_ref[2] * weight[2] + val_ref[3] * weight[3];
  EXPECT_NEAR(res2, -1.08916666668, 1e-8);
  delete interpo;
}

TEST(interpolator3D, 8P1)
{
  std::array<double, 3> xpos { { 0.75, 0.75, 0.75} };
  std::array<double, 3> q1 { { 0.5, 0.5, 0.5} };
  std::array<double, 3> q2 { { 0.5, 0.5, 1} };
  std::array<double, 3> q3 { { 0.5, 1, 0.5} };
  std::array<double, 3> q4 { { 0.5, 1, 1} };
  std::array<double, 3> q5 { { 1, 1, 0.5} };
  std::array<double, 3> q6 { { 1, 1, 1} };
  std::array<double, 3> q7 { { 1, 0.5, 0.5} };
  std::array<double, 3> q8 { { 1, 0.5, 1} };

  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  xref.push_back(q4);
  xref.push_back(q5);
  xref.push_back(q6);
  xref.push_back(q7);
  xref.push_back(q8);

  std::vector<double> val_ref;
  val_ref.push_back(1);
  val_ref.push_back(1.75);
  val_ref.push_back(3);
  val_ref.push_back(0.25);
  val_ref.push_back(2.763);
  val_ref.push_back(3.3336);
  val_ref.push_back(4.086);
  val_ref.push_back(1.253);

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::POLY3D);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);

  EXPECT_NEAR(res, 2.1795, 1e-4);
  delete interpo;
}


TEST(interpolator3D, 7P1)
{
  std::array<double, 3> xpos { { 0.75, 0.75, 0.75} };
  std::array<double, 3> q1 { { 0.5, 0.5, 0.5} };
  std::array<double, 3> q2 { { 0.5, 0.5, 1} };
  std::array<double, 3> q3 { { 0.5, 1, 0.5} };
  std::array<double, 3> q4 { { 0.5, 1, 1} };
  std::array<double, 3> q5 { { 1, 1, 0.5} };
  std::array<double, 3> q6 { { 1, 1, 1} };
  std::array<double, 3> q7 { { 1, 0.5, 0.5} };

  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  xref.push_back(q4);
  xref.push_back(q5);
  xref.push_back(q6);
  xref.push_back(q7);

  std::vector<double> val_ref;
  val_ref.push_back(1.0);
  val_ref.push_back(1.75);
  val_ref.push_back(3.0);
  val_ref.push_back(0.25);
  val_ref.push_back(2.763);
  val_ref.push_back(3.3336);
  val_ref.push_back(4.086);

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::POLY3D);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);

  EXPECT_NEAR(res, 3.0424, 1e-4);
  delete interpo;
}

TEST(interpolatorRBF, RBF1)
{
  std::array<double, 3> xpos { { 0.75, 0.75, 0.75} };
  std::array<double, 3> q1 { { 0.5, 0.5, 0.5} };
  std::array<double, 3> q2 { { 0.5, 0.5, 1} };
  std::array<double, 3> q3 { { 0.5, 1, 0.5} };
  std::array<double, 3> q4 { { 0.5, 1, 1} };
  std::array<double, 3> q5 { { 1, 1, 0.5} };
  std::array<double, 3> q6 { { 1, 1, 1} };
  std::array<double, 3> q7 { { 1, 0.5, 0.5} };

  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  xref.push_back(q4);
  xref.push_back(q5);
  xref.push_back(q6);
  xref.push_back(q7);

  std::vector<double> val_ref;
  val_ref.push_back(1.0);
  val_ref.push_back(1.75);
  val_ref.push_back(3.0);
  val_ref.push_back(0.25);
  val_ref.push_back(2.763);
  val_ref.push_back(3.3336);
  val_ref.push_back(4.086);

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::RBF);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);

  EXPECT_NEAR(res, 3.0424, 1e-4);
  delete interpo;
}

TEST(interpolatorRBF, RBF3)
{
  std::array<double, 3> xpos { { 5.0/3.0, 1.25} };
  std::array<double, 3> q1 { { 1.0, 1.0} };
  std::array<double, 3> q2 { { 2.0, 1.0} };
  std::array<double, 3> q3 { { 2.0, 2.0} };
  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  std::vector<double> val_ref;
  val_ref.push_back(2.0);
  val_ref.push_back(2.0);
  val_ref.push_back(1.0);


  IInterpolator *interpo = InterpolatorFactory::get(interpoType::RBF);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);

  EXPECT_NEAR(res, 1.75, 1e-5);
  delete interpo;
}

TEST(interpolatorRBF, RBF4)
{
  std::array<double, 3> xpos { { 5.0/3.0, 1.25} };
  std::array<double, 3> q1 { { 1.0, 1.0} };
  std::array<double, 3> q2 { { 2.0, 1.0} };
  std::array<double, 3> q3 { { 2.0, 2.0} };
  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  std::vector<double> val_ref;
  val_ref.push_back(2.25);
  val_ref.push_back(1.23);
  val_ref.push_back(0.45);


  IInterpolator *interpo = InterpolatorFactory::get(interpoType::RBF);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);

  EXPECT_NEAR(res, 1.375, 1e-5);
  delete interpo;
}

TEST(interpolatorRBF, RBF5)
{
  std::array<double, 3> xpos { { 1.0, 1.0} };
  std::array<double, 3> q1 { { 0.5, 0.25} };
  std::array<double, 3> q2 { { 0.75, 1.75} };
  std::array<double, 3> q3 { { 2.0, 1.0} };
  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  std::vector<double> val_ref;
  val_ref.push_back(-2.25);
  val_ref.push_back(1.23);
  val_ref.push_back(0.45);


  IInterpolator *interpo = InterpolatorFactory::get(interpoType::RBF);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);

  EXPECT_NEAR(res, -0.248178, 1e-5);
  delete interpo;
}


TEST(interpolatorRBF, RBF6)
{
  std::array<double, 3> xpos { { 5.0/3.0, 1.25} };
  std::array<double, 3> q1 { { 1.0, 1.0} };
  std::array<double, 3> q2 { { 2.0, 1.0} };
  std::array<double, 3> q3 { { 2.0, 2.0} };
  std::array<double, 3> q4 { { 1.0, 2.0} };
  std::vector<std::array<double, 3> > xref;
  xref.push_back(q1);
  xref.push_back(q2);
  xref.push_back(q3);
  xref.push_back(q4);
  std::vector<double> val_ref;
  val_ref.push_back(2.0);
  val_ref.push_back(2.0);
  val_ref.push_back(1.0);
  val_ref.push_back(1.0);


  IInterpolator *interpo = InterpolatorFactory::get(interpoType::RBF);
  double res = interpo->computeInterpolation(xpos, xref, val_ref);

  EXPECT_NEAR(res, 1.75, 1e-8);
  delete interpo;
}

TEST(interpolatorMLS, MLSVal)
{
  NPoint xpos = {0.,0.,0.};
  std::vector<NPoint> xref(8);

  xref[0] = { 1./10,  0, 0};
  xref[1] = { 1./10,  1./10, 0};
  xref[2] = { 0,  1./10, 0};
  xref[3] = {-1./10,  1./10, 0};
  xref[4] = {-1./10,  0, 0};
  xref[5] = {-1./10, -1./10, 0};
  xref[6] = { 0, -1./10, 0};
  xref[7] = { 1./10, -1./10, 0};

  std::vector<double> vref(8);
  size_t i = 0;
  for (auto &v: xref)
  {
    vref[i] = tan(v[0] + M_PI/4.);
    i++;
  }

  IInterpolator *interpo = InterpolatorFactory::get(interpoType::MLS);
  double res = interpo->computeInterpolation(xpos, xref, vref);

  EXPECT_NEAR(res, 1., 1e-3);
}
