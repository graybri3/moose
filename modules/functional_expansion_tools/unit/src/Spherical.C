//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include <array>

#include "Setup.h"

TEST(FunctionalExpansionTest, sphericalHarmonicsConstructor)
{
  const unsigned int order = 4;

  SphericalHarmonics sphericalharmonics({FunctionalBasisInterface::_domain_options = "x",
                                         FunctionalBasisInterface::_domain_options = "y",
                                         FunctionalBasisInterface::_domain_options = "z"},
                                        {order},
                                        expansion_type = "orthonormal",
                                        generation_type = "standard");
  EXPECT_EQ(sphericalharmonics.getOrder(0), order);
}

TEST(FunctionalExpansionTest, sphericalHarmonicsOrthonormalEvaluation)
{
  const unsigned int order = 4;
  const Point location(-0.90922108754014, 0.262698547343, 0.156796889218);
  expansion_type = "standard";
  generation_type = "orthonormal";
  SphericalHarmonics sphericalharmonics({FunctionalBasisInterface::_domain_options = "x",
                                         FunctionalBasisInterface::_domain_options = "y",
                                         FunctionalBasisInterface::_domain_options = "z"},
                                         {order},
                                         expansion_type,
                                         generation_type);
  sphericalharmonics.setLocation(location);
  const std::array<Real, 25> truth = {{0.28209479177388,
                                      0.019813587525427386,
                                      0.4718398868793876,
                                      0.1253274013807911,
                                      0.011364189298441322,
                                      0.042784560550813006,
                                      0.61455278861632501,
                                      0.27062629552428014,
                                      0.035042792743601552,
                                      0.0046833871214201237,
                                      0.029035307853304135,
                                      0.067886244775306651,
                                      0.59923369671782123,
                                      0.42940263272720819,
                                      0.089533731674951711,
                                      0.0092109712880876129,
                                      0.0016703447610675304,
                                      0.042625564207503415,
                                      0.17091573427046261,
                                      0.29038447774788734,
                                      0.57776583573696927,
                                      0.58466411946433028,
                                      0.16776154733386961,
                                      0.026684903647102437,
                                      0.0023045087667355533}};

auto & answer = sphericalharmonics.getAllGeneration();
for (std::size_t i = 0; i < sphericalharmonics.getNumberOfTerms(); ++i)
  EXPECT_NEAR(answer[i], truth[i], tol);
}

TEST(FunctionalExpansionTest, sphericalHarmonicsSeriesSqrtMuEvaluation)
{
  const unsigned int order = 4;
  const Point location(-0.90922108754014, 0.262698547343, 0.156796889218);
  SphericalHarmonics sphericalharmonics({FunctionalBasisInterface::_domain_options = "x",
                                         FunctionalBasisInterface::_domain_options = "y",
                                         FunctionalBasisInterface::_domain_options = "z"},
                                         {order},
                                         expansion_type,
                                         generation_type);
  sphericalharmonics.setLocation(location);
  const std::array<Real, 25> truth = {{0.28209479177388,
                                       0.019813587525427386,
                                       0.4718398868793876,
                                       0.1253274013807911,
                                       0.011364189298441322,
                                       0.042784560550813006,
                                       0.61455278861632501,
                                       0.27062629552428014,
                                       0.035042792743601552,
                                       0.0046833871214201237,
                                       0.029035307853304135,
                                       0.067886244775306651,
                                       0.59923369671782123,
                                       0.42940263272720819,
                                       0.089533731674951711,
                                       0.0092109712880876129,
                                       0.0016703447610675304,
                                       0.042625564207503415,
                                       0.17091573427046261,
                                       0.29038447774788734,
                                       0.57776583573696927,
                                       0.58466411946433028,
                                       0.16776154733386961,
                                       0.026684903647102437,
                                       0.0023045087667355533}};


auto & answer = sphericalharmonics.getAllGeneration();
for (std::size_t i = 0; i < sphericalharmonics.getNumberOfTerms(); ++i)
  EXPECT_NEAR(answer[i], truth[i], tol);
}

TEST(FunctionalExpansionTest, sphericalHarmonicsSeriesStandardEvaluation)
{
  const unsigned int order = 4;
  const Point location(-0.90922108754014, 0.262698547343, 0.156796889218);
  SphericalHarmonics sphericalharmonics({FunctionalBasisInterface::_domain_options = "x",
                                         FunctionalBasisInterface::_domain_options = "y",
                                         FunctionalBasisInterface::_domain_options = "z"},
                                         {order},
                                         expansion_type,
                                         generation_type);
  sphericalharmonics.setLocation(location);
  const std::array<Real, 25> truth = {{0.28209479177388,
                                       0.019813587525427386,
                                       0.4718398868793876,
                                       0.1253274013807911,
                                       0.011364189298441322,
                                       0.042784560550813006,
                                       0.61455278861632501,
                                       0.27062629552428014,
                                       0.035042792743601552,
                                       0.0046833871214201237,
                                       0.029035307853304135,
                                       0.067886244775306651,
                                       0.59923369671782123,
                                       0.42940263272720819,
                                       0.089533731674951711,
                                       0.0092109712880876129,
                                       0.0016703447610675304,
                                       0.042625564207503415,
                                       0.17091573427046261,
                                       0.29038447774788734,
                                       0.57776583573696927,
                                       0.58466411946433028,
                                       0.16776154733386961,
                                       0.026684903647102437,
                                       0.0023045087667355533}};

auto & answer = sphericalharmonics.getAllGeneration();
for (std::size_t i = 0; i < sphericalharmonics.getNumberOfTerms(); ++i)
{
  EXPECT_NEAR(answer[i], truth[i], tol);
}
}

// TEST(FunctionalExpansionsTest, SphericalConstructor)
// {
//
//   domains.push_back(FunctionalBasisInterface::_domain_options = "sphere");
//   orders = {4};
//   series.push_back(single_series_types_3D = "SphericalHarmonics");
//   Spherical spherical(domains, orders, series, name, expansion_type, generation_type);
//   EXPECT_EQ(spherical.getNumberOfTerms(), (orders[0] + 1)*(orders[0] + 1));
// }
