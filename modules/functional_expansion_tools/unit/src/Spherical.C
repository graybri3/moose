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

TEST(FunctionalExpansionsTest, sphericalHarmonicsConstructor)
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

 TEST(FunctionalExpansionsTest, sphericalHarmonicsSeriesOrthonormalEvaluation)
 {
   const std::vector<MooseEnum> domains = {FunctionalBasisInterface::_domain_options = "x",
                                           FunctionalBasisInterface::_domain_options = "y",
                                           FunctionalBasisInterface::_domain_options = "z"};
   const std::vector<std::size_t> orders = {4};
   std::vector<MooseEnum> series;
   expansion_type = "standard";
   generation_type = "orthonormal";
   series.push_back(single_series_types_3D = "SphericalHarmonics");
   Spherical spherical(domains, orders, series, name, expansion_type, generation_type);
   EXPECT_EQ(spherical.getNumberOfTerms(), (orders[0] + 1)*(orders[0] + 1));
 }

 TEST(FunctionalExpansionsTest, sphericalHarmonicsSeriesSqrtMuEvaluation)
 {
   const std::vector<MooseEnum> domains = {FunctionalBasisInterface::_domain_options = "x",
                                           FunctionalBasisInterface::_domain_options = "y",
                                           FunctionalBasisInterface::_domain_options = "z"};
   const std::vector<std::size_t> orders = {4};
   std::vector<MooseEnum> series;
   expansion_type = "standard";
   generation_type = "sqrt_mu";
   series.push_back(single_series_types_3D = "SphericalHarmonics");
   Spherical spherical(domains, orders, series, name, expansion_type, generation_type);
   EXPECT_EQ(spherical.getNumberOfTerms(), (orders[0] + 1)*(orders[0] + 1));
 }

 TEST(FunctionalExpansionsTest, sphericalHarmonicsSeriesStandardEvaluation)
 {
   const std::vector<MooseEnum> domains = {FunctionalBasisInterface::_domain_options = "x",
                                           FunctionalBasisInterface::_domain_options = "y",
                                           FunctionalBasisInterface::_domain_options = "z"};
   const std::vector<std::size_t> orders = {4};
   std::vector<MooseEnum> series;
   expansion_type = "standard";
   generation_type = "standard";
   series.push_back(single_series_types_3D = "SphericalHarmonics");
   Spherical spherical(domains, orders, series, name, expansion_type, generation_type);
   EXPECT_EQ(spherical.getNumberOfTerms(), (orders[0] + 1)*(orders[0] + 1));
 }

 TEST(FunctionalExpansionsTest, volumetricSphericalSeriesOrthonormalEvaluation)
 {
   const std::vector<MooseEnum> domains = {FunctionalBasisInterface::_domain_options = "x",
                                           FunctionalBasisInterface::_domain_options = "y",
                                           FunctionalBasisInterface::_domain_options = "z"};
   const std::vector<std::size_t> orders = {4};
   std::vector<MooseEnum> series;
   expansion_type = "standard";
   generation_type = "orthonormal";
   series.push_back(single_series_types_3D = "VolumetricSpherical");
   Spherical spherical(domains, orders, series, name, expansion_type, generation_type);
   EXPECT_EQ(spherical.getNumberOfTerms(), (orders[0] + 1)*(orders[0] + 2) / 2);
 }

 TEST(FunctionalExpansionsTest, volumetricSphericalSeriesSqrtMuEvaluation)
 {
   const std::vector<MooseEnum> domains = {FunctionalBasisInterface::_domain_options = "x",
                                           FunctionalBasisInterface::_domain_options = "y",
                                           FunctionalBasisInterface::_domain_options = "z"};
   const std::vector<std::size_t> orders = {4};
   std::vector<MooseEnum> series;
   expansion_type = "standard";
   generation_type = "sqrt_mu";
   series.push_back(single_series_types_3D = "VolumetricSpherical");
   Spherical spherical(domains, orders, series, name, expansion_type, generation_type);
   EXPECT_EQ(spherical.getNumberOfTerms(), (orders[0] + 1)*(orders[0] + 2) / 2);
 }

 TEST(FunctionalExpansionsTest, volumetricSphericalSeriesStandardEvaluation)
 {
   const std::vector<MooseEnum> domains = {FunctionalBasisInterface::_domain_options = "x",
                                           FunctionalBasisInterface::_domain_options = "y",
                                           FunctionalBasisInterface::_domain_options = "z"};
   const std::vector<std::size_t> orders = {4};
   std::vector<MooseEnum> series;
   expansion_type = "standard";
   generation_type = "standard";
   series.push_back(single_series_types_3D = "VolumetricSpherical");
   Spherical spherical(domains, orders, series, name, expansion_type, generation_type);
   EXPECT_EQ(spherical.getNumberOfTerms(), (orders[0] + 1)*(orders[0] + 2) / 2);
 }
//   const unsigned int order = 4;
//   SphericalHarmonics sphericalharmonics({FunctionalBasisInterface::_domain_options = "z",
//                                          FunctionalBasisInterface::_domain_options = "x",
//                                          FunctionalBasisInterface::_domain_options = "y"},
//                                          {order},
//                                          expansion_type = "standard",
//                                          generation_type = "orthonormal");
//   const Point location(-0.90922108754014, 0.262698547343, 0.156796889218);
//   sphericalharmonics.setLocation(Point(location));
//   const std::array<Real, 25> truth = {{0.28209479177388,
//                                       0.019813587525427386,
//                                       0.4718398868793876,
//                                       0.1253274013807911,
//                                       0.011364189298441322,
//                                       0.042784560550813006,
//                                       0.61455278861632501,
//                                       0.27062629552428014,
//                                       0.035042792743601552,
//                                       0.0046833871214201237,
//                                       0.029035307853304135,
//                                       0.067886244775306651,
//                                       0.59923369671782123,
//                                       0.42940263272720819,
//                                       0.089533731674951711,
//                                       0.0092109712880876129,
//                                       0.0016703447610675304,
//                                       0.042625564207503415,
//                                       0.17091573427046261,
//                                       0.29038447774788734,
//                                       0.57776583573696927,
//                                       0.58466411946433028,
//                                       0.16776154733386961,
//                                       0.026684903647102437,
//                                       0.0023045087667355533}};
//
// auto & answer = sphericalharmonics.getAllGeneration();
// for (std::size_t i = 0; i < sphericalharmonics.getNumberOfTerms(); i++)
//   EXPECT_NEAR(answer[i], truth[i], tol);
// }
//
// TEST(FunctionalExpansionsTest, sphericalHarmonicsSeriesSqrtMuEvaluation)
// {
//   const unsigned int order = 4;
//   const std::vector<MooseEnum> series;
//   series.push_back(single_series_types_3D = "SphericalHarmonics");
//   expansion_type = "standard";
//   generation_type = "sqrt_mu";
//   SphericalHarmonics sphericalharmonics({FunctionalBasisInterface::_domain_options = "x",
//                                          FunctionalBasisInterface::_domain_options = "y",
//                                          FunctionalBasisInterface::_domain_options = "z"},
//                                          {order},
//                                          series,
//                                          expansion_type,
//                                          generation_type);
//   const Point location(-0.90922108754014, 0.262698547343, 0.156796889218);
//   sphericalharmonics.setLocation(location);
//
//   const std::array<Real, 25> truth = {{0.28209479177388,
//                                        0.019813587525427386,
//                                        0.4718398868793876,
//                                        0.1253274013807911,
//                                        0.011364189298441322,
//                                        0.042784560550813006,
//                                        0.61455278861632501,
//                                        0.27062629552428014,
//                                        0.035042792743601552,
//                                        0.0046833871214201237,
//                                        0.029035307853304135,
//                                        0.067886244775306651,
//                                        0.59923369671782123,
//                                        0.42940263272720819,
//                                        0.089533731674951711,
//                                        0.0092109712880876129,
//                                        0.0016703447610675304,
//                                        0.042625564207503415,
//                                        0.17091573427046261,
//                                        0.29038447774788734,
//                                        0.57776583573696927,
//                                        0.58466411946433028,
//                                        0.16776154733386961,
//                                        0.026684903647102437,
//                                        0.0023045087667355533}};
//
// auto & answer = sphericalharmonics.getAllGeneration();
// for (std::size_t i = 0; i <= sphericalharmonics.getNumberOfTerms(); i++)
//   EXPECT_NEAR(answer[i], truth[i], tol);
// }
//
// TEST(FunctionalExpansionsTest, sphericalHarmonicsSeriesStandardEvaluation)
// {
//   const unsigned int order = 4;
//   SphericalHarmonics sphericalharmonics({FunctionalBasisInterface::_domain_options = "z",
//                                          FunctionalBasisInterface::_domain_options = "x",
//                                          FunctionalBasisInterface::_domain_options = "y"},
//                                          {order},
//                                          expansion_type = "standard",
//                                          generation_type = "standard");
//
//   const Point location(-0.90922108754014, 0.262698547343, 0.156796889218);
//   sphericalharmonics.setLocation(location);
//   const std::array<Real, 25> truth = {{0.28209479177388,
//                                        0.019813587525427386,
//                                        0.4718398868793876,
//                                        0.1253274013807911,
//                                        0.011364189298441322,
//                                        0.042784560550813006,
//                                        0.61455278861632501,
//                                        0.27062629552428014,
//                                        0.035042792743601552,
//                                        0.0046833871214201237,
//                                        0.029035307853304135,
//                                        0.067886244775306651,
//                                        0.59923369671782123,
//                                        0.42940263272720819,
//                                        0.089533731674951711,
//                                        0.0092109712880876129,
//                                        0.0016703447610675304,
//                                        0.042625564207503415,
//                                        0.17091573427046261,
//                                        0.29038447774788734,
//                                        0.57776583573696927,
//                                        0.58466411946433028,
//                                        0.16776154733386961,
//                                        0.026684903647102437,
//                                        0.0023045087667355533}};
//
// auto & answer = sphericalharmonics.getAllGeneration();
// for (std::size_t i = 0; i < sphericalharmonics.getNumberOfTerms(); i++)
//   EXPECT_NEAR(answer[i], truth[i], tol);
// }
// TEST(FunctionalExpansionsTest, volumetricsphericalSeriesEvaluation)
// {
//   const unsigned int order = 4;
//   const Point location(-0.90922108754014, 0.262698547343, 0.156796889218);
//   //const std::vector<Point> location = {Point(-0.90922108754014, 0.262698547343, 0.156796889218)};
//   expansion_type = "standard";
//   generation_type = "standard";
//
//   VolumetricSpherical volumetricspherical({FunctionalBasisInterface::_domain_options = "z",
//                                            FunctionalBasisInterface::_domain_options = "x",
//                                            FunctionalBasisInterface::_domain_options = "y"},
//                                            {order},
//                                            expansion_type,
//                                            generation_type);
//     volumetricspherical.setLocation(location);
//   const std::array<Real, 15> truth = {{1.000000000000000,
//                                         0.262698547343,
//                                        -0.90922108754014,
//                                        -0.23885105891041744,
//                                        0.791881706198328,
//                                        0.623918544603900,
//                                        1.365900806059890,
//                                        -1.357069098439213,
//                                        -0.288633095470502,
//                                        -1.073332058640328,
//                                        -1.349767446988918,
//                                        1.887803726543957,
//                                        0.404825692079851,
//                                        0.223343150486739,
//                                        0.698275682841869}};
//
//   auto & answer = volumetricspherical.getAllGeneration();
//   for (std::size_t i = 0; i < volumetricspherical.getNumberOfTerms(); i++)
//     EXPECT_NEAR(answer[i], truth[i], tol);
// }

TEST(FunctionalExpansionsTest, SphericalConstructor)
{
  const std::vector<MooseEnum> domains = {FunctionalBasisInterface::_domain_options = "x",
                                          FunctionalBasisInterface::_domain_options = "y",
                                          FunctionalBasisInterface::_domain_options = "z"};
  const std::vector<std::size_t> orders = {4};
  std::vector<MooseEnum> series;
  expansion_type = "standard";
  generation_type = "orthonormal";
  series.push_back(single_series_types_3D = "SphericalHarmonics");
  Spherical spherical(domains, orders, series, name, expansion_type, generation_type);
  EXPECT_EQ(spherical.getNumberOfTerms(), (orders[0] + 1)*(orders[0] + 1));
}
