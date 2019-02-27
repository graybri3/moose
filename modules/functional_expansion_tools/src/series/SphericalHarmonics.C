//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseUtils.h"
#include "SphericalHarmonics.h"
#include <functional>

#define MAX_DIRECT_CALCULATION_SPHERICALHARMONICS 4

SphericalHarmonics::SphericalHarmonics() : SingleSeriesBasisInterface() {}

SphericalHarmonics::SphericalHarmonics(const std::vector<MooseEnum> & domain,
                                       const std::vector<std::size_t> & order,
                                       MooseEnum expansion_type,
                                       MooseEnum generation_type)
      : SingleSeriesBasisInterface(domain, order, calculatedNumberOfTermsBasedOnOrder(order))
{
  if (expansion_type == "orthonormal")
    _evaluateExpansionWrapper = [this]() { this->evaluateOrthonormal(); };
  else if (expansion_type == "sqrt_mu")
    _evaluateExpansionWrapper = [this]() { this->evaluateSqrtMu(); };
  else if (expansion_type == "standard")
    _evaluateExpansionWrapper = [this]() { this->evaluateStandard(); };
  else
    mooseError("The specified type of normalization for expansion does not exist");

  if (generation_type == "orthonormal")
    _evaluateGenerationWrapper = [this]() { this->evaluateOrthonormal(); };
  else if (generation_type == "sqrt_mu")
    _evaluateGenerationWrapper = [this]() { this->evaluateSqrtMu(); };
  else if (generation_type == "standard")
    _evaluateGenerationWrapper = [this]() { this->evaluateStandard(); };
  else
    mooseError("The specified type of normalization for generation does not exist");
}

std::size_t
SphericalHarmonics::calculatedNumberOfTermsBasedOnOrder(const std::vector<std::size_t> & order) const
{
  return (order[0]+1)*(order[0]+1);
}

void
SphericalHarmonics::checkPhysicalBounds(const std::vector<Real> & bounds) const
{
  if (bounds.size() != 4)
      mooseError("Spherical Harmonics: Invalid number of bounds");
}
void
SphericalHarmonics::evaluateOrthonormal()
{
  const Real & theta = _standardized_location[1];
  const Real & phi = _standardized_location[2];
  const Real sin2 = sin(theta)*sin(theta);
  const Real sin4 = sin(theta)*sin(theta)*sin(theta)*sin(theta);//sin2*sin2;
  const Real cos2 = cos(theta)*cos(theta);
  const Real cos4 = cos2*cos2;
  switch (_orders[0])
  {
    default:
    case MAX_DIRECT_CALCULATION_SPHERICALHARMONICS: /*4*/
    save(24, 1.5*std::sqrt(0.546875/M_PI)*sin(theta)*sin(theta)*sin(theta)*sin(theta)*cos(4.0*phi));
    save(23, 1.5*std::sqrt(4.375/M_PI)*cos(theta)*sin(theta)*sin(theta)*sin(theta)*cos(3.0*phi));
    save(22, 1.5*std::sqrt(0.3125/M_PI)*sin(theta)*sin(theta)*(7.0*cos2-1.0)*cos(2.0*phi));
    save(21, 3.75*std::sqrt(0.1/M_PI)*cos(phi)*sin(theta)*(7.0*cos2*cos(theta)-3.0*cos(theta)));
    save(20, 0.1875/std::sqrt(M_PI)*(35.0*cos(theta)*cos(theta)*cos(theta)*cos(theta) - 30.0*cos(theta)*cos(theta) + 3.0));
    save(19, 3.75*std::sqrt(0.1*M_PI)*sin(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
    save(18, 1.5*std::sqrt(0.3125*M_PI)*sin(2*phi)*sin2*(7*cos2 -1));
    save(17, 1.5*std::sqrt(4.375*M_PI)*cos(theta)*sin2*sin(theta)*sin(3*phi));
    save(16, 1.5*std::sqrt(0.546875/M_PI)*sin4*sin(4*phi));
    libmesh_fallthrough();

    case 3:
    save(15, std::sqrt(1.09375/M_PI)*sin(theta)*sin(theta)*sin(theta)*cos(3.0*phi));
    save(14, std::sqrt(6.5625/M_PI)*cos(theta)*sin2*cos(2.0*phi));
    save(13, std::sqrt(0.65625/M_PI)*(5.0*cos2-1.0)*sin(theta)*cos(phi));
    save(12, std::sqrt(0.4375/M_PI)*(5.0*cos(theta)*cos(theta)*cos(theta)-3.0*cos(theta)));
    save(11, std::sqrt(0.65625/M_PI)*(5.0*cos(theta)*cos(theta)-1.0)*sin(theta)*sin(phi));
    save(10, std::sqrt(6.5625/M_PI)*sin(2.0*phi)*cos(theta)*sin2);
    save(9, std::sqrt(1.09375/M_PI)*sin(3.0*phi)*sin(theta)*sin(theta)*sin(theta));
    libmesh_fallthrough();

    case 2:
    save(8, std::sqrt(0.9375/M_PI)*sin(theta)*sin(theta)*cos(2.0*phi));
    save(7, std::sqrt(3.75/M_PI)*cos(phi)*sin(theta)*cos(theta));
    save(6, std::sqrt(0.3125/M_PI)*(0.5+1.5*cos(theta)));
    save(5, std::sqrt(3.75/M_PI)*sin(phi)*sin(theta)*cos(theta));
    save(4, std::sqrt(0.9375/M_PI)*sin(2*phi)*sin2);
    libmesh_fallthrough();

    case 1:
    save(3, std::sqrt(0.75/M_PI)*cos(phi)*sin(theta));
    save(2, std::sqrt(0.75/M_PI)*cos(theta));
    save(1, std::sqrt(0.75/M_PI)*sin(phi)*sin(theta));
    libmesh_fallthrough();

    case 0:
    save(0, std::sqrt(1/(4*M_PI)));
  }
  // const Real & xn = _standardized_location[1];
  // const Real & yn = _standardized_location[2];
  // const Real & zn = _standardized_location[3];
  // const Real & rho = _standardized_location[0];
  // const Real & rho2 = rho*rho;
  // const Real & rho4 = rho2*rho2;
  // switch (_orders[0])
  // {
  //   default:
  //   case MAX_DIRECT_CALCULATION_SPHERICALHARMONICS: /*2*/
  //   save(24, (xn*xn*(xn*xn - 3.0*yn*yn) - yn*yn*(3.0*xn*xn - yn*yn))/(rho4));//sin4*cos(4*phi));
  //   save(23, (xn*xn-3*yn*yn)*xn*zn/(rho4));//cos(theta)*sin2*sin(theta)*cos(3*phi));
  //   save(22, ((xn*xn-yn*yn)*(7*zn*zn-rho2))/(rho4));//sin2*(7*cos2-1)*cos(2*phi));
  //   save(21, (xn*zn*(7*zn*zn - 3*rho2))/(rho4));//cos(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
  //   save(20, (35*zn*zn*zn*zn - 30*zn*zn*rho2 + 3*rho4)/(rho4));//(35*cos4-30*cos2 + 3));
  //   save(19, (yn*zn*(3*xn*xn - yn*yn))/(rho4));//sin(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
  //   save(18, (xn*yn*(7*zn*zn - 3*rho2))/(rho4));//sin(2*phi)*sin2*(7*cos2 -1));
  //   save(17, (3*xn*xn - yn*yn)*(yn*zn)/(rho4));//cos(theta)*sin2*sin(theta)*sin(3*phi));
  //   save(16, (xn*yn*(xn*xn - yn*yn))/(rho4));//sin4*sin(4*phi));
  //   libmesh_fallthrough();
  //
  //   case 3:
  //   save(15, (xn*(xn*xn - 3.0*yn*yn))/(rho2*rho));//sin2*sin(theta)*cos(3*phi));
  //   save(14, (zn*(xn*xn - yn*yn))/(rho2*rho));//cos(theta)*sin2*cos(2*phi));
  //   save(13, (xn*(4*zn*zn - xn*xn - yn*yn))/(rho2*rho));//(5*cos2-1)*sin(theta)*cos(phi));
  //   save(12, zn*(2*zn*zn - 3*xn*xn - 3*yn*yn)/(rho2*rho));//(5*cos2*cos(theta)-3*cos(theta)));
  //   save(11, (yn*(4*zn*zn - xn*xn - yn*yn))/(rho2*rho));//(5*cos2-1)*sin(theta)*sin(phi));
  //   save(10, (xn*yn*zn)/(rho2*rho));//sin(2*phi)*cos(theta)*sin2);
  //   save(9, ((3*xn*xn - yn*yn)*yn)/(rho2*rho));//sin(3*phi)*sin2*sin(theta));
  //   libmesh_fallthrough();
  //
  //   case 2:
  //   save(8, (xn*xn-yn*yn)/(rho2));//sin2*cos(2.0*phi));
  //   save(7, zn*xn/(rho2));//cos(phi)*sin(theta)*cos(theta));
  //   save(6, (-xn*xn - yn*yn + 2.0*zn*zn)/(rho2));
  //   save(5, yn*zn/(rho2));
  //   save(4, xn*yn/(rho2));
  //   libmesh_fallthrough();
  //
  //   case 1:
  //   save(3, xn/rho);
  //   save(2, zn/rho);
  //   save(1, yn/rho);
  //   libmesh_fallthrough();
  //
  //   case 0:
  //   save(0, 1);
  // }
}

void
SphericalHarmonics::evaluateSqrtMu()
{
  const Real & theta = _standardized_location[1];
  const Real & phi = _standardized_location[2];
  const Real sin2 = sin(theta)*sin(theta);
  const Real sin4 = sin2*sin2;
  const Real cos2 = cos(theta)*cos(theta);
  const Real cos4 = cos2*cos2;
  switch (_orders[0])
  {
    default:
    case MAX_DIRECT_CALCULATION_SPHERICALHARMONICS: /*4*/
    save(24, 1.5*std::sqrt(35/(64*M_PI))*sin4*cos(4*phi));
    save(23, 1.5*std::sqrt(35/(8*M_PI))*cos(theta)*sin2*sin(theta)*cos(3*phi));
    save(22, 1.5*(std::sqrt(5/(16*M_PI)))*sin2*(7*cos2-1)*cos(2*phi));
    save(21, 15/4*std::sqrt(1/(10*M_PI))*cos(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
    save(20, 3/(16*std::sqrt(M_PI))*(35*cos4-30*cos2 + 3.0));
    save(19, 15/4*std::sqrt(1/(10*M_PI))*sin(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
    save(18, 3/2*std::sqrt(5/(16*M_PI))*sin(2*phi)*sin2*(7*cos2 -1));
    save(17, 3/2*std::sqrt(35/(8*M_PI))*cos(theta)*sin2*sin(theta)*sin(3*phi));
    save(16, 1.5*std::sqrt(35/(64*M_PI))*sin4*sin(4*phi));
    libmesh_fallthrough();

    case 3:
    save(15, std::sqrt(32/(35*M_PI))*sin2*sin(theta)*cos(3*phi));
    save(14, std::sqrt(105/(16*M_PI))*cos(theta)*sin2*cos(2*phi));
    save(13, std::sqrt(21/(32*M_PI))*(5*cos2-1.0)*sin(theta)*cos(phi));
    save(12, std::sqrt(7/(16*M_PI))*(5*cos2*cos(theta)-3*cos(theta)));
    save(11, std::sqrt(21/(16*M_PI))*(5*cos2-1.0)*sin(theta)*sin(phi));
    save(10, std::sqrt(105/(16*M_PI))*sin(2*phi)*cos(theta)*sin2);
    save(9, std::sqrt(35/(32*M_PI))*sin(3*phi)*sin2*sin(theta));
    libmesh_fallthrough();

    case 2:
    save(8, std::sqrt(15/(16*M_PI))*sin2*cos(2*phi));
    save(7, std::sqrt(15/(4*M_PI))*cos(phi)*sin(theta)*cos(theta));
    save(6, std::sqrt(5/(16*M_PI))*(3*cos2-1));
    save(5, std::sqrt(15/(4*M_PI))*sin(phi)*sin(theta)*cos(theta));
    save(4, std::sqrt(15/(4*M_PI))*sin(2*phi)*sin2);
    libmesh_fallthrough();

    case 1:
    save(3, std::sqrt(3/(4*M_PI))*cos(phi)*sin(theta));
    save(2, std::sqrt(3/(4*M_PI))*cos(theta));
    save(1, std::sqrt(3/(4*M_PI))*sin(phi)*sin(theta));
    libmesh_fallthrough();

    case 0:
    save(0, std::sqrt(1/(4*M_PI)));
  }
  // const Real & x = _standardized_location[1];
  // const Real & y = _standardized_location[2];
  // const Real & z = _standardized_location[3];
  // const Real & r = _standardized_location[0];
  // switch (_orders[0])
  // {
  //   default:
  //   case MAX_DIRECT_CALCULATION_SPHERICALHARMONICS: /*2*/
  //   save(24, (x*x*(x*x - 3.0*y*y) - y*y*(3.0*x*x - y*y))/(r*r*r*r));//sin4*cos(4*phi));
  //   save(23, (x*x-3*y*y)*x*z/(r*r*r*r));//cos(theta)*sin2*sin(theta)*cos(3*phi));
  //   save(22, ((x*x-y*y)*(7*z*z-r*r))/(r*r*r*r));//sin2*(7*cos2-1)*cos(2*phi));
  //   save(21, (x*z*(7*z*z - 3*r*r))/(r*r*r*r));//cos(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
  //   save(20, (35*z*z*z*z - 30*z*z*r*r + 3*r*r*r*r)/(r*r*r*r));//(35*cos4-30*cos2 + 3));
  //   save(19, (y*z*(3*x*x - y*y))/(r*r*r*r));//sin(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
  //   save(18, (x*y*(7*z*z - 3*r*r))/(r*r*r*r));//sin(2*phi)*sin2*(7*cos2 -1));
  //   save(17, (3*x*x - y*y)*(y*z)/(r*r*r*r));//cos(theta)*sin2*sin(theta)*sin(3*phi));
  //   save(16, (x*y*(x*x - y*y))/(r*r));//sin4*sin(4*phi));
  //   libmesh_fallthrough();
  //
  //   case 3:
  //   save(15, (x*(x*x - 3.0*y*y))/(r*r*r));//sin2*sin(theta)*cos(3*phi));
  //   save(14, (z*(x*x - y*y))/(r*r*r));//cos(theta)*sin2*cos(2*phi));
  //   save(13, (x*(4*z*z - x*x - y*y))/(r*r*r));//(5*cos2-1)*sin(theta)*cos(phi));
  //   save(12, z*(2*z*z - 3*x*x - 3*y*y)/(r*r*r));//(5*cos2*cos(theta)-3*cos(theta)));
  //   save(11, (y*(4*z*z - x*x - y*y))/(r*r*r));//(5*cos2-1)*sin(theta)*sin(phi));
  //   save(10, (x*y*z)/(r*r*r));//sin(2*phi)*cos(theta)*sin2);
  //   save(9, ((3*x*x - y*y)*y)/(r*r*r));//sin(3*phi)*sin2*sin(theta));
  //   libmesh_fallthrough();
  //
  //   case 2:
  //   save(8, (x*x-y*y)/(r*r));//sin2*cos(2.0*phi));
  //   save(7, z*x/(r*r));//cos(phi)*sin(theta)*cos(theta));
  //   save(6, (-x*x - y*y + 2.0*z*z)/(r*r));
  //   save(5, y*z/(r*r));
  //   save(4, x*y/(r*r));
  //   libmesh_fallthrough();
  //
  //   case 1:
  //   save(3, x/r);
  //   save(2, z/r);
  //   save(1, y/r);
  //   libmesh_fallthrough();
  //
  //   case 0:
  //   save(0, 1);
  // }
}

void
SphericalHarmonics::evaluateStandard()
{

  const Real & x = _standardized_location[1];
  const Real & y = _standardized_location[2];
  const Real & z = _standardized_location[3];
  const Real & r = _standardized_location[0];
  const Real & r2 = r*r;
//   switch (_orders[0])
//   {
//     default:
//     case MAX_DIRECT_CALCULATION_SPHERICALHARMONICS: /*2*/
//     save(24, (x*x*(x*x - 3.0*y*y) - y*y*(3.0*x*x - y*y))/(r*r*r*r));//sin4*cos(4*phi));
//     save(23, (x*x-3*y*y)*x*z/(r*r*r*r));//cos(theta)*sin2*sin(theta)*cos(3*phi));
//     save(22, ((x*x-y*y)*(7*z*z-r*r))/(r*r*r*r));//sin2*(7*cos2-1)*cos(2*phi));
//     save(21, (x*z*(7*z*z - 3*r*r))/(r*r*r*r));//cos(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
//     save(20, (35*z*z*z*z - 30*z*z*r*r + 3*r*r*r*r)/(r*r*r*r));//(35*cos4-30*cos2 + 3));
//     save(19, (y*z*(3*x*x - y*y))/(r*r*r*r));//sin(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
//     save(18, (x*y*(7*z*z - 3*r*r))/(r*r*r*r));//sin(2*phi)*sin2*(7*cos2 -1));
//     save(17, (3*x*x - y*y)*(y*z)/(r*r*r*r));//cos(theta)*sin2*sin(theta)*sin(3*phi));
//     save(16, (x*y*(x*x - y*y))/(r*r));//sin4*sin(4*phi));
//     libmesh_fallthrough();
//
//     case 3:
//     save(15, (x*(x*x - 3.0*y*y))/(r*r*r));//sin2*sin(theta)*cos(3*phi));
//     save(14, (z*(x*x - y*y))/(r*r*r));//cos(theta)*sin2*cos(2*phi));
//     save(13, (x*(4*z*z - x*x - y*y))/(r*r*r));//(5*cos2-1)*sin(theta)*cos(phi));
//     save(12, z*(2*z*z - 3*x*x - 3*y*y)/(r*r*r));//(5*cos2*cos(theta)-3*cos(theta)));
//     save(11, (y*(4*z*z - x*x - y*y))/(r*r*r));//(5*cos2-1)*sin(theta)*sin(phi));
//     save(10, (x*y*z)/(r*r*r));//sin(2*phi)*cos(theta)*sin2);
//     save(9, ((3*x*x - y*y)*y)/(r*r*r));//sin(3*phi)*sin2*sin(theta));
//     libmesh_fallthrough();
//
//     case 2:
//     save(8, (x*x-y*y)/(r*r));//sin2*cos(2.0*phi));
//     save(7, z*x/(r*r));//cos(phi)*sin(theta)*cos(theta));
//     save(6, (-x*x - y*y + 2.0*z*z)/(r*r));
//     save(5, (y * z)/ r2);
//     save(4, (x * y) / r2);
//     libmesh_fallthrough();
//
//     case 1:
//     save(3, x/r);
//     save(2, z/r);
//     save(1, y/r);
//     libmesh_fallthrough();
//
//     case 0:
//     save(0, 1);
// }
  const Real & theta = _standardized_location[1];
  const Real & phi = _standardized_location[2];
  const Real sin2 = sin(theta)*sin(theta);
  const Real sin4 = sin2*sin2;
  const Real cos2 = cos(theta)*cos(theta);
  const Real cos4 = cos2*cos2;
  switch (_orders[0])
  {
    default:
    case MAX_DIRECT_CALCULATION_SPHERICALHARMONICS: /*2*/
    save(24, sin4*cos(4*phi));
    save(23, cos(theta)*sin2*sin(theta)*cos(3*phi));
    save(22, sin2*(7*cos2-1)*cos(2*phi));
    save(21, cos(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
    save(20, 35*cos4-30*cos2 + 3);
    save(19, sin(phi)*sin(theta)*(7*cos2*cos(theta)-3*cos(theta)));
    save(18, sin(2*phi)*sin2*(7*cos2 -1));
    save(17, cos(theta)*sin2*sin(theta)*sin(3*phi));
    save(16, sin4*sin(4*phi));
    libmesh_fallthrough();

    case 3:
    save(15, sin2*sin(theta)*cos(3*phi));
    save(14, cos(theta)*sin2*cos(2*phi));
    save(13, (5*cos2-1)*sin(theta)*cos(phi));
    save(12, (5*cos2*cos(theta)-3*cos(theta)));
    save(11, (5*cos2-1)*sin(theta)*sin(phi));
    save(10, sin(2*phi)*cos(theta)*sin2);
    save(9, sin(3*phi)*sin2*sin(theta));
    libmesh_fallthrough();

    case 2:
    save(8, sin2*cos(2.0*phi));
    save(7, cos(phi)*sin(theta)*cos(theta));
    save(6, 3.0*cos(theta)*cos(theta)-1.0);
    save(5, sin(phi)*sin(theta)*cos(theta));
    save(4, sin(2.0*phi)*sin2);
    libmesh_fallthrough();

    case 1:
    save(3, cos(phi)*sin(theta));
    save(2, cos(theta));
    save(1, sin(phi)*sin(theta));
    libmesh_fallthrough();

    case 0:
    save(0, 1);
  }
}

const std::vector<Real> &
SphericalHarmonics::getStandardizedFunctionLimits() const
{
  static const std::vector<Real> getStandardizedFunctionLimits = {0.0, 1.0, 0, M_PI, 0, 2*M_PI};//0, M_PI, 0, 2*M_PI};
//  static const std::vector<Real> getStandardizedFunctionLimits = {0, 1, -1, 1, -1, 1};
  return getStandardizedFunctionLimits;
}

Real
SphericalHarmonics::getStandardizedFunctionVolume() const
{
  return 1.33333*M_PI; //Volume of a unit sphere
}

std::vector<Real>
SphericalHarmonics::getStandardizedLocation(const std::vector<Real> & location) const
{
  // Get the offset corresponding to the 'x' direction
  const Real offset1 = location[0] - _physical_bounds[0];
  // Get the offset corresponding to the 'y' direction
  const Real offset2 = location[1] - _physical_bounds[1];
  // Get the user-provided radius bound
  const Real & radius = _physical_bounds[3];
  //Get the offset corresponding to the 'z' direction
  const Real & offset3 = location[2] - _physical_bounds[2];

  const Real & norm_x = offset1 / radius;
  const Real & norm_y = offset2 / radius;
  const Real & norm_z = offset3 / radius;
  // Covert to a radis and normalize
  const Real standardizedRadius = std::sqrt(offset1 * offset1 + offset2 * offset2 + offset3*offset3) / radius;
  // Get the angle
  const Real theta = atan2(offset2, offset1);
  //Get phi
  const Real phi = acos(offset3/radius);

  return {standardizedRadius, theta, phi};
//  return {standardizedRadius, norm_x, norm_y, norm_z};
}

bool
SphericalHarmonics::isInPhysicalBounds(const Point & point) const
{
  const std::vector<Real> location = extractLocationFromPoint(point);
  const std::vector<Real> standardized_location = getStandardizedLocation(location);

  if (standardized_location[0] > 1.0 || standardized_location[2] > M_PI || standardized_location[2] < 0)
    return false;
  else
    return true;
}
