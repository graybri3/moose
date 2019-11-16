//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseUtils.h"
#include "Spherical.h"
#include <functional>

#define MAX_DIRECT_CALCULATION_VOLUMETRICSPHERICAL 4

VolumetricSpherical::VolumetricSpherical() : SingleSeriesBasisInterface() {}

VolumetricSpherical::VolumetricSpherical(const std::vector<MooseEnum> & domain,
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
VolumetricSpherical::calculatedNumberOfTermsBasedOnOrder(const std::vector<std::size_t> & order) const
{
  return ((order[0] + 1) * (order[0] + 2)) / 2;
}

void
VolumetricSpherical::checkPhysicalBounds(const std::vector<Real> & bounds) const
{
  if (bounds.size() != 4)
    mooseError("VolumetricSpherical: Invalid number of bounds specified for a single series!");
}
void
VolumetricSpherical::evaluateOrthonormal()
{
  const Real & rho = _standardized_location[0];
//  const Real & theta = _standardized_location[1];
//  const Real & phi = _standardized_location[2];
  const Real rho2 = rho * rho;
  const Real & x = _standardized_location[1];
  const Real & y = _standardized_location[2];
  const Real & z = _standardized_location[3];
//  const Real sin2 = sin(theta) * sin(theta);
//  const Real cos2 = cos(theta) * cos(theta);

  // switch(_orders[0])
  // {
  //   default:
  //   case MAX_DIRECT_CALCULATION_VOLUMETRICSPHERICAL:
  //     save(14, std::sqrt(3465 / (256 * M_PI)) * rho2 * rho2 * sin2 * sin2 * (8*cos(phi)*cos(phi)*cos(phi)*cos(phi) - 8*cos(phi)*cos(phi) + 1));
  //     save(13, std::sqrt(495 / (64 * M_PI)) * (4.5 * rho2 * rho2 - 3.5 * rho2) * sin2 * (2 * cos(phi) - 1) * (7 * cos2 - 1));
  //     save(12, std::sqrt(99 / (256 * M_PI)) * (7.875 * rho2 * rho2 - 8.75 * rho2 + 1.875) * (35 * cos2 * cos2 - 30 * cos2 + 3));
  //     save(11, std::sqrt(495 / (64 * M_PI)) * (4.5 * rho2 * rho2 - 3.5 * rho2) * sin2 * sin(2 * phi) * (7 * cos2 - 1));
  //     save(10, std::sqrt(3465 / (64 * M_PI)) * rho2 * rho2 * sin2 * sin2 * sin(2 * phi) * (1 - 2 * sin(phi) * sin(phi)));
  //     libmesh_fallthrough();
  //   case 3:
  //     save(9, std::sqrt(315 / (32 * M_PI)) * rho2 * rho * sin2 * sin(theta) * cos(phi) * (4 * cos(phi) * cos(phi) - 3));
  //     save(8, std::sqrt(189 / (32 * M_PI)) * (3.5 * rho2 * rho - 2.5 * rho) * sin(theta) * cos(phi) * (5*cos2 - 1));
  //     save(7, std::sqrt(189 / (32 * M_PI)) * (3.5 * rho2 * rho - 2.5 * rho) * sin(theta) * sin(phi) * (5*cos2 - 1));
  //     save(6, std::sqrt(315 / (32 * M_PI)) * rho2 * rho * sin2*sin(theta) * sin(phi) * (3 - 4 * sin(phi) * sin(phi)));
  //     libmesh_fallthrough();
  //   case 2:
  //     save(5, std::sqrt(105 / (16 * M_PI)) * rho2 * sin2 * cos(2*phi));
  //     save(4, std::sqrt(35 / (16 * M_PI))*(2.5 * rho2 - 1.5)*(3 * cos2 - 1));
  //     save(3, std::sqrt(105 / (16 * M_PI))*rho2*sin2*sin(2 * phi));
  //     libmesh_fallthrough();
  //   case 1:
  //     save(2, std::sqrt(15 / (4 * M_PI))*rho*sin(theta)*sin(phi));
  //     save(1, std::sqrt(15 / (4 * M_PI))*rho*sin(theta)*sin(phi));
  //     libmesh_fallthrough();
  //   case 0:
  //     save(0, std::sqrt(3 / (4 * M_PI)));
  // }
  switch(_orders[0])
  {
    default:
    case MAX_DIRECT_CALCULATION_VOLUMETRICSPHERICAL:
      save(14, std::sqrt(3465 / (256 * M_PI)) * x * x * (x*x - 3*y*y) - y*y * (3*x*x - y*y));
      save(13, std::sqrt(495 / (64 * M_PI)) * (4.5 * rho2 * rho2 - 3.5 * rho2) * ((x*x - y*y)*(7*z*z - rho*rho))/(rho*rho*rho*rho));
      save(12, std::sqrt(99 / (256 * M_PI)) * (35 * z *z*z*z - 30 *z*z*rho2 + 3*rho2*rho2)/(rho2*rho2));
      save(11, std::sqrt(495 / (64 * M_PI)) * (4.5 * rho2 * rho2 - 3.5 * rho2) * x * y *(7*z*z - rho2)/(rho2*rho2));
      save(10, std::sqrt(3465 / (64 * M_PI)) * x * y * (x*x - y*y));
      libmesh_fallthrough();
    case 3:
      save(9, std::sqrt(315 / (32 * M_PI)) * (x*x - 3 * y * y) * x);
      save(8, std::sqrt(189 / (32 * M_PI)) * (3.5 * rho2 * rho - 2.5 * rho) * x * (4*z*z - x*x -y*y) / (rho2 * rho));
      save(7, std::sqrt(189 / (32 * M_PI)) * (3.5 * rho2 * rho - 2.5 * rho) * y * (4*z*z - x*x -y*y) / (rho2 * rho));
      save(6, std::sqrt(315 / (32 * M_PI)) * (3 * x * x - y * y) * y);
      libmesh_fallthrough();
    case 2:
      save(5, std::sqrt(105 / (16 * M_PI)) * (x*x - y*y));
      save(4, std::sqrt(35 / (16 * M_PI))*(2.5 * rho2 - 1.5) * (-x * x - y * y + 2*z*z)/(rho2));
      save(3, std::sqrt(105 / (16 * M_PI)) * 2 * x * y);
      libmesh_fallthrough();
    case 1:
      save(2, std::sqrt(15 / (4 * M_PI))*x);
      save(1, std::sqrt(15 / (4 * M_PI))*y);
      libmesh_fallthrough();
    case 0:
      save(0, std::sqrt(3 / (4 * M_PI)));
  }
}
void
VolumetricSpherical::evaluateSqrtMu()
{
  const Real & rho = _standardized_location[0];
//  const Real & theta = _standardized_location[1];
//  const Real & phi = _standardized_location[2];
  const Real rho2 = rho * rho;
//  const Real sin2 = sin(theta) * sin(theta);
//  const Real cos2 = cos(theta) * cos(theta);
  const Real & x = _standardized_location[1];
  const Real & y = _standardized_location[2];
  const Real & z = _standardized_location[3];
  switch(_orders[0])
  {
    default:
    case MAX_DIRECT_CALCULATION_VOLUMETRICSPHERICAL:
      save(14, std::sqrt(3465 / (256 * M_PI)) * x * x * (x*x - 3*y*y) - y*y * (3*x*x - y*y));
      save(13, std::sqrt(495 / (64 * M_PI)) * (4.5 * rho2 * rho2 - 3.5 * rho2) * ((x*x - y*y)*(7*z*z - rho*rho))/(rho*rho*rho*rho));
      save(12, std::sqrt(99 / (256 * M_PI)) *((7.875 * rho2 - 8.75) * rho2 + 1.875)* (35 * z *z*z*z - 30 *z*z*rho2 + 3*rho2*rho2)/(rho2*rho2));
      save(11, std::sqrt(495 / (64 * M_PI)) * (4.5 * rho2 * rho2 - 3.5 * rho2) * x * y * (7*z*z - rho2)/(rho2*rho2));
      save(10, std::sqrt(3465 / (64 * M_PI)) * x * y * (x*x - y*y));
      libmesh_fallthrough();
    case 3:
      save(9, std::sqrt(315 / (32 * M_PI)) * (x*x - 3 * y * y) * x);
      save(8, std::sqrt(189 / (32 * M_PI)) * (3.5 * rho2 * rho - 2.5 * rho) * x * (4*z*z - x*x -y*y) / (rho2 * rho));
      save(7, std::sqrt(189 / (32 * M_PI)) * (3.5 * rho2 * rho - 2.5 * rho) * y * (4*z*z - x*x -y*y) / (rho2 * rho));
      save(6, std::sqrt(315 / (32 * M_PI)) * (3 * x * x - y * y) * y);
      libmesh_fallthrough();
    case 2:
      save(5, std::sqrt(105 / (16 * M_PI)) * (x*x - y*y));
      save(4, std::sqrt(35 / (16 * M_PI))*(2.5 * rho2 - 1.5) * (-x * x - y * y + 2*z*z)/(rho2));
      save(3, std::sqrt(105 / (16 * M_PI)) * 2 * x * y);
      libmesh_fallthrough();
    case 1:
      save(2, std::sqrt(15 / (4 * M_PI))*x);
      save(1, std::sqrt(15 / (4 * M_PI))*y);
      libmesh_fallthrough();
    case 0:
      save(0, std::sqrt(3 / (4 * M_PI)));
  }
  // switch(_orders[0])
  // {
  //   default:
  //   case MAX_DIRECT_CALCULATION_VOLUMETRICSPHERICAL:
  //     save(14, std::sqrt(3465 / (256 * M_PI)) * x * x * (x*x - 3*y*y) - y*y * (3*x*x - y*y));
  //     save(13, std::sqrt(495 / (64 * M_PI)) * (4.5 * rho2 * rho2 - 3.5 * rho2) * ((x*x - y*y)*(7*z*z - rho*rho))/(rho*rho*rho*rho));
  //     save(12, std::sqrt(99 / (256 * M_PI)) * (35 * z *z*z*z - 30 *z*z*rho2 + 3*rho2*rho2)/(rho2*rho2));
  //     save(11, std::sqrt(495 / (64 * M_PI)) * (4.5 * rho2 * rho2 - 3.5 * rho2) * ((x*x - y*y) - (7*z*z - rho2))/(rho2*rho2)));
  //     save(10, std::sqrt(3465 / (64 * M_PI)) * x * y * (x*x - y*y));
  //     libmesh_fallthrough();
  //   case 3:
  //     save(9, std::sqrt(315 / (32 * M_PI)) * (x*x - 3 * y * y) * y);
  //     save(8, std::sqrt(189 / (32 * M_PI)) * (3.5 * rho2 * rho - 2.5 * rho) * x * (4*z*z - x*x -y*y) / (rho2 * rho));
  //     save(7, std::sqrt(189 / (32 * M_PI)) * (3.5 * rho2 * rho - 2.5 * rho) * y * (4*z*z - x*x -y*y) / (rho2 * rho));
  //     save(6, std::sqrt(315 / (32 * M_PI)) * (3 * x * x - y * y) * y);
  //     libmesh_fallthrough();
  //   case 2:
  //     save(5, std::sqrt(105 / (16 * M_PI)) * (x*x - y*y));
  //     save(4, std::sqrt(35 / (16 * M_PI))*(2.5 * rho2 - 1.5) * (-x * x + y * y + 2*z*z)/(rho2));
  //     save(3, std::sqrt(105 / (16 * M_PI)) * 2 * x * y);
  //     libmesh_fallthrough();
  //   case 1:
  //     save(2, std::sqrt(15 / (4 * M_PI))*x);
  //     save(1, std::sqrt(15 / (4 * M_PI))*y);
  //     libmesh_fallthrough();
  //   case 0:
  //     save(0, std::sqrt(3 / (4 * M_PI)));
  // }
}
void
VolumetricSpherical::evaluateStandard()
{
  const Real & rho = _standardized_location[0];
  const Real rho2 = rho * rho;
  const Real & x = _standardized_location[1];
  const Real & y = _standardized_location[2];
  const Real & z = _standardized_location[3];

  // const Real & theta = _standardized_location[1];
  // const Real & phi = _standardized_location[2];
  // const Real sin2 = sin(theta) * sin(theta);
  // const Real cos2 = cos(theta) * cos(theta);

  // switch(_orders[0])
  // {
  //   default:
  //   case MAX_DIRECT_CALCULATION_VOLUMETRICSPHERICAL:
  //     save(14, rho2 * rho2 * sin2 * sin2 * (8*cos(phi)*cos(phi)*cos(phi)*cos(phi) - 8*cos(phi)*cos(phi) + 1));
  //     save(13, (4.5 * rho2 * rho2 - 3.5 * rho2) * sin2 * (2 * cos(phi) - 1) * (7 * cos2 - 1));
  //     save(12, (7.875 * rho2 * rho2 - 8.75 * rho2 + 1.875) * (35 * cos2 * cos2 - 30 * cos2 + 3));
  //     save(11, (4.5 * rho2 * rho2 - 3.5 * rho2) * sin2 * sin(2 * phi) * (7 * cos2 - 1));
  //     save(10, rho2 * rho2 * sin2 * sin2 * sin(2 * phi) * (1 - 2 * sin(phi) * sin(phi)));
  //     libmesh_fallthrough();
  //   case 3:
  //     save(9, rho2 * rho * sin2 * sin(theta) * cos(phi) * (4 * cos(phi) * cos(phi) - 3));
  //     save(8, (3.5 * rho2 * rho - 2.5 * rho) * sin(theta) * cos(phi) * (5*cos2 - 1));
  //     save(7, (3.5 * rho2 * rho - 2.5 * rho) * sin(theta) * sin(phi) * (5*cos2 - 1));
  //     save(6, rho2 * rho * sin2*sin(theta) * sin(phi) * (3 - 4 * sin(phi) * sin(phi)));
  //     libmesh_fallthrough();
  //   case 2:
  //     save(5, rho2 * sin2 * cos(2*phi));
  //     save(4, (2.5 * rho2 - 1.5)*(3 * cos2 - 1));
  //     save(3, rho2 * sin2 * sin(2 * phi));
  //     libmesh_fallthrough();
  //   case 1:
  //     save(2, rho*sin(theta)*sin(phi));
  //     save(1, rho*sin(theta)*sin(phi));
  //     libmesh_fallthrough();
  //   case 0:
  //     save(0, 1);
  // }
  switch(_orders[0])
  {
    default:
    case MAX_DIRECT_CALCULATION_VOLUMETRICSPHERICAL:
      save(14, x * x * (x*x - 3*y*y) - y*y * (3*x*x - y*y));
      save(13, (4.5 * rho2 * rho2 - 3.5 * rho2) * ((x*x - y*y)*(7*z*z - rho*rho))/(rho*rho*rho*rho));
      save(12, ((7.875 * rho2 - 8.75) * rho2 + 1.875) * (35 * z *z*z*z - 30 *z*z*rho2 + 3*rho2*rho2)/(rho2*rho2));
      save(11, (4.5 * rho2 * rho2 - 3.5 * rho2) * x * y * (7*z*z - rho2)/(rho2*rho2));
      save(10, x * y * (x*x - y*y));
      libmesh_fallthrough();
    case 3:
      save(9, (x*x - 3 * y * y) * x);
      save(8, (3.5 * rho2 * rho - 2.5 * rho) * x * (4*z*z - x*x -y*y) / (rho2 * rho));
      save(7, (3.5 * rho2 * rho - 2.5 * rho) * y * (4*z*z - x*x -y*y) / (rho2 * rho));
      save(6, (3 * x * x - y * y) * y);
      libmesh_fallthrough();
    case 2:
      save(5, (x*x - y*y));
      save(4, (2.5 * rho2 - 1.5) * (-x * x - y * y + 2*z*z)/(rho2));
      save(3, x * y);
      libmesh_fallthrough();
    case 1:
      save(2, x);
      save(1, y);
      libmesh_fallthrough();
    case 0:
      save(0, 1);
  }
}

const std::vector<Real> &
VolumetricSpherical::getStandardizedFunctionLimits() const
{
  //static const std::vector<Real> standardizedFunctionLimits = {0, 1, -M_PI, M_PI, 0, M_PI};
  static const std::vector<Real> standardizedFunctionLimits = {0, 1, -1, 1, -1, 1, -1, 1};
  return standardizedFunctionLimits;
}

Real
VolumetricSpherical::getStandardizedFunctionVolume() const
{
  return 1.3333 * M_PI;
}

std::vector<Real>
VolumetricSpherical::getStandardizedLocation(const std::vector<Real> & location) const
{
  const Real offset1 = location[0] - _physical_bounds[0];
  const Real offset2 = location[1] - _physical_bounds[1];
  const Real offset3 = location[2] - _physical_bounds[2];

  const Real & radius = _physical_bounds[3];

  const Real norm_x = offset1 / radius;
  const Real norm_y = offset2 / radius;
  const Real norm_z = offset3 / radius;

  const Real standardizedRadius = sqrt(offset1 * offset1 + offset2 * offset2 + offset3 * offset3) / radius;

  const Real phi = atan2(offset2, offset1);

  const Real theta = acos(offset3/radius);
  const Real norm_r = sqrt(norm_x*norm_x + norm_y*norm_y + norm_z*norm_z);

//  return{standardizedRadius, theta, phi};
/*
* Radius is from 0 to 1 and legendre space needs converted to -1 to 1
*/
  return{standardizedRadius, norm_x, norm_y, norm_z};
}
bool
VolumetricSpherical::isInPhysicalBounds(const Point & point) const
{

  const std::vector<Real> location = extractLocationFromPoint(point);
  const std::vector<Real> standardized_location = getStandardizedLocation(location);

  /*
   * The radius (standardized_location[0]) is always positive, so only check
   * against the maximum radius (1). The theta components should always be in
   * bounds.
   */
  if (standardized_location[0] > 1.0) //|| standardized_location[1] > 1.0 || standardized_location[2] > 1.0|| standardized_location[3] > 1.0)
    return false;
  else
    return true;
}
std::size_t
VolumetricSpherical::simpleDoubleToSingle(std::size_t n, long m) const
{
  return (n * (n + 2) + m) / 2;
}
