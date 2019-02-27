//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Spherical.h"
#include "SphericalHarmonics.h"

Spherical::Spherical(const std::string & who_is_using_me)
  : CompositeSeriesBasisInterface(who_is_using_me)
  {
  }

Spherical::Spherical(const std::vector<MooseEnum> & domains,
                     const std::vector<std::size_t> & orders,
                     const std::vector<MooseEnum> & series_types,
                     const std::string & who_is_using_me,
                     MooseEnum expansion_type,
                     MooseEnum generation_type)
  : CompositeSeriesBasisInterface(orders, series_types, who_is_using_me)
{

if (_series_types[0] == "SphericalHarmonics")
{
  std::vector<MooseEnum> local_domain = {domains[0], domains[1], domains[2]};
  std::vector<std::size_t> local_order = {orders[0]};
  _series.push_back(libmesh_make_unique<SphericalHarmonics>(
      local_domain, local_order, expansion_type, generation_type));
}
else
{
  mooseError("Spherical: No other sphere series implemented except Spherical Harmonics!");
}
  /*
  * Set the _number_of_terms for the composite series by looping over each of the single series.
  * This also initializes _basis_evaluation with zero values and the appropriate length.
  */
  setNumberOfTerms();
}

void
Spherical::setPhysicalBounds(const std::vector<Real> & bounds)
{
  //each single series
  if(bounds.size() != 4)
    mooseError("Spherical: Need x-center, y-center, z-center, radius");

  _series[0]->setPhysicalBounds({bounds[0], bounds[1], bounds[2], bounds[3]});
}
