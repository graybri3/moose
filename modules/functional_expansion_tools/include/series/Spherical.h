//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPHERICAL_H
#define SPHERICAL_H

#include "CompositeSeriesBasisInterface.h"

/**
  * This class constructsa functional expansion using a separate series for each
  * Spherical dimension. 3D
**/
class Spherical final : public CompositeSeriesBasisInterface
{
public:
  Spherical(const std::string & who_is_using_me);
  Spherical(const std::vector<MooseEnum> & domain,
            const std::vector<std::size_t> & order,
            const std::vector<MooseEnum> & series_types,
            const std::string & who_is_using_me,
            MooseEnum expansion_type,
            MooseEnum generation_type);

  // Overrides from FunctionalBasisInterface
  virtual void setPhysicalBounds(const std::vector<Real> & bounds) final;
};

#endif // SPHERICAL_H
