//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADMATHEATSOURCE_H
#define ADMATHEATSOURCE_H

#include "ADKernel.h"

template <ComputeStage>
class ADMatHeatSource;

declareADValidParams(ADMatHeatSource);

template <ComputeStage compute_stage>
class ADMatHeatSource : public ADKernel<compute_stage>
{
public:
  ADMatHeatSource(const InputParameters & parameters);

protected:
  virtual ADResidual computeQpResidual();

  const Real _scalar;
  const ADMaterialProperty(Real) & _material_property;

  usingKernelMembers;
};

#endif // ADMATHEATSOURCE_H
