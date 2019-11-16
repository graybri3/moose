//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Cartesian.h"
#include "CylindricalDuo.h"
#include "Legendre.h"
#include "Zernike.h"
#include "Spherical.h"
#include "SphericalHarmonics.h"
#include "VolumetricSpherical.h"

#ifndef SETUP_H
#define SETUP_H

// Set the global tolerances
const double tol = 1e-13;

// Set the name
const std::string name = "UnitTesting";

// Recreate the MooseEnum types used in validParams<FunctionSeries>()
extern MooseEnum single_series_types_1D;
extern MooseEnum single_series_types_2D;
extern MooseEnum single_series_types_3D;
extern MooseEnum expansion_type;
extern MooseEnum generation_type;

#endif // SETUP_H
