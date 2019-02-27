//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef LINEARINTERPOLATION_H
#define LINEARINTERPOLATION_H

#include <vector>
#include <string>

#include "Moose.h"

/**
 * This class interpolates values given a set of data pairs and an abscissa.
 */
class LinearInterpolation
{
public:
  /* Constructor, Takes two vectors of points for which to apply the fit.  One should be of the
   * independent variable while the other should be of the dependent variable.  These values should
   * correspond to one and other in the same position.
   */
  LinearInterpolation(const std::vector<Real> & X, const std::vector<Real> & Y);
  LinearInterpolation() : _x(std::vector<Real>()), _y(std::vector<Real>()) {}

  virtual ~LinearInterpolation() = default;

  /**
   * Set the x and y values.
   */
  void setData(const std::vector<Real> & X, const std::vector<Real> & Y)
  {
    _x = X;
    _y = Y;
    errorCheck();
  }

  void errorCheck();

  /**
   * This function will take an independent variable input and will return the dependent variable
   * based on the generated fit
   */
  Real sample(Real x) const;

  /**
   * This function will take an independent variable input and will return the derivative of the
   * dependent variable
   * with respect to the independent variable based on the generated fit
   */
  Real sampleDerivative(Real x) const;

  /**
   * This function returns the size of the array holding the points, i.e. the number of sample
   * points
   */
  unsigned int getSampleSize();

  /**
   * This function returns the integral of the function
   */
  Real integrate();

  Real domain(int i) const;
  Real range(int i) const;

private:
  std::vector<Real> _x;
  std::vector<Real> _y;

  static int _file_number;
};

#endif // LINEARINTERPOLATION_H
