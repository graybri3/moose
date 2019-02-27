/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef NODALTRANSLATIONALINERTIA_H
#define NODALTRANSLATIONALINERTIA_H

#include "TimeNodalKernel.h"

// Forward Declarations
class NodalTranslationalInertia;

template <>
InputParameters validParams<NodalTranslationalInertia>();

/**
 * Calculates the inertial force and mass proportional damping for a nodal mass
 */
class NodalTranslationalInertia : public TimeNodalKernel
{
public:
  NodalTranslationalInertia(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  /// Mass associated with the node
  const Real _mass;

  /// Old value of displacement
  const VariableValue * _u_old;

  /// Newmark time integration parameter
  const Real _beta;

  /// Newmark time integration parameter
  const Real _gamma;

  /// Mass proportional Rayliegh damping
  const Real & _eta;

  /// HHT time integration parameter
  const Real & _alpha;

  /// Auxiliary system object
  AuxiliarySystem * _aux_sys;

  /// Variable number corresponding to the velocity aux variable
  unsigned int _vel_num;

  /// Variable number corresponding to the acceleration aux variable
  unsigned int _accel_num;

  /// Map between boundary nodes and nodal mass
  std::map<dof_id_type, Real> _node_id_to_mass;

  /// Velocity variable value
  const MooseArray<Number> * _vel;

  /// Old velocity variable value
  const MooseArray<Number> * _vel_old;

  /// Acceleration variable value
  const MooseArray<Number> * _accel;

  /// du_dot_du variable value
  const MooseArray<Number> * _du_dot_du;

  /// du_dotdot_du variable value
  const MooseArray<Number> * _du_dotdot_du;
};

#endif /* NODALTRANSLATIONALINERTIA_H */
