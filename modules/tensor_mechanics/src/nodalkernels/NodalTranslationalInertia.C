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

#include "NodalTranslationalInertia.h"
#include "MooseVariable.h"
#include "AuxiliarySystem.h"
#include "MooseUtils.h"
#include "DelimitedFileReader.h"

registerMooseObject("TensorMechanicsApp", NodalTranslationalInertia);

template <>
InputParameters
validParams<NodalTranslationalInertia>()
{
  InputParameters params = validParams<TimeNodalKernel>();
  params.addClassDescription("Computes the interial forces and mass proportional damping terms "
                             "corresponding to nodal mass.");
  params.addCoupledVar("velocity", "velocity variable");
  params.addCoupledVar("acceleration", "acceleration variable");
  params.addRangeCheckedParam<Real>(
      "beta", "beta>0.0", "beta parameter for Newmark Time integration");
  params.addRangeCheckedParam<Real>(
      "gamma", "gamma>0.0", "gamma parameter for Newmark Time integration");
  params.addRangeCheckedParam<Real>("eta",
                                    0.0,
                                    "eta>=0.0",
                                    "Constant real number defining the eta parameter for "
                                    "Rayleigh damping.");
  params.addRangeCheckedParam<Real>("alpha",
                                    0.0,
                                    "alpha >= -0.3333 & alpha <= 0.0",
                                    "Alpha parameter for mass dependent numerical damping induced "
                                    "by HHT time integration scheme");
  params.addParam<Real>("mass", "Mass associated with the node");
  params.addParam<FileName>(
      "nodal_mass_file",
      "The file containing the nodal positions and the corresponding nodal masses.");
  return params;
}

NodalTranslationalInertia::NodalTranslationalInertia(const InputParameters & parameters)
  : TimeNodalKernel(parameters),
    _mass(isParamValid("mass") ? getParam<Real>("mass") : 0.0),
    _beta(isParamValid("beta") ? getParam<Real>("beta") : 0.1),
    _gamma(isParamValid("gamma") ? getParam<Real>("gamma") : 0.1),
    _eta(getParam<Real>("eta")),
    _alpha(getParam<Real>("alpha"))
{
  if (isParamValid("beta") && isParamValid("gamma") && isParamValid("velocity") &&
      isParamValid("acceleration"))
  {
    _u_old = &(_var.dofValuesOld());
    _aux_sys = &(_fe_problem.getAuxiliarySystem());

    MooseVariable * vel_variable = getVar("velocity", 0);
    _vel_num = vel_variable->number();
    MooseVariable * accel_variable = getVar("acceleration", 0);
    _accel_num = accel_variable->number();
  }
  else if (!isParamValid("beta") && !isParamValid("gamma") && !isParamValid("velocity") &&
           !isParamValid("acceleration"))
  {
    _vel = &(_var.dofValuesDot());
    _vel_old = &(_var.dofValuesDotOld());
    _accel = &(_var.dofValuesDotDot());
    _du_dot_du = &(_var.duDotDu());
    _du_dotdot_du = &(_var.duDotDotDu());
  }
  else
    mooseError("NodalTranslationalInertia: Either all or none of `beta`, `gamma`, `velocity` and "
               "`acceleration` should be provided as input.");

  if (!isParamValid("nodal_mass_file") && !isParamValid("mass"))
    mooseError(
        "NodalTranslationalInertia: Please provide either mass or nodal_mass_file as input.");
  else if (isParamValid("nodal_mass_file") && isParamValid("mass"))
    mooseError("NodalTranslationalInertia: Please provide either mass or nodal_mass_file as input, "
               "not both.");

  if (isParamValid("nodal_mass_file"))
  {
    MooseUtils::DelimitedFileReader nodal_mass_file(getParam<FileName>("nodal_mass_file"));
    nodal_mass_file.setHeaderFlag(MooseUtils::DelimitedFileReader::HeaderFlag::OFF);
    nodal_mass_file.read();
    std::vector<std::vector<Real>> data = nodal_mass_file.getData();
    if (data.size() != 4)
      mooseError("NodalTranslationalInertia: The number of columns in ",
                 getParam<FileName>("nodal_mass_file"),
                 " should be 4.");

    unsigned int node_found = 0;
    const std::set<BoundaryID> bnd_ids = BoundaryRestrictable::boundaryIDs();
    for (auto & bnd_id : bnd_ids)
    {
      const std::vector<dof_id_type> & bnd_node_set = _mesh.getNodeList(bnd_id);
      for (auto & bnd_node : bnd_node_set)
      {
        const Node & node = _mesh.nodeRef(bnd_node);
        _node_id_to_mass[bnd_node] = 0.0;

        for (unsigned int i = 0; i < data[0].size(); ++i)
        {
          if (MooseUtils::absoluteFuzzyEqual(data[0][i], node(0), 1e-6) &&
              MooseUtils::absoluteFuzzyEqual(data[1][i], node(1), 1e-6) &&
              MooseUtils::absoluteFuzzyEqual(data[2][i], node(2), 1e-6))
          {
            _node_id_to_mass[bnd_node] = data[3][i];
            node_found += 1;
            break;
          }
        }
      }
    }
    if (node_found != data[0].size())
      mooseError("NodalTranslationalInertia: Out of ",
                 data[0].size(),
                 " nodal positions in ",
                 getParam<FileName>("nodal_mass_file"),
                 " only ",
                 node_found,
                 " nodes were found in the boundary.");
  }
}

Real
NodalTranslationalInertia::computeQpResidual()
{
  if (_dt == 0)
    return 0;
  else
  {
    Real mass = 0.0;
    if (isParamValid("mass"))
      mass = _mass;
    else
    {
      if (_node_id_to_mass.find(_current_node->id()) != _node_id_to_mass.end())
        mass = _node_id_to_mass[_current_node->id()];
      else
        mooseError("NodalTranslationalInertia: Unable to find an entry for the current node in the "
                   "_node_id_to_mass map.");
    }

    if (isParamValid("beta"))
    {
      mooseAssert(_beta > 0.0, "NodalTranslationalInertia: Beta parameter should be positive.");
      const NumericVector<Number> & aux_sol_old = _aux_sys->solutionOld();

      const Real vel_old = aux_sol_old(_current_node->dof_number(_aux_sys->number(), _vel_num, 0));
      const Real accel_old =
          aux_sol_old(_current_node->dof_number(_aux_sys->number(), _accel_num, 0));

      const Real accel =
          1. / _beta *
          (((_u[_qp] - (*_u_old)[_qp]) / (_dt * _dt)) - vel_old / _dt - accel_old * (0.5 - _beta));
      const Real vel = vel_old + (_dt * (1 - _gamma)) * accel_old + _gamma * _dt * accel;
      return mass * (accel + vel * _eta * (1 + _alpha) - _alpha * _eta * vel_old);
    }
    else
      return mass * ((*_accel)[_qp] + (*_vel)[_qp] * _eta * (1.0 + _alpha) -
                     _alpha * _eta * (*_vel_old)[_qp]);
  }
}

Real
NodalTranslationalInertia::computeQpJacobian()
{
  mooseAssert(_beta > 0.0, "NodalTranslationalInertia: Beta parameter should be positive.");

  if (_dt == 0)
    return 0;
  else
  {
    Real mass = 0.0;
    if (isParamValid("mass"))
      mass = _mass;
    else
    {
      if (_node_id_to_mass.find(_current_node->id()) != _node_id_to_mass.end())
        mass = _node_id_to_mass[_current_node->id()];
      else
        mooseError("NodalTranslationalInertia: Unable to find an entry for the current node in the "
                   "_node_id_to_mass map.");
    }

    if (isParamValid("beta"))
      return mass / (_beta * _dt * _dt) + _eta * (1 + _alpha) * mass * _gamma / _beta / _dt;
    else
      return mass * (*_du_dotdot_du)[_qp] + _eta * (1.0 + _alpha) * mass * (*_du_dot_du)[_qp];
  }
}
