# 1phase, heat advecting with a moving fluid
# Using the PorousFlowFullySaturated Action with KT stabilization
# This should produce an identical result to heat_advection_1D_KT
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 50
  xmin = 0
  xmax = 1
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [./temp]
    initial_condition = 200
  [../]
  [./pp]
  [../]
[]

[ICs]
  [./pp]
    type = FunctionIC
    variable = pp
    function = '1-x'
  [../]
[]

[BCs]
  [./pp0]
    type = PresetBC
    variable = pp
    boundary = left
    value = 1
  [../]
  [./pp1]
    type = PresetBC
    variable = pp
    boundary = right
    value = 0
  [../]
  [./spit_heat]
    type = PresetBC
    variable = temp
    boundary = left
    value = 300
  [../]
  [./suck_heat]
    type = PresetBC
    variable = temp
    boundary = right
    value = 200
  [../]
[]

[PorousFlowFullySaturated]
  porepressure = pp
  temperature = temp
  coupling_type = ThermoHydro
  fp = simple_fluid
  add_darcy_aux = false
  stabilization = KT
  flux_limiter_type = superbee
[]

[Modules]
  [./FluidProperties]
    [./simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 100
      density0 = 1000
      viscosity = 4.4
      thermal_expansion = 0
      cv = 2
    [../]
  [../]
[]

[Materials]
  [./porosity]
    type = PorousFlowPorosityConst
    porosity = 0.2
  [../]
  [./zero_thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '0 0 0  0 0 0  0 0 0'
  [../]
  [./rock_heat]
    type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 1.0
    density = 125
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1.1 0 0 0 2 0 0 0 3'
  [../]
[]

[Preconditioning]
  [./andy]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres bjacobi 1E-15 1E-10 10000'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 0.01
  end_time = 0.6
[]

[Outputs]
  file_base = heat_advection_1d_KT_out
  exodus = true
  interval = 10
[]
