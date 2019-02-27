[Mesh]
  type = FileMesh
  file = 'sphere.e'
[]

[Variables]
  [./m]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./s_in]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./diff_m]
    type = Diffusion
    variable = m
  [../]
[]

[ICs]
  [./start_m]
    type = ConstantIC
    variable = m
    value = 1000
  [../]
[]

[BCs]
  [./interface_value]
    type = FXFluxBC
    variable = m
    boundary = 1
    function = FX_Basis_Value_Main
  [../]
[]

[Functions]
  [./FX_Basis_Value_Main]
    type = FunctionSeries
    series_type = Spherical
    orders = '3'
    physical_bounds = '0.0 0.0 0.0 1.0'
    sphere = SphericalHarmonics
    generation_type = 'standard'
    expansion_type = 'standard'
    print_when_set = true # Print coefficients when a MultiAppFXTransfer is executed
  [../]
[]

[UserObjects]
  [./FX_Value_UserObject_Main]
    type = FXVolumeUserObject
    function = FX_Basis_Value_Main
    variable = m
    print_state = true # Print after the FX coefficients are computer
    print_when_set = true # Print coefficients when a MultiAppFXTransfer is executed
  [../]
[]

[Postprocessors]
  [./average_value]
    type = ElementAverageValue
    variable = m
  [../]
  [./peak_value]
    type = ElementExtremeValue
    value_type = max
    variable = m
  [../]
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
# petsc_options_iname = '-pc_type -pc_hypre_type'
#  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
