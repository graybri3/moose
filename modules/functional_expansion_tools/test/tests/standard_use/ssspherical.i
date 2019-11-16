[Mesh]
  type = FileMesh
  file = 'sphere.e'
[]

[Variables]
  [./m]
  [../]
[]

[Kernels]
  [./diff_m]
    type = Diffusion
    variable = m
  [../]
  [./source_m]
    type = BodyForce
    variable = m
    value = 100
  [../]
[]

[ICs]
  [./start_m]
    type = ConstantIC
    variable = m
    value = 1
  [../]
[]

[BCs]
  [./interface_value]
    type = DirichletBC
    variable = m
    boundary = 1
    value = 100
  [../]
[]

[Functions]
  [./FX_Basis_Value_Main]
    type = FunctionSeries
    series_type = Spherical
    orders = '4'
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

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
  csv = true
[]
