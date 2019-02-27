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
    value = 1.00001
  [../]
[]

[BCs]
[./bottom]
    type = DirichletBC
    variable = m
    boundary = 1
    value = 1
  [../]
[]

[Functions]
  [./FX_Basis_Value_Main]
    type = FunctionSeries
    series_type = Spherical
    orders = '4'
    physical_bounds = '0.0 0.0 0.0 10.0'
    sphere = SphericalHarmonics
    generation_type = 'sqrt_mu'
    expansion_type = 'sqrt_mu'
    print_when_set = true # Print coefficients when a MultiAppFXTransfer is executed
  [../]
[]

[UserObjects]
  [./FX_Value_UserObject_Main]
    type = FXVolumeUserObject
    function = FX_Basis_Value_Main
    variable = m
    boundary = 1
    print_state = true # Print after the FX coefficients are computer
    print_when_set = true # Print coefficients when a MultiAppFXTransfer is executed
  [../]
[]

[Executioner]
  type = Steady
  num_steps = 4
  dt = 1
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
