[Mesh]
  type = FileMesh
  file = 'sphere_sector_3d.e'
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
    boundary = '1 2 3'
  [../]
[]

[BCs]
  [./interface_value]
    type = FXValueBC
    variable = m
    boundary = '1 2 3'
    function = FX_Basis_Value_Main
  [../]
[]

[Functions]
  [./FX_Basis_Value_Main]
    type = FunctionSeries
    series_type = Spherical
    orders = '4'
    physical_bounds = '0.0 0.0 0.0 10.0'
    sphere = VolumetricSpherical
    generation_type = 'orthonormal'
    expansion_type = 'orthonormal'
    print_when_set = true # Print coefficients when a MultiAppFXTransfer is executed
  [../]
[]

[UserObjects]
  [./FX_Value_UserObject_Main]
    type = FXVolumeUserObject
    function = FX_Basis_Value_Main
    variable = m
    boundary = '5'
    print_state = true # Print after the FX coefficients are computer
    print_when_set = true # Print coefficients when a MultiAppFXTransfer is executed
  [../]
[]

[Executioner]
  type = Steady
 # num_steps = 4
 # dt = 1
 # solve_type = PJFNK
 # petsc_options_iname = '-pc_type -pc_hypre_type'
 # petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
