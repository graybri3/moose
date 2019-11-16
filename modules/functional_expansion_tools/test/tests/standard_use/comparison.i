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
