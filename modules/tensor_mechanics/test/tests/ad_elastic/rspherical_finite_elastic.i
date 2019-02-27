[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 5
[]

[Problem]
  coord_type = RSPHERICAL
[]

[GlobalParams]
  displacements = 'disp_r'
[]

[Variables]
  # scale with one over Young's modulus
  [./disp_r]
    scaling = 1e-10
  [../]
[]

[Kernels]
  [./stress_r]
    type = ADStressDivergenceRSphericalTensors
    component = 0
    variable = disp_r
    use_displaced_mesh = true
  [../]
[]

[BCs]
  [./center]
    type = PresetBC
    variable = disp_r
    boundary = left
    value = 0
  [../]
  [./rdisp]
    type = PresetBC
    variable = disp_r
    boundary = right
    value = 0.1
  [../]
[]

[Materials]
  [./elasticity]
    type = ComputeIsotropicElasticityTensor
    poissons_ratio = 0.3
    youngs_modulus = 1e10
  [../]
[]

[Materials]
  [./strain]
    type = ADComputeRSphericalFiniteStrain
  [../]
  [./stress]
    type = ADComputeFiniteStrainElasticStress
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.05
  solve_type = 'NEWTON'

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomeramg

  dtmin = 0.05
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
