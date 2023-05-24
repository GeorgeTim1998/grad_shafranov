## Project name: Grad-Shafranov
- Utilizes fenics python package for solving partial differential equation using finite element method and mshr pachage for creating arbitrary calculations domains and finite elements meshes 
- The main purpose of this project is to solve Grad-Shafranov partial differential equation in various tokamak installations (MEPhIST, ITER, FRC, Spheromak etc)
- Allows modelling of: 
```
1. Soloviev equilibrium
2. Plasma equilibrium inside plasma vessel and outside controlling poloidal coils
3. Plasma sudden shifts from the equilibrium inside a conducting vessel
```

## How to set up project
- Clone the repository and navigate to created folder
- `conda env create -f environment.yml`

## How to launch
- `conda activate fenicsproject`
