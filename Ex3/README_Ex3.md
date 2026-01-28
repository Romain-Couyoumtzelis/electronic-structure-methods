# Exercise 3: Basis Sets, Dissociation Energy, and Geometry Optimization

This exercise explores the effect of basis set size on H₂ energies, compares RHF vs UHF methods through dissociation curves, and introduces geometry optimization.

## Learning Goals

- Understand the influence of basis sets on energy calculations
- Learn the geometry optimization procedure
- Understand the difference between RHF and UHF methods

## Notebooks

### IESM_Ex3.ipynb
Main introduction notebook providing an overview of the exercise, learning goals, and references to slides and resources.

### Ex3_H2molecule.ipynb
Effect of basis set size on molecular energies:
- Computing H₂ equilibrium energy with different basis sets (6-31G, 6-311G, aug-cc-pVTZ)
- Understanding SCF convergence behavior
- Comparing calculated energies to the exact value (E_exact = -1.174474 a.u.)
- Plotting SCF convergence across basis sets
- Understanding DIIS (Direct Inversion in the Iterative Subspace)

### Ex3_Dissociation.ipynb
Recording dissociation curves for H₂:
- **RHF (Restricted Hartree-Fock)**: Both electrons share the same spatial orbital
- **UHF (Unrestricted Hartree-Fock)**: α and β electrons can have different spatial orbitals
- Calculating energies at various internuclear distances (0.5-6.0 Å)
- Computing and plotting interaction energies
- Understanding why RHF fails at dissociation (ionic character)
- Physical interpretation of the RHF vs UHF difference

### Ex3_GeometryOpt.ipynb
Geometry optimization of H₂O:
- Understanding the potential energy surface (PES) concept
- Starting from suboptimal geometries (wrong bond lengths and angles)
- Using `psi4.optimize()` instead of `psi4.energy()`
- Analyzing convergence criteria (Delta E, Max Force, Max Disp)
- Visualizing the optimization trajectory on a 3D PES
- Understanding why multiple SCF cycles occur per optimization step

### optimizations1.ipynb & optimizations2.ipynb
Additional geometry optimization examples:
- Starting from different initial geometries (angles: 120° and 165°)
- Comparing optimization paths on the PES
- Visualizing how different starting points converge to the same minimum

## Key Concepts

- **Potential Energy Surface (PES)**: Energy as a function of nuclear coordinates
- **Geometry Optimization**: Finding the minimum energy geometry by following gradients
- **RHF vs UHF**: RHF enforces paired electrons in same spatial orbital; UHF allows spin polarization
- **Dissociation limit**: At infinite separation, H₂ should give 2×E(H atom)
- **Convergence criteria**: Energy change, forces, and displacements must fall below thresholds

## References

- Jensen, F. (2017). *Introduction to Computational Chemistry*. John Wiley & Sons.
  - Chapter 5: Basis Sets
- Sherrill Group Resources:
  - [Hartree-Fock Introduction (slides)](http://vergil.chemistry.gatech.edu/courses/chem6485/pdf/Hartree-Fock-Intro.pdf)
  - [Geometry Optimization (slides)](http://vergil.chemistry.gatech.edu/courses/chem6485/pdf/geom-opt.pdf)
