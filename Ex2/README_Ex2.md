# Exercise 2: First Steps in Psi4 - Basis Sets and Hartree-Fock

This exercise introduces the Psi4 quantum chemistry software and explores the fundamentals of basis sets through practical calculations.

## Learning Goals

- Understand what basis functions are
- Get familiar with basis function notation
- Learn the basics of Psi4
- Get familiar with molecular geometries (Z-matrix and Cartesian coordinates)

## Notebooks

### IESM_Ex2.ipynb
Main introduction notebook providing an overview of the exercise, links to slides and report templates.

### Ex2_theory.ipynb
Theoretical background on basis sets:
- **LCAO**: Linear Combination of Atomic Orbitals approach
- **Slater-Type Orbitals (STOs)**: Mathematical form and properties
- **Gaussian-Type Orbitals (GTOs)**: Why Gaussians are computationally advantageous
- **Contracted GTOs (CGTOs)**: STO-3G minimal basis set
- **Pople-Type Split-Valence Basis Sets**: Understanding notation (3-21G, 6-31G, etc.)
- Polarization and diffuse functions

### Ex2_hydrogen_atom.ipynb
First hands-on tutorial with Psi4:
- Setting up the Psi4 environment in Python
- Defining molecular geometry (single atoms)
- Running single-point energy calculations
- Comparing basis sets (STO-3G, 6-31G, 6-311G) for H and He atoms
- Comparing computed energies to analytical values
- Understanding output log files

### Ex2_water-molecule.ipynb
Working with polyatomic molecules:
- Building molecular geometries using Z-matrix format
- Calculating H₂O energies with various basis sets (6-31G, 6-311G, 6-31+G, 6-31G**, 6-31++G**)
- Understanding the effect of polarization and diffuse functions
- Comparing bent vs. linear water molecule geometries
- Comparing H₂O and BeH₂ molecular shapes (Walsh diagram concepts)
- Analyzing Psi4 output files with `grep` and `tail`

## Key Concepts

- **Single-point calculation**: Computing energy at a fixed geometry
- **Spin multiplicity**: 2S+1, where S is total spin
- **SCF convergence**: Iterative optimization of wavefunction coefficients
- **Basis set hierarchy**: Minimal → Split-valence → Polarized → Diffuse

## References

- Jensen, F. (2017). *Introduction to Computational Chemistry*. John Wiley & Sons.
  - Chapter 5: Basis Sets
- Smith DG, et al. (2020). PSI4 1.4: Open-source software for high-throughput quantum chemistry. *J. Chem. Phys.*
