# Introduction to Electronic Structure Methods (IESM)

A comprehensive computational chemistry course covering quantum chemistry methods from basic linear algebra to transition state theory, using Psi4 as the primary software package.

## Course Overview

This repository contains 9 exercises progressing from fundamental quantum mechanics concepts to advanced topics like DFT troubleshooting and transition state searches.

---

## Exercise 1: Linear Algebra in Quantum Mechanics

**Folder:** [Ex1/](Ex1/)

### Learning Goals
- Review basic concepts of linear algebra
- Review basic notation of quantum mechanics
- Basic vector operations using NumPy

### Notebooks

| Notebook | Description |
|----------|-------------|
| IESM_Ex1.ipynb | Main introduction with links to slides and resources |
| Ex1_LinearAlgebra.ipynb | Vector products, matrix multiplication, determinants, matrix exponentials |
| Ex1_Numpy.ipynb | Wavefunctions as vectors, Dirac notation, normalization, orthogonality |
| Ex1_QM.ipynb | Operators, eigenfunctions, commutators, Hermiticity, Hilbert space, Schrödinger equation |

### Key Concepts
- Bra-ket (Dirac) notation
- Hermitian operators and observables
- Slater determinants for fermions

---

## Exercise 2: First Steps in Psi4 - Basis Sets and Hartree-Fock

**Folder:** [Ex2/](Ex2/)

### Learning Goals
- Understand what basis functions are
- Get familiar with basis function notation
- Learn the basics of Psi4
- Get familiar with molecular geometries (Z-matrix and Cartesian)

### Notebooks

| Notebook | Description |
|----------|-------------|
| IESM_Ex2.ipynb | Main introduction to Psi4 and basis sets |
| Ex2_theory.ipynb | LCAO, STOs, GTOs, Pople-type basis sets |
| Ex2_hydrogen_atom.ipynb | First Psi4 tutorial with H and He atoms |
| Ex2_water-molecule.ipynb | H₂O calculations, effect of basis sets, molecular geometry |

### Key Concepts
- Single-point calculations
- SCF convergence
- Basis set hierarchy: Minimal → Split-valence → Polarized → Diffuse

---

## Exercise 3: Basis Sets, Dissociation Energy, and Geometry Optimization

**Folder:** [Ex3/](Ex3/)

### Learning Goals
- Understand the influence of basis sets on energy calculations
- Learn the geometry optimization procedure
- Understand the difference between RHF and UHF methods

### Notebooks

| Notebook | Description |
|----------|-------------|
| IESM_Ex3.ipynb | Overview of dissociation and optimization |
| Ex3_H2molecule.ipynb | Effect of basis set size on H₂ energy |
| Ex3_Dissociation.ipynb | RHF vs UHF dissociation curves |
| Ex3_GeometryOpt.ipynb | H₂O geometry optimization on PES |
| optimizations1/2.ipynb | Additional optimization examples |

### Key Concepts
- Potential Energy Surface (PES)
- RHF vs UHF: spin-restricted vs unrestricted
- Convergence criteria: ΔE, forces, displacements

---

## Exercise 4: The Hartree-Fock Procedure in Detail

**Folder:** [Ex4/](Ex4/)

### Learning Goals
- Build and diagonalize the Fock matrix
- Understand the steps in a HF SCF calculation
- Practice Einstein summation notation in NumPy

### Notebooks

| Notebook | Description |
|----------|-------------|
| IESM_Ex4b.ipynb | Introduction to implementing HF |
| HF.ipynb | Step-by-step HF implementation using Psi4 integrals |

### Key Equations
```
FC = SCE                    (Fock equation)
F = H + 2J - K              (Fock matrix)
D_pq = Σᵢ c*_pi·c_qi        (Density matrix)
E = Σ_pq (H_pq + F_pq)·D_pq (SCF energy)
```

### Implementation Steps
1. Build overlap matrix S
2. Construct orthogonalization matrix A = S^(-1/2)
3. Initial guess: F = H (core Hamiltonian)
4. Transform and diagonalize: F' = AFA
5. Build density matrix from occupied orbitals
6. Compute J and K from electron repulsion integrals
7. Iterate until convergence

---

## Exercise 5: Post-Hartree-Fock Methods - CI and MPn

**Folder:** [Ex5/](Ex5/)

### Learning Goals
- Understand what correlation energy is
- Get familiar with post-HF methods and their advantages/disadvantages

### Notebooks

| Notebook | Description |
|----------|-------------|
| IESM_Ex5.ipynb | Introduction to electron correlation |
| Ex5_Boron.ipynb | Correlation energy recovery: HF, MP2-5, CISD, FCI, CCSD(T) |
| Ex5_BDE.ipynb | Bond dissociation energy: CH₃F → CH₃• + F• |
| Ex5_HNO3.ipynb | Effect of correlation on HNO₃ geometry |

### Methods Comparison

| Method | Type | Description |
|--------|------|-------------|
| HF | Mean-field | No correlation |
| MP2/MP3/MP5 | Perturbation | Perturbative correlation |
| CISD | Variational | Singles & Doubles excitations |
| FCI | Variational | Exact within basis |
| CCSD(T) | Non-variational | "Gold standard" |

### Key Concepts
- Correlation energy: E_corr = E_exact - E_HF
- Size consistency
- MPn convergence behavior

---

## Exercise 6: DFT vs (Post) HF Methods

**Folder:** [Ex6/](Ex6/)

### Learning Goals
- Compare accuracy and efficiency of DFT vs wavefunction methods
- Compare exchange-correlation functionals
- Learn frontier orbital visualization

### Notebooks

| Notebook | Description |
|----------|-------------|
| IESM_Ex6.ipynb | Introduction to DFT |
| Ex6_aval.ipynb | Methylcyclohexane A-value calculation |
| Ex6_geo.ipynb | NO₃• radical geometry prediction |

### DFT Functionals Covered

| Functional | Type | Notes |
|------------|------|-------|
| BLYP | GGA | Becke exchange + LYP correlation |
| PBE | GGA | Parameter-free |
| TPSS | meta-GGA | Includes kinetic energy density |
| B3LYP | Hybrid | 20% HF exchange |
| M06-2X | Hybrid meta-GGA | 54% HF exchange |
| B3LYP-D3BJ | Hybrid + Dispersion | With Grimme's D3 |

### Key Concepts
- Jacob's Ladder of DFT
- Kohn-Sham orbitals
- Chemical accuracy target: < 1 kcal/mol

---

## Exercise 7: Troubleshooting, Pitfalls, and Traps

**Folder:** [Ex7/](Ex7/)

### Learning Goals
- Become familiar with common calculation errors
- Determine when dispersion corrections are needed
- Evaluate integration grid effects

### Notebooks

| Notebook | Description |
|----------|-------------|
| IESM_Ex7.ipynb | Overview of common issues |
| Ex7_errors.ipynb | Debugging: allyl radical, Cu-H₂O, N₂, ³O |
| Ex7_DFT_Hard_Easy.ipynb | Hard/easy cases: CaO, ethane dimer, methanol |
| Ex7_DFT_integration_grid.ipynb | Grid convergence: Ar₂ with B3LYP vs M06-HF |

### Common Error Types
- Incorrect charge/multiplicity
- Wrong reference (RHF vs UHF vs UKS)
- SCF convergence failures
- Insufficient integration grid

### Key Concepts
- Dispersion corrections (D3, NL)
- Integration grid quality (SG1 vs 75,302 vs 99,590)
- Meta-GGA functionals need finer grids

---

## Exercise 8: Potential Energy Scans and Trajectory Visualization

**Folder:** [Ex8/](Ex8/)

### Learning Goals
- Differentiate between rigid and relaxed PES scans
- Connect PES scans to molecular properties

### Notebooks

| Notebook | Description |
|----------|-------------|
| IESM_Ex8.ipynb | Introduction to PES scans |
| Ex8.ipynb | Butane dihedral scan (180° to -180°) |

### Butane Conformations

| Conformation | Dihedral | Relative Energy |
|--------------|----------|-----------------|
| Anti | 180° | 0 (reference) |
| Gauche | ±60° | ~0.9 kcal/mol |
| Eclipsed (CH₃-CH₃) | 0° | ~5-6 kcal/mol |

### Key Concepts
- Rigid vs relaxed scans
- Z-matrix parameterization
- Newman projections
- Rotational barriers

---

## Exercise 9: Finding Transition States and Barrier Heights

**Folder:** [Ex9/](Ex9/)

### Learning Goals
- Navigate the PES to transition states
- Visualize chemical reactions

### Notebooks

| Notebook | Description |
|----------|-------------|
| IESM_Ex9.ipynb | Introduction to transition state theory |
| Ex9.ipynb | TS search for chloropropanol → propylene oxide |
| Forward/Backward.ipynb | IRC calculations |
| Visualize.ipynb | Reaction pathway animation |

### TS Search Workflow
1. **Guess construction**: Modify geometry toward expected TS
2. **Constrained optimization**: Fix reaction coordinate
3. **Hessian calculation**: Check for one negative eigenvalue
4. **TS optimization**: `opt_type: ts`
5. **Verification**: Single imaginary frequency
6. **IRC calculation**: Verify correct reactants/products

### Key Psi4 Options
```python
psi4.set_options({
    "opt_type": "ts",
    "full_hess_every": 0,
    "frozen_distance": "1 3"
})
```

### Key Concepts
- First-order saddle points
- Imaginary frequencies
- Early vs late transition states (Hammond postulate)
- Intrinsic Reaction Coordinate (IRC)

---

## Software Requirements

- **Psi4**: Quantum chemistry package
- **Python 3**: With NumPy, Matplotlib, Pandas
- **py3Dmol**: Molecular visualization
- **Jupyter**: Notebook environment

## References

- Jensen, F. (2017). *Introduction to Computational Chemistry*. John Wiley & Sons.
- Cohen-Tannoudji, C., Diu, B., & Laloe, F. (1986). *Quantum Mechanics*.
- Sherrill Group resources: [vergil.chemistry.gatech.edu](http://vergil.chemistry.gatech.edu)
- Psi4 documentation: [psicode.org](https://psicode.org)
