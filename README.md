# Introduction to Electronic Structure Methods (IESM)

A comprehensive computational chemistry course covering quantum chemistry methods from basic linear algebra to transition state theory, using Psi4 as the primary software package.

---

# Exercise 1: Linear Algebra in Quantum Mechanics

This exercise reviews fundamental concepts from linear algebra and quantum mechanics that are essential for computational chemistry.

## Learning Goals

- Review basic concepts of linear algebra
- Review basic notation of quantum mechanics
- Basic vector operations using NumPy

## Notebooks

### IESM_Ex1.ipynb
Main introduction notebook that provides an overview of linear algebra concepts relevant to quantum mechanics. Contains links to slides and report templates.

### Ex1_LinearAlgebra.ipynb
Review of basic linear algebra operations:
- Vector products (cross product and scalar/dot product)
- Matrix multiplication
- Determinant evaluation
- Matrix exponentials and operator algebra

### Ex1_Numpy.ipynb
Introduction to working with vectors using NumPy:
- Representing wavefunctions as vectors
- Dirac notation (bra-ket) and its connection to vector operations
- Checking normalization and orthogonality of basis vectors
- Using `numpy` functions for inner products and complex conjugates

### Ex1_QM.ipynb
Comprehensive review of basic quantum mechanics concepts:
- **Operators**: Definition, eigenfunctions, and eigenvalues
- **Commutators**: Definition and connection to the Heisenberg uncertainty principle
- **Hilbert Space**: Bra-ket (Dirac) notation, state vectors, and basis expansions
- **Hermiticity**: Properties of Hermitian operators, linearity, and matrix representations
- **Wavefunctions**: Position representation, probability distributions, observables
- **Fermions and Bosons**: Symmetry constraints and Slater determinants
- **Schrödinger Equation**: Time-independent eigenvalue problem and its matrix formulation

## References

- Cohen-Tannoudji, C., Diu, B., & Laloe, F. (1986). *Quantum Mechanics, Volume 1*
  - Chapter II B-E and Complements

---

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

---

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

---

# Exercise 4: The Hartree-Fock Procedure in Detail

This exercise provides a hands-on implementation of the Hartree-Fock SCF procedure, building and diagonalizing the Fock matrix step by step while using Psi4 to evaluate the integrals.

## Learning Goals

- Build and diagonalize the Fock matrix
- Understand the steps in a Hartree-Fock SCF calculation
- Practice implicit summation (Einstein notation) in NumPy

## Notebooks

### IESM_Ex4b.ipynb
Main introduction notebook providing an overview, learning goals, and references.

### HF.ipynb
Complete step-by-step implementation of the Hartree-Fock procedure:

#### Part 1: Orthogonalizing the AO Basis Set
- Understanding overlap integrals and the overlap matrix S
- Building the overlap matrix from basis vectors: S = B†B
- Using Einstein summation notation (`np.einsum`)
- Constructing the orthogonalization matrix A = S^(-1/2)
- Transforming to an orthonormal basis: S' = ASA

#### Part 2: Gaussian AO Basis Set
- Setting up H₂O molecule in Psi4
- Understanding number of electrons and occupied orbitals
- Building the overlap matrix using `mints.ao_overlap()`
- Checking if the AO basis is orthonormal

#### Part 3: Initial Guess and Fock Matrix
- Using the core Hamiltonian (H = T + V) as initial guess
- Transforming Fock matrix to orthonormal basis: F' = AFA
- Diagonalizing F' to get eigenvalues and eigenvectors
- Transforming coefficients back to AO basis: C = A·C'

#### Part 4: Density Matrix
- Building the density matrix from occupied orbitals: D_pq = Σᵢ c*_pi·c_qi
- Understanding the physical meaning of D

#### Part 5: Coulomb and Exchange Integrals
- Building the electron repulsion integral tensor I_pqrs
- Computing Coulomb integral matrix: J_pq = Σ_rs D_rs·I_pqrs
- Computing Exchange integral matrix: K_ps = Σ_rq D_rq·I_pqrs
- Building the Fock matrix: F = H + 2J - K

#### Part 6: SCF Energy and Iteration
- Computing SCF energy: E = Σ_pq (H_pq + F_pq)·D_pq + E_nuc
- Implementing the full SCF loop with convergence criteria
- Comparing results with Psi4's built-in SCF
- Exploring the effect of convergence threshold on iterations

## Key Equations

- **Fock equation**: FC = SCE
- **Orthogonalization**: A†SA = 1
- **Density matrix**: D_pq = Σᵢ c*_pi·c_qi
- **Fock matrix**: F = H + 2J - K
- **SCF energy**: E = Σ_pq (H_pq + F_pq)·D_pq + E_nuc

## Data Files

- `ijk.dat` - Supporting data file for calculations

## References

- Sherrill Group: [Hartree-Fock Introduction](http://vergil.chemistry.gatech.edu/courses/chem6485/pdf/Hartree-Fock-Intro.pdf)
- Psi4Education: Original notebook by McDonald, Sirianni, Shepherd & Garrett-Roe

---

# Exercise 5: Post-Hartree-Fock Methods - CI and MPn

This exercise explores electron correlation and post-Hartree-Fock methods including Configuration Interaction (CI), Møller-Plesset Perturbation Theory (MPn), and Coupled Cluster (CC) methods.

## Learning Goals

- Understand what correlation energy is
- Get familiar with post-HF methods and their advantages/disadvantages
- Compare performance of different methods for correlation recovery and bond dissociation energies

## Notebooks

### IESM_Ex5.ipynb
Main introduction notebook explaining electron correlation and post-HF methods.

### Ex5_Boron.ipynb
**Recovering Correlation Energy: Boron Atom**

Comparing different methods for recovering correlation energy:
- Computing HF energy and comparing to "exact" experimental energy
- Running calculations with various post-HF methods:
  - **CISD** (Configuration Interaction Singles and Doubles)
  - **FCI** (Full Configuration Interaction)
  - **MP2, MP3, MP5** (Møller-Plesset Perturbation Theory)
  - **CCSD, CCSD(T)** (Coupled Cluster)
- Calculating percentage of correlation energy recovered by each method
- Understanding why some methods can "over-recover" correlation (basis set limitations)

### Ex5_BDE.ipynb
**Homolytic Cleavage of the C-F Bond**

Bond dissociation energy (BDE) calculations:
- Calculating BDE for CH₃F → CH₃• + F•
- Comparing HF, MP2, and MP3 results to experimental value (109.2 kcal/mol)
- Understanding the workflow: HF geometry optimization → post-HF single point
- Analyzing convergence of MPn series
- Discussing accuracy of different methods for radical systems

### Ex5_HNO3.ipynb
**Influence of Correlation on Geometry: HNO₃ Molecule**

Structural parameters with different methods:
- Optimizing HNO₃ geometry at HF, MP2, and MP3 levels
- Comparing computed bond lengths and angles to experimental values
- Understanding how electron correlation affects:
  - N-O bond lengths
  - O-N-O bond angles
- Discussing sources of discrepancy between theory and experiment

## Key Concepts

- **Correlation energy**: E_corr = E_exact - E_HF
- **Size consistency**: Property that E(A+B) = E(A) + E(B) for non-interacting systems
- **MPn convergence**: MP2 often overshoots, MP3 corrects back
- **Coupled Cluster**: Includes higher excitations through exponential ansatz
- **Basis set effects**: Finite basis limits the correlation that can be recovered

## Methods Covered

| Method | Type | Description |
|--------|------|-------------|
| HF | Mean-field | No correlation |
| MP2/MP3/MP5 | Perturbation | Perturbative correlation |
| CISD | Variational | Singles & Doubles excitations |
| FCI | Variational | Exact within basis |
| CCSD | Non-variational | Coupled Cluster S+D |
| CCSD(T) | Non-variational | "Gold standard" of quantum chemistry |

## References

- Sherrill, D. [An Introduction to Configuration Interaction Theory](http://vergil.chemistry.gatech.edu/notes/ci.pdf)
- Jensen, F. (2017). *Introduction to Computational Chemistry*. Chapter 4.8.1: Møller–Plesset perturbation theory

---

# Exercise 6: DFT vs (Post) HF Methods

This exercise compares Density Functional Theory (DFT) methods to wavefunction-based methods (HF, MP2) for predicting reaction enthalpies and geometric properties.

## Learning Goals

- Compare accuracy and efficiency of electron density-based methods to wavefunction-based methods
- Compare exchange-correlation functionals used in DFT calculations
- Learn how frontier orbital visualization supports electronic structure analysis

## Notebooks

### IESM_Ex6.ipynb
Main introduction notebook explaining the basics of DFT and the comparison with wavefunction methods.

### Ex6_aval.ipynb
**Methylcyclohexane A-value Calculation**

Computing conformational energy differences:
- Comparing equatorial vs axial methyl group conformations
- Downloading structures from PubChem
- Computing A-values with multiple methods:
  - Wavefunction-based: HF, MP2
  - DFT functionals: BLYP, MPW1PW91, B97-2, PBE, TPSS, M06-L, M06-2X, B3LYP-D3BJ
- Testing two basis sets: 6-31+G* and 6-31+G**
- Comparing to experimental A-value (1.74 kcal/mol)
- Analyzing computational cost vs accuracy trade-offs
- Understanding "chemical accuracy" target (<1 kcal/mol error)

### Ex6_geo.ipynb
**Geometric Properties: NO₃• Radical**

Predicting structural parameters of a challenging system:
- NO₃• radical: experimentally D₃ₕ symmetric (all N-O = 1.24 Å, all angles = 120°)
- Optimizing geometry with:
  - HF (with SOSCF for convergence)
  - MP2
  - Various DFT functionals (BLYP, BP86, PBE, B3LYP, B97-2, M06-L, M06-2X, TPSS, mPW1PW91, ωB97X-D)
- Visualizing frontier molecular orbitals (SOMO, HOMO, LUMO)
- Using `cubeprop` for orbital visualization
- Analyzing symmetry breaking in calculations
- Discussing static vs dynamic correlation

## Key Concepts

- **Exchange-correlation functional**: The unknown part of DFT that must be approximated
- **Jacob's Ladder of DFT**: LDA → GGA → meta-GGA → hybrid → double-hybrid
- **Kohn-Sham orbitals**: Auxiliary orbitals in DFT (useful for visualization despite lacking strict physical meaning)
- **A-value**: Energy difference between axial and equatorial conformers
- **SOSCF**: Second-order SCF for difficult convergence cases

## DFT Functionals Covered

| Functional | Type | Notes |
|------------|------|-------|
| BLYP | GGA | Becke exchange + LYP correlation |
| PBE | GGA | Parameter-free GGA |
| TPSS | meta-GGA | Includes kinetic energy density |
| B3LYP | Hybrid | 20% HF exchange |
| M06-2X | Hybrid meta-GGA | 54% HF exchange, good for main-group |
| B3LYP-D3BJ | Hybrid + Dispersion | Includes Grimme's D3 correction |

## References

- Sherrill, D. [Introduction to DFT (video)](https://www.youtube.com/watch?v=QGyfGCZT110)
- Becke, A.D. (2014). Perspective: Fifty years of density-functional theory in chemical physics. [DOI](https://doi.org/10.1063/1.4869598)

---

# Exercise 7: Troubleshooting, Pitfalls, and Traps

This exercise covers common errors in electronic structure calculations, dispersion corrections in DFT, and the importance of integration grids.

## Learning Goals

- Become familiar with common errors in computational chemistry calculations
- Determine which systems require dispersion corrections with DFT
- Evaluate how integration grid size affects calculation accuracy

## Notebooks

### IESM_Ex7.ipynb
Main introduction notebook providing an overview of common issues and troubleshooting strategies.

### Ex7_errors.ipynb
**Fixing Errors in Calculations**

Debugging common computational chemistry errors:
- **Allyl radical**: Fixing spin multiplicity and reference issues
- **Cu[H₂O]²⁺ complex**: Handling transition metal calculations
- **N₂ geometry optimization**: Fixing symmetry and convergence issues
- **³O single point**: Setting correct spin state and reference

Key error types covered:
- Incorrect charge/multiplicity
- Wrong reference (RHF vs UHF vs UKS)
- SCF convergence failures
- Symmetry-related issues

### Ex7_DFT_Hard_Easy.ipynb
**Hard and Easy Cases for DFT**

Understanding when DFT succeeds or fails:

**Dispersion Correction Section:**
- London dispersion interactions and their origin
- Why standard DFT functionals miss dispersion
- DFT-D3 empirical corrections (Grimme)
- Non-local correlation functionals (VV10, vdW-DF)

**Three Test Cases:**
1. **CaO bond dissociation** (Hard): Multireference character
2. **Ethane dimer interaction** (Medium): Dispersion-dominated
3. **Methanol eclipsed/staggered** (Easy): Conformational energy

Methods compared: PBE, B3LYP, MN15, MP2, SCAN, with/without dispersion

Key insights:
- Dispersion corrections (D3, NL) crucial for non-bonded interactions
- Some functionals trained on specific datasets perform well on those systems
- Error cancellation can artificially improve results

### Ex7_DFT_integration_grid.ipynb
**Integration Grids in DFT**

Understanding numerical integration quality:
- DFT requires numerical integration for XC terms
- Grid defined by radial and angular points (e.g., 75,302 or 99,590)
- Standard grids (SG-0, SG-1) may be insufficient

**Test case**: Argon dimer dissociation curve (3-6 Å)

Functionals tested:
- **B3LYP**: Well-behaved with coarse grids
- **M06-HF**: Sensitive to grid quality, shows artifacts

Grid options compared:
- SG1 (small, pruned)
- 75,302 (medium)
- 99,590 (fine)

Key findings:
- Meta-GGA functionals often need finer grids
- PES curves for weak interactions are grid-sensitive
- Always check grid convergence for production calculations

## Key Concepts

- **Dispersion**: Long-range electron correlation, decays as R⁻⁶
- **Integration grid**: Numerical quadrature for XC integrals
- **Chemical accuracy**: Error < 1 kcal/mol
- **Reference wavefunction**: RHF/UHF/ROHF for HF; RKS/UKS/ROKS for DFT

## References

- Morgante, P. & Peverati, R. (2020). The devil is in the details. [DOI](https://doi.org/10.1002/qua.26332)
- Bursch, M. et al. (2022). Best-Practice DFT Protocols. [DOI](https://doi.org/10.1002/ange.202205735)

---

# Exercise 8: Potential Energy Scans and Trajectory Visualization

This exercise explores potential energy surface (PES) scans along specific degrees of freedom, using butane dihedral rotation as the main example.

## Learning Goals

- Differentiate between rigid and relaxed PES scans
- Connect results of PES scans to known molecular properties
- Visualize conformational changes along reaction coordinates

## Notebooks

### IESM_Ex8.ipynb
Main introduction notebook explaining PES scan concepts and methodology.

### Ex8.ipynb
**Dihedral Scan of Butane**

Complete workflow for scanning a potential energy surface:

#### Part 1: Setting Up the Scan
- Building butane geometry using Z-matrix format
- Understanding internal coordinates (bonds, angles, dihedrals)
- Geometry optimization as starting point

#### Part 2: Rigid vs Relaxed Scans
- **Rigid scan**: Only the scanned coordinate changes, all others fixed
- **Relaxed scan**: Other coordinates optimized at each point (constrained optimization)
- Trade-offs between computational cost and accuracy

#### Part 3: Performing the Scan
- Creating geometries at different dihedral angles (180° to -180°)
- Single-point energy calculations at each geometry
- Using Python loops and f-strings for automated scans

#### Part 4: Analysis and Visualization
- Plotting the 1D potential energy profile
- Identifying minima and maxima (conformers and barriers)
- Interactive 3D visualization with sliders
- Monitoring geometric parameters (bond lengths) during the scan

#### Part 5: Conformational Analysis
- **Anti conformation** (180°): Global minimum, staggered
- **Gauche conformations** (±60°): Local minima
- **Eclipsed conformations** (0°, ±120°): Energy maxima
- Newman projections and organic chemistry connections

## Key Results

| Conformation | Dihedral | Relative Energy |
|--------------|----------|-----------------|
| Anti | 180° | 0 (reference) |
| Gauche | ±60° | ~0.9 kcal/mol |
| Eclipsed (CH₃-CH₃) | 0° | ~5-6 kcal/mol |
| Eclipsed (H-H) | ±120° | ~3-4 kcal/mol |

## Key Concepts

- **Potential Energy Surface (PES)**: Energy as function of nuclear coordinates (3N-6 dimensions)
- **Dihedral angle**: Torsion angle between four atoms
- **Rotational barrier**: Energy required to rotate around a bond
- **Conformational isomers**: Same connectivity, different 3D arrangement

## Practical Tips

- Use `psi4.optimize()` for initial structure optimization
- Use f-strings to parameterize Z-matrix templates
- Store results in dictionaries/lists for analysis
- Convert energies: `psi4.constants.hartree2kcalmol`
- Use `np.linspace()` for evenly spaced scan points

## References

- Murcko, M.A., Castejon, H., Wiberg, K.B. (1996). *J. Phys. Chem.* 100, 16162

---

# Exercise 9: Finding Transition States and Barrier Heights

This exercise covers transition state theory, including locating first-order saddle points and computing reaction barriers for the cyclization of deprotonated chloropropanol to propylene oxide.

## Learning Goals

- Understand how to navigate the PES to transition states
- Visualize chemical reactions and reaction pathways
- Connect transition state structures to reaction mechanisms

## Notebooks

### IESM_Ex9.ipynb
Main introduction notebook explaining transition state theory and saddle point optimization.

### Ex9.ipynb
**Locating Transition States: Constrained Optimizations**

Complete workflow for finding transition states:

#### Part 1: Constructing a TS Guess
- Starting from reactant geometry (chloropropanoate)
- Modifying Z-matrix to approximate TS:
  - Elongating C-Cl bond (1.79 → 2.40 Å)
  - Reducing O-C-C angle (110° → 80°) for ring formation
  - Removing H from oxygen (deprotonation)
- Physical intuition for TS structure

#### Part 2: Constrained Optimization
- Freezing the reaction coordinate (C-Cl distance)
- Relaxing all other degrees of freedom
- Using `frozen_distance` option in Psi4
- Importance of being "close enough" to the TS

#### Part 3: Transition State Optimization
- Setting `opt_type: ts` for saddle point search
- Computing Hessian with `full_hess_every: 0`
- Verifying TS by checking for exactly one imaginary frequency

#### Part 4: Normal Mode Analysis
- Visualizing vibrational modes with `show_normal_modes()`
- Identifying the imaginary mode (reaction coordinate)
- Understanding low vs high frequency motions
- Writing normal modes to Molden format

#### Part 5: Intrinsic Reaction Coordinate (IRC)
Additional notebooks explore:
- **Forward.ipynb**: Following IRC toward products
- **Backward.ipynb**: Following IRC toward reactants
- **Visualize.ipynb**: Animating the reaction pathway

## Reaction Studied

```
Chloropropanoate (deprotonated) → Propylene oxide + Cl⁻
```

An intramolecular Sₙ2 reaction forming an epoxide ring.

## Key Concepts

- **Transition State (TS)**: First-order saddle point on PES
- **Imaginary frequency**: Indicates unstable mode (reaction coordinate)
- **Constrained optimization**: Fix some coordinates, optimize others
- **Early vs Late TS**: Hammond postulate - TS resembles higher energy species
- **Intrinsic Reaction Coordinate (IRC)**: Minimum energy path connecting TS to minima

## TS Search Workflow

1. **Guess construction**: Modify reactant geometry toward expected TS
2. **Constrained optimization**: Fix reaction coordinate, relax rest
3. **Hessian calculation**: Check for one negative eigenvalue
4. **TS optimization**: Walk uphill along reaction coordinate
5. **Verification**: Confirm single imaginary frequency
6. **IRC calculation**: Verify TS connects correct reactants/products

## Psi4 Options for TS Search

```python
psi4.set_options({
    "opt_type": "ts",           # Optimize to saddle point
    "geom_maxiter": 500,        # Max optimization steps
    "full_hess_every": 0,       # Compute initial Hessian
    "normal_modes_write": True, # Save normal modes
    "frozen_distance": "1 3"    # Freeze distance between atoms 1 and 3
})
```

## References

- Jensen, F. (2017). *Introduction to Computational Chemistry*. Chapter 12.8, p.416ff
