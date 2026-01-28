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
