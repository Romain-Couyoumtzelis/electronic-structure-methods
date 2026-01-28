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
