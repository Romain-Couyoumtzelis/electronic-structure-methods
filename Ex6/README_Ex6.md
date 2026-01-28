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
