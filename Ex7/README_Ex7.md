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
