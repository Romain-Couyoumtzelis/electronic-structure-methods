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
