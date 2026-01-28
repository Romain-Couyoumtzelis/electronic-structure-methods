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
