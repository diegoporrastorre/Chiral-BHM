# Initial numerical workflow

The equations and conventions are in `derivations/mean_field_bogoliubov.md`.
`src/mean_field.py` implements the model, analytic Bogoliubov matrix, and
a local stationary root solver with its analytic real Jacobian.

In a Python environment with the packages in `requirements.txt`, run:

```powershell
python src/validate_mean_field.py
```

For the project-local dependencies installed during initial development:

```powershell
& 'C:\Users\usuario\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -c "import sys; sys.path.insert(0, '.local-deps'); sys.path.insert(0, 'src'); import runpy; runpy.run_path('src/validate_mean_field.py', run_name='__main__')"
```

Checks cover finite-difference linearization for all three boundary/profile
choices, symmetry of the pairing block, Hermiticity of the conservative
normal block, homogeneous PBC stationary densities and momentum-space
spectra for V = 0 and V = U, and an exactly solvable linear driven chain.

The default damped width of five sites is an initial configurable choice,
not a verified parameter from the paper. Construct a new Chain when changing
parameters, since hopping and profile arrays are built on initialization.

This is a foundation, not a completed phase-diagram calculation. A single
local root solve does not find all branches. Next use multiple seeds and
continuation in both sweep directions, recording residuals and growth rates,
and compare OBC profiles with PBC. Test selected solutions by time evolution.
Refine system sizes and tolerances before interpreting phase boundaries.
Topological diagnostics have not yet been implemented.

## Initial validation (2026-11-03)

All checks passed using NumPy 2.5.4 and SciPy 1.18.1. Maximum directional
finite-difference error was 6.3e-10; the real-space and momentum-space
PBC spectra agreed within 2.5e-14. The local installation required running
outside the sandbox because the sandbox could not read its package files.

At N=12, J=1, Delta=0.5, U=-0.0002, kappa=1, phi=pi/3, epsilon=38,
the homogeneous PBC density equation at V=0 has three branches:

| Density | Maximum Im(omega) | Linear stability |
| --- | --- | --- |
| 974.003409 | -0.178868 | Stable |
| 2088.609970 | +0.317738 | Unstable |
| 4436.386621 | -0.500000 | Stable |

At V=U the same benchmark has one uniform branch, density 2572.261005,
with maximum Im(omega) = -0.5. These are small-system validation points,
not a phase diagram or evidence about topology.
