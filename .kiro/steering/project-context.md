# Composipy — Project Context & Steering

## What this project is

Composipy is an open-source Python library for composite plate structural analysis and optimization in aerospace engineering. It was developed by Rafael Pereira da Silva as part of his Master's dissertation at ITA (Instituto Tecnológico de Aeronáutica), titled *"Gradient-Based Buckling Optimization of Composite Plates Combining a Semi-Analytical Model and Lamination Parameters"* (2024), advised by Prof. Flávio Bussamra and co-advised by Dr. Saullo Castro (TU Delft).

The library has **11,000+ downloads** on PyPI and is licensed under MIT.

- **GitHub:** https://github.com/rafaelpsilva07/composipy
- **Docs:** https://rafaelpsilva07.github.io/composipy
- **Version:** 1.6.0
- **PyPI:** `pip install composipy`

---

## Scientific Foundation

The library implements three interconnected methods from the thesis:

### 1. Rayleigh-Ritz with Bardell Shape Functions
Composite plate buckling is solved as a generalized eigenvalue problem `(K - λ·KG)·c = 0`. The shape functions are the **Bardell hierarchical polynomials** (Cubic Hermite Splines + Legendre orthogonal polynomials), which allow arbitrary boundary conditions to be configured simply by including or excluding specific basis function indices. This is the main differentiator from simple analytical solutions.

### 2. Lamination Parameters (LP)
Instead of designing in ply-angle space (discrete, combinatorially explosive), the library supports the **lamination parameter space**: four continuous parameters (ξ1…ξ4 for bending, or W1…W4 in the thesis notation) that encode the bending stiffness matrix D through material invariants (Tsai-Pagano formulation). The feasible region is bounded by a parabola `W3 ≥ 2W1² − 1`. This turns the problem into a differentiable continuous optimization.

### 3. Gradient-Based Optimization (SLSQP)
SciPy's SLSQP optimizer drives two objective functions:
- **maximize_buckling_load**: fixed thickness T, optimize (W1, W3) to maximize λ_crit
- **minimize_panel_weight**: optimize (T, W1, W3) to minimize T³ subject to λ_crit ≥ 1

The Lamination Database approach (Section 3.3 of the thesis) converts continuous LP optima back to discrete, manufacturable stacking sequences by nearest-neighbour lookup in a pre-built database.

---

## Architecture

```
OrthotropicMaterial / IsotropicMaterial   ← ply material properties + allowables
        ↓
LaminateProperty                          ← Classical Laminate Theory: ABD matrices,
                                             lamination parameters xiA/xiD
                                             supports BOTH angle-stacking and LP-dict construction
        ↓ (used by both)
PlateStructure          LaminateStrength
(Ritz-Bardell buckling) (ply-by-ply stress/strain/margins)
        ↓
  optimize/
  maximize_buckling_load()
  minimize_panel_weight()
```

```
composipy/
├── core/
│   ├── material.py        OrthotropicMaterial, IsotropicMaterial
│   ├── property.py        LaminateProperty (ABD, xiA, xiD, lamination params)
│   ├── structure.py       PlateStructure (Ritz buckling, mode shape plot)
│   └── strength.py        LaminateStrength (strains, stresses, max-stress margins)
├── optimize/
│   ├── _maximize_buckling.py   maximize_buckling_load()
│   ├── _minimize_panel_weight.py  minimize_panel_weight()
│   └── utils.py           Ncr_from_lp(), constraint helpers, plot_optimization()
├── nastranapi/
│   └── pcomp_generator.py  build_sequence() (LaTeX parser), build_pcomp() (PCOMP card)
├── pre_integrated_component/
│   ├── _S.py              Bardell shape functions as eval-strings dict
│   ├── _ii_F.py           Pre-integrated integrals as lookup dicts (the perf secret)
│   ├── functions.py       Wrappers: ii_ff, ii_fxi_fxi, etc.
│   ├── build_k.py         Stiffness matrix assembly: calc_k33_ijkl, calc_kG33_ijkl, etc.
│   ├── write_pre_integrated_terms.py  (generator script, run once with sympy)
│   └── write_shape_function.py        (generator script, run once with sympy)
└── utils/
    └── validators.py      ComposipyValidator: _float_number, _int_number, _is_instance
```

---

## Key Technical Details

### Pre-integration Scheme (Critical Performance Decision)
All integrals of Bardell basis function products over [-1, 1] are computed **symbolically once** (using sympy, via the write_* generator scripts) and stored as hardcoded numeric strings in `_ii_F.py`. At runtime, stiffness matrix entries are assembled using dictionary lookups + `eval()`. This avoids all numerical quadrature and gives sub-second buckling calculations (m=n=7 runs in ~0.3s vs. >1 hour with runtime symbolic integration).

Six integral types are stored:
- `ii_FF`: ∫f_i·f_k dξ
- `ii_FXI_F`: ∫f_i'·f_k dξ  
- `ii_FXI_FXI`: ∫f_i'·f_k' dξ
- `ii_FXIXI_F`: ∫f_i''·f_k dξ
- `ii_FXIXI_FXI`: ∫f_i''·f_k' dξ
- `ii_FXIXI_FXIXI`: ∫f_i''·f_k'' dξ

### Boundary Conditions
Configured by including/excluding Bardell basis function indices from the shape function expansion. Indices 0-3 = Cubic Hermite Splines (translations/rotations at each end); indices ≥4 = hierarchical Legendre polynomials (zero at edges). The `_compute_constraints()` method in `PlateStructure` builds the index sets accordingly. Supports mixed BCs (e.g., root panels clamped on one axis, simply-supported on the other — exactly as used in Case Study 3 of the thesis).

### Stacking Convention
Index 0 of the stacking list = **bottom-most ply** (first manufactured). This differs from some US textbooks. z_position[0] = -T/2.

### Lamination Parameters via LaminateProperty
Two construction modes:
1. `LaminateProperty([0, 45, -45, 90], ply)` — angle stacking, full ABD computed from Q_layup integration
2. `LaminateProperty({'xiD': [xi1, 0, xi3, 0], 'T': thickness}, ply)` — LP mode, D computed algebraically from invariants. Used by the optimizer for continuous optimization. B=0 (valid only for symmetric+balanced laminates).

### Sparse Eigensolver
`scipy.sparse.linalg.eigsh` with Cayley transform mode solves `KG·v = λ·K·v`. Eigenvalues recovered as `-1/λ`. Returns the `num_eigvalues` smallest critical load multipliers.

---

## Dependencies

- `numpy` — all matrix/array operations
- `scipy` — sparse eigensolver (`eigsh`), optimizer (`minimize` with `method='SLSQP'`)
- `matplotlib` — buckling mode shape plots, optimization contour plots
- `pandas` — stress/strain results as DataFrames
- Python 3.x (standard library: `itertools`, `time`, `warnings`)

Build: `setuptools` + `wheel`
Tests: `pytest`

---

## Testing

Run with: `pytest tests/`

Test files:
- `test_ABD.py` — ABD matrices vs. Mendonça (2005) textbook
- `test_lamination_parameters.py` — xiA, xiD values
- `test_laminate_strains.py` — ply strains/stresses vs. NASA TM-1995-009349
- `test_Material_Q0.py` — Q_0 matrix
- `test_plate_buckling.py` / `test_plate_buckling_LP.py` — eigenvalues for PINNED/CLAMPED plates
- `optimization_test/test_maximize_buckling_load.py` — vs. Gürdal & Haftka (1991/1993) NATO ASI results
- `optimization_test/test_minimize_panel_weight.py`

Tolerances: tight for pure algebra (machine precision), generous for buckling/optimization (rtol=0.1 — reflects the ~3% difference between Ritz and closed-form due to D16/D26 coupling terms which the analytical solution neglects).

---

## Thesis Results Reproduced by Composipy

- **Case Study 1**: Rayleigh-Ritz verified against NASTRAN SOL 105. Max difference 0.7% (PINNED), 0.4% (CLAMPED), 0.5% (SSSF). D16/D26 explain the ~3.5% gap vs. closed-form for quasi-isotropic laminates.
- **Case Study 2**: Reproduced Gürdal & Haftka (1991) results within 3.2%. Average loss of efficiency (continuous → discrete via Lamination Database): 8% overall, 16% for 8-ply, 9% for 24-ply.
- **Case Study 3**: Wing-box optimization (Liu, Haftka & Akgun 2000 benchmark). Total plies: 1082 vs. 1128 (Liu 2000), 1198 (Liu 2009), 1100 (Liu 2015). Full run time: 16 minutes on AMD Ryzen 7 3800X.

---

## Rafael's Goals & Vision

- This library is Rafael's scientific and engineering contribution, grown from his MSc thesis
- It is intended to be a **practical, engineer-friendly tool** for composite plate design in industry
- Future directions include: non-conventional ply angles, torsional spring BCs, blending constraints, Machine Learning dataset generation, Design of Experiments
- As an MIT-licensed open-source project, new contributors are welcome
- The library philosophy mirrors how pandas/numpy/scipy feel to use: composable, object-oriented, returns standard Python/numpy/pandas objects

---

## Code Style & Conventions

- Python 3, OOP-first (class hierarchy with inheritance from validator base)
- Properties are lazy-computed and cached (pattern: `self._X = None`, compute on first `@property` access)
- All public classes and functions have docstrings with Parameters / Returns / Examples sections (NumPy docstring style)
- `__all__` is defined in every module
- Validators are applied consistently in `__init__` methods via `ComposipyValidator` helpers
- Pre-integrated lookup tables use `eval()` on string expressions — do not refactor this without careful benchmarking
- Tests use `pytest`, reference published aerospace papers/textbooks for verification values
- CI: GitHub Actions runs tests on push/PR; separate workflow for PyPI publish

---

## Important: What NOT to break

1. The pre-integration lookup tables in `_ii_F.py` and `_S.py` — these were generated by symbolic computation and must not be hand-edited
2. The `eval()` pattern in `functions.py` — it is intentional for performance
3. The stacking convention (index 0 = bottom ply) — changing it would break all downstream CLT results
4. The LP-dict construction mode of `LaminateProperty` — used by the optimizer
5. The boundary condition index exclusion logic in `_compute_constraints()` — directly implements Bardell's formulation
