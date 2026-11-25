# Debris-Flow Runout Explorer

A browser-based 1D debris-flow runout simulator. It solves a depth-averaged, shallow-water-form system with Coulomb/Voellmy basal friction using a finite-volume method, and runs entirely in your browser: Python and NumPy execute client-side through [Pyodide](https://pyodide.org) in a Web Worker. There is no server and nothing to install.

**Live app:** https://lillian-hanson.github.io/gpgn436_debris_flow/

The app is a tool for building intuition about what controls runout (release volume, channel width, slope, friction). **It is not a hazard-prediction tool.** It has not been validated against field data, and its friction parameters cannot be identified from a runout distance (see [Model limits](#model-limits)).

Final project for GPGN 436/536 (Colorado School of Mines).

## What it simulates

A block of debris (3 m deep, length set by the release volume and channel width) is released at rest at the top of a channel. The channel has two reaches: a steep upper reach with a slope you set, and a 5° runout reach starting 350 m down the bed. The flow accelerates down the upper reach, decelerates on the runout reach, and deposits. The run ends when the front stops advancing.

Inputs: release volume V, channel width b, upper-reach slope θ₁, friction angle φ, turbulence coefficient ξ.
Outputs: side view and overhead schematic of the flow over time, depth profile, front position vs. time, and metrics computed from the solver output (runout, peak front speed, peak Froude number, maximum deposit thickness, mass-balance error).

### Governing equations

With x the distance along the bed, h the flow depth normal to the bed, and u the depth-averaged velocity:

```
mass:      h_t  + (hu)_x                           = 0
momentum: (hu)_t + (hu^2 + 0.5 g cos(θ) h^2)_x     = g h sin(θ) − τ_b/ρ
friction:  τ_b/ρ = μ g h cos(θ) sgn(u) + g u|u| / ξ,    μ = tan(φ)
```

Coulomb friction is the ξ → ∞ limit. The pressure term uses an isotropic earth-pressure coefficient (K = 1), so this is a shallow-water-form model with Coulomb/Voellmy friction, not the full Savage–Hutter model. The system is hyperbolic wherever h > 0, with wave speeds u ± √(g h cos θ).

The model is width-averaged: the channel width b only converts depth to volume (V = ∫ h b dx). The overhead view in the app is a schematic extrusion of the 1D solution, not a 2D simulation.

### Numerical method

- **Finite volume, flux-limited Lax-Wendroff with a Rusanov fallback.** At each cell interface the two-step (Richtmyer) Lax-Wendroff flux is blended with the Rusanov (local Lax-Friedrichs) flux:
  `F = F_Rusanov + φ (F_LW − F_Rusanov)`. A van Leer-type limiter sets φ ∈ [0, 1] from the smoothness of h and hu. φ = 1 (Lax-Wendroff) where the flow is smooth; φ = 0 (Rusanov) at shocks, extrema and next to dry cells.
- **Why limited:** pure Lax-Wendroff was tested and fails on this problem. It oscillates at the bore that forms at the slope break, drives the depth negative, and blew up at grid spacings of 0.5 m and 0.25 m.
- **Time stepping:** explicit, with Δt = 0.45 Δx / max(|u| + √(g h cos θ)) recomputed every step.
- **Source terms** by operator splitting: gravity (explicit), Coulomb friction as a projection toward u = 0 (so material at rest stays at rest, with no sign chatter), and Voellmy drag by backward Euler (stable as h → 0).
- **Dry front:** cells shallower than 10⁻¹⁰ m carry no momentum; their depth is kept so mass is conserved exactly.
- **Grid:** 1200 cells over 600 m (Δx = 0.5 m). Runout for the default inputs is 437 / 433 / 431 / 431 m at Δx = 2 / 1 / 0.5 / 0.25 m.

Every function in [`solver/core.py`](solver/core.py) has a docstring stating the equation it implements, the exact discretization, why that choice was made, units and array shapes, and its known weakness.

## Verification

The web app and the verification code import the same solver module ([`solver/`](solver/)); the scheme is never reimplemented for testing. The notebook [`verification/ritter_dam_break.ipynb`](verification/ritter_dam_break.ipynb) shows each test with plots and a pass/fail table, and [`tests/test_verification.py`](tests/test_verification.py) runs the same checks in CI before every deploy.

| Test | What it checks | Result |
|---|---|---|
| Ritter dry-bed dam break (exact solution) | flux, CFL, dry front | L2 error 5.3 × 10⁻³ at N = 3200; convergence rate rising 0.72 → 0.90 |
| Stoker wet-bed dam break (exact solution) | shock capturing (the limiter) | L1 rate ≈ 1; shock within 0.8 cells of the exact position; no oscillations (overshoot 1.9 × 10⁻⁴ m on a 0.30 m jump) |
| Uniform Voellmy flow (exact ODE solution) | source terms | terminal velocity to 2 × 10⁻¹⁴; transient error first order in Δt |
| Mass conservation | conservative update | relative error ≤ 6 × 10⁻¹⁶ |

Convergence rates stay below 1 in L2 because the front and the rarefaction-head kink are handled at first order; this is expected and explained in the notebook. Not covered by an exact solution: the bore at the slope break with friction acting, which is checked only by mass conservation and grid convergence of runout.

## Model limits

- Width-averaged 1D, so no lateral spreading or overbank flow.
- Constant effective friction, with no pore-pressure evolution and no entrainment.
- Single-phase mixture.
- The scheme drops to first order at the front (flux limiter), so the front is smeared.
- φ and ξ trade off against each other, so matching a runout distance does not identify either one.
- The φ and ξ ranges shown in the app come from the RAMMS::DEBRISFLOW manual, a model with yield-stress and curvature terms this one lacks; they are a guide, not calibrated values for this model.

## Repository layout

```
solver/                 NumPy-only solver, used by both the app and the verification
  core.py               finite-volume scheme: fluxes, limiter, update, sources, boundaries
  runout.py             two-reach bed, release block, stopping rule, metrics
tests/
  test_verification.py  Ritter, Stoker, uniform Voellmy flow, mass conservation
verification/
  ritter_dam_break.ipynb  verification notebook with plots
web/
  index.html, app.js, style.css   UI and plotting (main thread)
  worker.js             Web Worker: loads Pyodide 0.28.3 + NumPy, runs the solver
.github/workflows/pages.yml       runs the tests, then deploys web/ + solver/ to GitHub Pages
```

## Running locally

Requires Python 3 with NumPy. The notebook also uses SciPy, Matplotlib and Jupyter.

Run the verification tests:

```bash
python tests/test_verification.py
```

Re-run the verification notebook:

```bash
cd verification && jupyter nbconvert --to notebook --execute --inplace ritter_dam_break.ipynb
```

Serve the site with the same layout as GitHub Pages, then open http://localhost:8765/:

```bash
rm -rf _site && mkdir _site && cp -r web/* _site/ && mkdir _site/solver && cp solver/*.py _site/solver/ && python -m http.server 8765 -d _site
```

The app loads Pyodide and NumPy (about 10 MB) from the jsDelivr CDN on first load, so it needs an internet connection.

## Deployment

Pushing to `main` triggers [`.github/workflows/pages.yml`](.github/workflows/pages.yml), which runs the verification tests and, if they pass, publishes `web/` plus a copy of `solver/` to GitHub Pages. All asset paths are relative, so the site works under the `/gpgn436_debris_flow/` prefix. One-time setup: repository Settings → Pages → Source = "GitHub Actions".

## References

- Savage, S. B. and Hutter, K. (1989). The motion of a finite mass of granular material down a rough incline. *J. Fluid Mech.* 199, 177–215. doi:10.1017/S0022112089000340
- Voellmy, A. (1955). Über die Zerstörungskraft von Lawinen. *Schweizerische Bauzeitung* 73.
- George, D. L. and Iverson, R. M. (2014). A depth-averaged debris-flow model that includes the effects of evolving dilatancy. II. Numerical predictions and experimental tests. *Proc. R. Soc. A* 470, 20130820. doi:10.1098/rspa.2013.0820
- Lax, P. D. and Wendroff, B. (1960). Systems of conservation laws. *Commun. Pure Appl. Math.* 13.
- Richtmyer, R. D. (1963). A survey of difference methods for non-steady fluid dynamics. NCAR Technical Note 63-2.
- van Leer, B. (1974). Towards the ultimate conservative difference scheme II. *J. Comput. Phys.* 14.
- Rusanov, V. V. (1961). The calculation of the interaction of non-stationary shock waves with barriers. *USSR Comput. Math. Math. Phys.* 1(2), 304–320 (1962). doi:10.1016/0041-5553(62)90062-9
- Delestre, O. et al. SWASHES: a compilation of shallow water analytic solutions for hydraulic and environmental studies. arXiv:1110.0288 (Ritter and Stoker dam-break solutions).
- LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*. Cambridge University Press.
- Toro, E. F. *Riemann Solvers and Numerical Methods for Fluid Dynamics*. Springer.
- RAMMS AG / WSL Institute for Snow and Avalanche Research SLF (2024). *RAMMS::DEBRISFLOW User Manual* v1.8.0 (friction parameter guidance, §5.1).
