# CLAUDE.md

## Project

GPGN 436/536 (Colorado School of Mines) final project: a browser-based 1D debris-flow runout simulator. It solves a depth-averaged, shallow-water-form hyperbolic PDE system with Voellmy friction (Savage-Hutter-type, but with isotropic pressure; see References and provenance) using a finite-volume method. Python runs client-side via Pyodide in a Web Worker. Deployed as a static site on GitHub Pages.

Purpose: a zero-install tool for building intuition about what controls runout (volume, width, slope, friction). It is NOT a hazard-prediction tool. Do not write UI copy, docstrings, or README text that claims predictive or hazard-zoning capability.

Deployment requirements (from the assignment):
- Use Pyodide in a Web Worker to run Python directly in the browser.
- Deploy the website on GitHub Pages using GitHub Actions, without a backend server.
- The final submission is a URL that can directly run the software in a browser. The deployed Pages URL must work with no install, login, or local setup.

Course rubric constraints:
- The numerical method must be one learned in this course. Do not swap in a scheme (WENO, HLLC, Godunov with exact Riemann solver, etc.) without asking.
- The solver must be verified against an analytic solution, and the verification must test the same code the app runs.
- A focused 1D app is acceptable. 2D is out of scope.

## Governing equations

One coupled hyperbolic system, U_t + F(U)_x = S(U), with U = [h, hu].

- Mass: h_t + (hu)_x = 0
- Momentum: (hu)_t + (hu^2 + 0.5 g cos(theta) h^2)_x = g h sin(theta) - tau_b/rho
- Friction (Voellmy): tau_b/rho = mu g h cos(theta) sgn(u) + g u|u| / xi, with mu = tan(phi)
- Coulomb is the xi -> infinity limit. Do not add a separate Coulomb code path.

Variables: x = distance along the bed (m); h = flow depth normal to the bed (m); u = depth-averaged velocity along the bed (m/s), a SOLVED variable and never an input; theta = bed slope (rad internally, degrees in the UI); phi = effective friction angle; xi = turbulence coefficient (m/s^2); g = 9.81 m/s^2. Density cancels.

Wave speeds: lambda = u +/- sqrt(g h cos(theta)). The system is hyperbolic wherever h > 0 and degenerates at the dry front (h -> 0).

## References and provenance

The equations above were written from the standard depth-averaged shallow-flow formulation. They were not copied from the course notes or from one paper. Cite these in docstrings, the notebook, and the presentation. Items marked (unverified) were not checked against the source.

Cited sources:
- Savage, S. B. and Hutter, K. (1989). The motion of a finite mass of granular material down a rough incline. J. Fluid Mech. 199, 177-215. doi:10.1017/S0022112089000340. Depth-averaged equations with a Coulomb-like basal friction law, resembling the nonlinear shallow-water equations.
- Voellmy, A. (1955). Ueber die Zerstoerungskraft von Lawinen. Schweizerische Bauzeitung 73, 159-162, 212-217, 246-249, 280-285. Basal resistance = dry Coulomb friction + a term proportional to u^2 ("turbulent" friction). The xi notation is from the Voellmy-Salm model (Salm, Burkard and Gubler, 1990, Swiss avalanche guidelines; full reference unverified).
- Iverson, R. M. and George, D. L. (2014). A depth-averaged debris-flow model that includes the effects of evolving dilatancy. I. Physical basis. Proc. R. Soc. A 470, 20130819. doi:10.1098/rspa.2013.0819. Generalized Coulomb basal friction (D-Claw).
- George, D. L. and Iverson, R. M. (2014). A depth-averaged debris-flow model that includes the effects of evolving dilatancy. II. Numerical predictions and experimental tests. Proc. R. Soc. A 470, 20130820. doi:10.1098/rspa.2013.0820. Eigenstructure analysis showing the system is hyperbolic. Use this to support "debris flow is hyperbolic."
- Rusanov, V. V. (1961). The calculation of the interaction of non-stationary shock waves with barriers. Zh. Vychisl. Mat. Mat. Fiz. 1(2), 267-279. English: USSR Comput. Math. Math. Phys. 1(2), 304-320 (1962). doi:10.1016/0041-5553(62)90062-9. Builds on Lax, P. D. (1954), Weak solutions of nonlinear hyperbolic equations and their numerical computation, Commun. Pure Appl. Math. 7, 159.
- Lax, P. D. and Wendroff, B. (1960). Systems of conservation laws. Commun. Pure Appl. Math. 13, 217-237 (unverified: confirm pages and doi). The Lax-Wendroff scheme.
- Richtmyer, R. D. (1963). A survey of difference methods for non-steady fluid dynamics. NCAR Technical Note 63-2 (unverified). The two-step form of Lax-Wendroff used here; also in Richtmyer and Morton, Difference Methods for Initial-Value Problems (1967).
- van Leer, B. (1974). Towards the ultimate conservative difference scheme II. Monotonicity and conservation combined in a second-order scheme. J. Comput. Phys. 14, 361-370 (unverified). The limiter function psi(r) = (r + |r|)/(1 + |r|).
- RAMMS AG / WSL Institute for Snow and Avalanche Research SLF (2024). RAMMS::DEBRISFLOW User Manual v1.8.0, July 2024 (read 2026-11-25). Eq. 3.1, p. 16: S = mu N + rho g u^2 / xi, N = rho h g cos(slope), the same convention as this model. Section 5.1, p. 83: mu normally 0.05-0.4, values above 0.4 rarely useful; Table 5.3: xi 100-200 m/s^2 for granular (solid-dominated) flows, 200-1000 m/s^2 for muddy (fluid-like) flows, "only suggestions, not fixed definitions"; Table 5.4: mu = 0.2, xi = 200 m/s^2 if the flow type is unknown; Table 5.5, p. 84: Dorfbach best fit mu = 0.225, xi = 130 m/s^2. Gives no values for rock avalanches or snow. RAMMS adds a yield stress and terrain curvature that this model lacks.
- Ritter dry-bed dam-break formulas, as compiled in SWASHES: a compilation of shallow water analytic solutions for hydraulic and environmental studies, arXiv:1110.0288.
- Ritter, A. (1892). Die Fortpflanzung der Wasserwellen. Z. Ver. Dtsch. Ing. 36 (unverified: confirm volume and pages before citing).
- Toro, E. F., Riemann Solvers and Numerical Methods for Fluid Dynamics (Springer) and LeVeque, R. J., Finite Volume Methods for Hyperbolic Problems (Cambridge, 2002) (both unverified). Either is a single source for the flux, CFL, and source-term treatment. LeVeque ch. 6 covers flux-limited (high-resolution) Lax-Wendroff methods.

Derived here, not cited (show the derivation in the notebook):
- Wave speeds lambda = u +/- sqrt(g h cos(theta)) from the eigenvalues of the flux Jacobian.
- CFL time step dt = C dx / max(|u| + c).
- Terminal velocity u_t = sqrt(xi h (sin(theta) - mu cos(theta))) from the steady balance g h sin(theta) = mu g h cos(theta) + g u^2 / xi.
- Placement and units of xi: the Voellmy drag is written as g u|u| / xi per unit area (after dividing by rho) in the hu equation, and g u|u| / (xi h) in the velocity equation. Voellmy conventions for xi differ between papers, so state this one explicitly.

Known gaps:
- Pressure term: from memory (unverified), the original Savage-Hutter equations multiply 0.5 g cos(theta) h^2 by an earth-pressure coefficient K. This model sets K = 1 (isotropic pressure). Call it a "shallow-water-form model with Coulomb/Voellmy friction," not "Savage-Hutter."
- Course coverage: Dr. Jin taught Lax-Wendroff by name (per Lillian, 2026-11-23), so Lax-Wendroff is the course method. The scheme is a flux-limited Lax-Wendroff, and two of its parts were NOT confirmed as course content: the van Leer-type flux limiter and the Rusanov (local Lax-Friedrichs) flux it falls back to. Pure Lax-Wendroff was tested and fails on this problem (h reached -27 m at dx = 1 m; NaN at dx = 0.5 m and 0.25 m), so these additions are needed, not optional. Describe the scheme as "flux-limited Lax-Wendroff with a Rusanov fallback," never as plain "Lax-Wendroff," and do not state in code, comments, or slides that the limiter or Rusanov was covered in class unless Dr. Jin confirms it.

Cross-check (run 2026-11-19; symbolic with SymPy, numeric with a prototype NumPy solver):
- Ritter solution satisfies h_t + (hu)_x = 0 and (hu)_t + (hu^2 + 0.5 g h^2)_x = 0 exactly (residual 0). It matches h = h0 at x = -c0 t, u = 0 at the head, and h = 0 with front speed 2 c0 at x = 2 c0 t.
- Flux Jacobian eigenvalues are u +/- sqrt(g h cos(theta)).
- The conservative momentum equation minus u times the mass equation equals h (u_t + u u_x + g cos(theta) h_x), confirming the velocity form.
- Steady Voellmy balance gives u^2 = xi h (sin(theta) - mu cos(theta)). All terms in the hu equation have units m^2/s^2.
- Terminal Froude number is independent of h: Fr^2 = xi (sin(theta) - mu cos(theta)) / (g cos(theta)).

## Design decisions (current; change only after discussion)

- Width-averaged 1D along the channel. Channel width b is a parameter: V = integral of h*b dx. There is no cross-channel dimension, so no overflow, superelevation, or lateral spreading. The "overhead view" in the UI is a schematic extrusion of the 1D solution and must be labeled as such.
- Two-reach bed: steep upper reach (user-set theta_1) and a fixed 5 degree runout reach starting at s = 350 m. Without a slope break, a constant slope with phi < theta never produces a deposit. Because cos(theta) jumps at the break, the pressure flux needs care there (smooth the transition or discretize consistently). Flag this if you touch the flux.
- Release: initial block of depth h0 = 3 m with length L0 = V / (b h0), at rest (u = 0). Wall at the upstream end, open outflow at the downstream end.
- Grain size and "debris type" enter only through phi and xi. Do not add other material inputs.
- Material context is text only (decided 2026-11-25): the phi and xi sliders carry the RAMMS guidance ranges as annotations, with the citation in the model-limits panel. No material presets or dropdown. Named materials ("mud", "gravel", "rock slide") with fixed phi and xi have no source and would imply that material identifies friction, which contradicts the phi-xi trade-off limit. The only sourced split is granular vs. muddy, and it is a split in xi alone. The phi slider stays at 6-40 degrees, wider than the RAMMS range (about 3-22 degrees), so that phi >= theta_1 (material that does not move) is reachable across the theta_1 range; the annotation states the RAMMS range instead (decided 2026-11-25).
- Fixed constants: g, h0, domain length 600 m, CFL = 0.45.

## Numerics

- Finite volume, flux-limited Lax-Wendroff: F[i+1/2] = F_rusanov + phi (F_LW - F_rusanov).
  - F_LW: two-step (Richtmyer) Lax-Wendroff, U* = 0.5 (U_L + U_R) - 0.5 dt/dx (F(U_R) - F(U_L)), F_LW = F(U*).
  - F_rusanov: local Lax-Friedrichs, a = max(|u| + sqrt(g h cos(theta))) over the interface states.
  - phi in [0, 1]: min over h and hu of the van Leer function of the jump ratios on both sides, capped at 1. phi = 0 where min(h_L, h_R) < 1e-3 m (dry-front safeguard).
  - Two ghost cells per side (the limiter needs the next cell out).
- Pure (unlimited) Lax-Wendroff is not an option: it oscillates at the bore that forms at the slope break and drives h negative. Do not remove the limiter or the dry-front safeguard without rerunning the runout at dx = 1, 0.5 and 0.25 m.
- Explicit update with operator splitting for the source terms.
- Time step: dt = CFL * dx / max(|u| + c). Recompute every step.
- Dry front: threshold h_min below which h and hu are set to zero. Never divide by h without a guard.
- Voellmy drag in the velocity equation is g u|u| / (xi h), which blows up as h -> 0. Treat it semi-implicitly (or skip it below h_min).
- Coulomb stopping rule at u = 0: if |gravity term| <= friction term, the material stays at rest. No sign chatter around u = 0.
- Expect smearing at the debris front (the limiter is first order there) and a convergence rate below 1 near it. Expect a small glitch at the sonic point of a transonic rarefaction (x = 0 in Ritter), a known Lax-Wendroff artifact. State these rather than hiding them.

## Architecture

- `solver/` is a pure NumPy module with no DOM or Pyodide imports. The web app and the verification notebook both import this same module. Never fork the solver for either.
- Pyodide only ships certain packages, so keep dependencies to NumPy. Ask before adding SciPy or anything else.
- The solver runs in a Web Worker and posts snapshots to the main thread. The main thread never runs the solver loop.
- No backend server: everything is static files served by GitHub Pages. Do not add API routes, serverless functions, or anything that needs a running process.
- Deploy with a GitHub Actions workflow that publishes to Pages on push to the main branch. Use relative asset paths, since Pages serves the site under `/<repo-name>/`. Test the deployed URL, not only localhost.
- Load Pyodide and the `solver/` module through paths that work on Pages (a pinned Pyodide version from its CDN, or vendored files). Pin the version.
- The UI shows h(s, t) side view, overhead schematic, final deposit profile, front position vs. time, and summary metrics (runout, peak front speed, peak Froude number, max deposit thickness, mass error).
- Metrics must be computed from solver output. The design mockup's numbers (540 m runout, 18.6 m/s, etc.) are PLACEHOLDERS and must never appear as hard-coded values.

## Code comments (presentation-ready)

This repo may be opened live during the presentation. Anyone reading `solver/` should be able to see exactly how each equation was discretized and why, without a separate explanation. This overrides the usual "keep comments minimal" default. It does not override "simplicity first": comment the code, do not add code to be comment-worthy.

Every solver function gets a docstring covering:
- **Equation:** the continuous equation or term it implements, written in plain-text math (e.g. `(hu)_t + (hu^2 + 0.5 g cos(theta) h^2)_x = S`).
- **Discretization:** the exact formula used, with cell indexing (`U[i]` is the cell average, `F[i+1/2]` is the interface flux between cells i and i+1). Name the scheme (finite volume, flux-limited Lax-Wendroff with a Rusanov / local Lax-Friedrichs fallback).
- **Why:** the reason for each choice, and the alternative rejected. Examples: why Lax-Wendroff is limited (pure Lax-Wendroff goes negative at the bore and the dry front); why Rusanov is the fallback (simple and robust at the dry front); why explicit time stepping with a CFL limit; why the drag term is semi-implicit; why the dry-front threshold exists.
- **Units and shapes:** units of every argument and the array shape and indexing convention.
- **Known weakness:** what the choice costs (e.g. first-order smearing of the front where the limiter switches off, the sonic-point glitch, the pressure-flux jump at the slope break).

Inline comments:
- Put a short comment above each non-obvious line or block saying what it does and why. Prefer explaining intent over restating the code.
- Tie each code block to its equation (`# mass flux: (hu)`, `# hydrostatic pressure term: 0.5 g cos(theta) h^2`).
- Mark every numerical safeguard (h_min guard, u = 0 stopping rule, CFL cap, positivity fix) with a comment stating the failure it prevents.
- Put the update formula next to the loop that applies it: `U[i] -= dt/dx * (F[i+1/2] - F[i-1/2])`.
- Give magic numbers a named constant with a comment (CFL = 0.45, h_min, grid defaults).

Top of each module: a block comment summarizing the whole algorithm in order (initialize, compute dt, fluxes, update, sources, boundaries, repeat) and listing the decisions in "Design decisions" above that affect it.

Keep comments accurate. When code changes, update the comment in the same edit. A comment that contradicts the code is worse than none.

## Verification notebook

A separate Jupyter notebook (`verification/ritter_dam_break.ipynb`) verifies the solver against the exact dam-break solution. It is a deliverable and must be readable in a presentation.

- Import the real solver from `solver/`. Never reimplement the scheme in the notebook.
- Cases: flat bed (theta = 0), no friction (mu = 0 and xi -> infinity). First the dry-bed Ritter solution, then the wet-bed Stoker solution (formulas below). Stoker is required, not optional: Ritter has no shock, so it cannot test the flux limiter (pure, unlimited Lax-Wendroff passes Ritter but fails on the app's runout).
- Overlay numerical h(x, t) and u(x, t) on the exact curves at several times. Mark the rarefaction head at x = -c0 t and the front at x = 2 c0 t.
- Report the L2 (and max) error of h against the exact solution at a fixed time. Run a refinement study: halve dx at fixed CFL and plot error vs. dx on log-log axes with a reference slope-1 line (a guide only, not a pass criterion). State the measured rate and explain any deviation (front smearing lowers the rate).
- Report mass conservation error.
- Second section: steady Voellmy flow against u_t = sqrt(xi h (sin(theta) - mu cos(theta))), to cover the source terms Ritter turns off.
- Markdown cells before each plot say what is being tested, the expected result, and what a failure would mean. End with a short pass/fail summary table with the tolerances used.
- The notebook runs top to bottom in a fresh kernel with no hidden state. Fix random seeds if any are used (there should be none).
- The notebook runs under conda, not Pyodide, so it may use SciPy and Matplotlib for the exact solution and plots. The solver it imports must still depend only on NumPy.

## Verification (required)

1. Ritter dry-bed dam break: theta = 0, no friction. Exact solution: c0 = sqrt(g h0), h = (2 c0 - x/t)^2 / (9 g), u = (2/3)(c0 + x/t), front at x = 2 c0 t. The solution is h0 for x < -c0 t, the parabola above for -c0 t <= x <= 2 c0 t, and 0 beyond (dam at x = 0, reservoir on the left). Report L2 error vs. dx at fixed CFL (log-log). Reference (CFL 0.45, h_min = 1e-10, g = 9.81, h0 = 1, domain [-10, 10], t = 1 s, N = 100 to 3200), measured 2026-11-23:
   - Current scheme (flux-limited Lax-Wendroff): L1 rate rises 0.81 to 0.94, L2 rate 0.72 to 0.90, L2 error 5.3e-3 at N = 3200, mass error 0.
   - Earlier scheme (first-order Rusanov), kept for comparison: L1 0.56 to 0.79, L2 0.54 to 0.72, L2 error 1.2e-2 at N = 3200, mass error <= 4e-16.
   - Both trend toward but stay below 1, because the front and the head kink are handled at first order. Do not make "slope = 1" or "slope = 2" a pass criterion. Report the measured rate and explain it.
2. Stoker wet-bed dam break: theta = 0, no friction, h_L = 1 m, h_R = 0.1 m, t = 1 s. Exact solution: rarefaction with u + 2c constant, so u_m = 2 (c_L - c_m); shock by Rankine-Hugoniot, u_m = (h_m - h_R) sqrt(g (h_m + h_R) / (2 h_m h_R)), speed s = h_m u_m / (h_m - h_R); solve for h_m by root finding (h_m = 0.3962 m, s = 3.1051 m/s). Derived here; SWASHES compiles the same case. Pass criteria: L1 rates > 0.7, shock position within 2 cells, no undershoot ahead of the shock, plateau overshoot < 1% of the jump, mass error < 1e-13. Measured 2026-11-23 (N = 100 to 3200): L1 rates 1.02, 0.95, 0.98, 1.01, 0.83; shock 0.5 to 0.8 cells behind the exact position; overshoot 1.9e-4 m at N = 3200. Do not use L2 as a criterion here: a captured shock limits L2 to rate 1/2 on average, and at N = 1600 to 3200 the L2 error stalls (rate 0.05) because it depends on where the jump falls inside a cell.
3. Steady uniform flow with Voellmy friction: u_t = sqrt(xi h (sin(theta) - mu cos(theta))). Tests the source terms separately, since both dam breaks turn them off.
4. Mass conservation to near machine precision on closed domains.

A change to the flux, source terms, or dry-front handling is not done until all four pass. No exact solution covers the bore at the slope break with friction acting; that is checked only by mass conservation and grid convergence of runout.

## Coding rules

### 1. Think before coding
Don't assume. Don't hide confusion. Surface tradeoffs.
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them. Don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

### 2. Simplicity first
Minimum code that solves the problem. Nothing speculative.
- No features beyond what was asked.
- No abstractions for single-use code.
- No flexibility or configurability that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it. Ask: would a senior engineer say this is overcomplicated? If yes, simplify.

### 3. Surgical changes
Touch only what you must. Clean up only your own mess.
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it. Don't delete it.
- Remove imports, variables, and functions that YOUR changes made unused. Leave pre-existing dead code unless asked.
- Test: every changed line should trace directly to the request.

### 4. Goal-driven execution
Define success criteria. Loop until verified.
- "Add validation" means: write tests for invalid inputs, then make them pass.
- "Fix the bug" means: write a test that reproduces it, then make it pass.
- "Refactor X" means: ensure tests pass before and after.
- For multi-step tasks, state a brief plan, each step with its check:
  1. [Step] -> verify: [check]
  2. [Step] -> verify: [check]
- In this repo the verification suite (Ritter, Stoker, steady Voellmy flow, mass conservation) is the default success criterion for any solver change.

## Decision-making rules (critical thinking partner)

Default mode is constructive disagreement, not agreement. This applies to design decisions, plans, interpretations, and reviews of my work or of your own.

Behavior:
1. Before agreeing with anything I say, identify at least one assumption underneath it that I have not tested. State it plainly.
2. When I propose a decision, idea, plan, or interpretation, first argue the strongest opposing case. Do not soften it. Do not append "but you might be right." Make me defend my position.
3. If I push back, do not retreat because I objected. Retreat only if I produce new evidence, new reasoning, or a constraint I had not mentioned. "Fair point" without new information is not enough.
4. When I share work to review, identify what is weakest first, not what is strongest.
5. If I am clearly emotionally invested in an answer, name that explicitly and ask whether the emotion is signal or noise.
6. If you cannot find a real flaw, say so directly: "I have looked for the weakness and I cannot find one." Do not invent a flaw to perform thoroughness.
7. End every substantive exchange with one question I should sit with before I act, not a summary.

Tone:
- Direct, not aggressive. Specific, not abstract.
- One disagreement at a time, not a list.
- Cite my own words when challenging me.

Do not:
- Open with praise before disagreeing.
- Use "great question," "interesting point," or any opener that reads as flattery.
- Hedge with "I could be wrong but."
- Add a closing reassurance like "your instinct is good."

Project-specific application:
- Any claim in a presentation, README, or UI string must be defensible against the model's limits (see below). Challenge claims that overreach, mine included.
- Treat the open design decisions (Voellmy-only, two-reach bed, flux-limited Lax-Wendroff) as hypotheses to attack, not settled facts, until the verification suite supports them.

## Model limits (keep consistent in docs and UI)

- Width-averaged 1D, so no lateral spreading or overbank flow
- Constant effective friction, with no pore-pressure evolution and no entrainment
- Single-phase mixture
- The scheme drops to first order at the front (flux limiter), so the front is smeared
- Phi and xi trade off against each other, so matching a runout distance does not identify either one
- The phi and xi ranges shown in the UI come from RAMMS, a model with yield-stress and curvature terms this one lacks; they are a guide, not calibrated values for this model

## Commands

All Python runs in the `geophys` conda env (NumPy, SciPy, Matplotlib, Jupyter; no pytest).

- Solver tests: `conda run -n geophys python tests/test_verification.py`
- Verification notebook: `cd verification && conda run -n geophys jupyter nbconvert --to notebook --execute --inplace ritter_dam_break.ipynb`
- Build and serve the site locally (same layout as Pages): `rm -rf _site && mkdir _site && cp -r web/* _site/ && mkdir _site/solver && cp solver/*.py _site/solver/ && conda run -n geophys python -m http.server 8765 -d _site`, then open http://localhost:8765/
- Deploy: push to `main`; `.github/workflows/pages.yml` runs the tests, then publishes `_site` to Pages. One-time setup: repo Settings → Pages → Source = "GitHub Actions".
