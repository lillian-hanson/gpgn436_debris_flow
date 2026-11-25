"""
Finite-volume core for a 1D shallow-water-form flow with Coulomb/Voellmy friction.

System solved (one coupled hyperbolic system, U_t + F(U)_x = S(U), U = [h, hu]):

    mass:      h_t  + (hu)_x                                 = 0
    momentum: (hu)_t + (hu^2 + 0.5 g cos(theta) h^2)_x       = g h sin(theta) - tau_b/rho
    friction:  tau_b/rho = mu g h cos(theta) sgn(u) + g u|u| / xi,   mu = tan(phi)

x is distance along the bed (m), h is flow depth normal to the bed (m), u = hu/h is the
depth-averaged bed-parallel velocity (m/s, a solved variable). Coulomb friction is the
xi -> infinity limit of the same code path (xi = np.inf is allowed), not a separate branch.
The pressure term uses an isotropic earth-pressure coefficient K = 1, so this is a
"shallow-water-form model with Coulomb/Voellmy friction", not the full Savage-Hutter model.

Algorithm, one call to `step` (Lie / Godunov operator splitting):
    1. dt = CFL * dx / max(|u| + sqrt(g h cos(theta)))       (`cfl_dt`, recomputed every step)
    2. add two ghost cells per side from the boundary condition (`add_ghosts`)
    3. at every interface i+1/2, two candidate fluxes:
         high order: two-step (Richtmyer) Lax-Wendroff flux     (`lax_wendroff_flux`)
         low order:  Rusanov (local Lax-Friedrichs) flux        (`rusanov_flux`)
    4. flux limiter phi in [0, 1] from the smoothness of h and hu (`flux_limiter`);
       F = F_rusanov + phi (F_LW - F_rusanov), with phi = 0 next to dry cells
    5. conservative update U[i] -= dt/dx (F[i+1/2] - F[i-1/2]) (`hyperbolic_update`)
    6. sources over the same dt: gravity, then Coulomb (with the u = 0 stopping rule),
       then semi-implicit Voellmy drag                         (`apply_sources`)
    7. dry-front guard: hu = 0 where h < H_MIN
    repeat.

Design decisions from CLAUDE.md that affect this module:
    - Finite volume, flux-limited Lax-Wendroff (Lax-Wendroff where the solution is smooth,
      Rusanov at shocks and the dry front), explicit in time.
    - Source terms by operator splitting; Voellmy drag semi-implicit; Coulomb stopping rule.
    - Dry-front threshold H_MIN; never divide by h without that guard.
    - Bed slope theta can vary along x (two-reach bed). cos(theta) inside the pressure flux
      is evaluated at the interface and used for BOTH states there (see `rusanov_flux`).
    - NumPy only: this exact module runs in the browser (Pyodide) and in the notebook.

Array conventions used throughout:
    N cells, cell i covers [x_{i-1/2}, x_{i+1/2}], width dx (uniform).
    h, hu, cos_c, sin_c : shape (N,), cell averages / cell-centre values, U[i].
    cos_f               : shape (N+1,), interface values; cos_f[i] sits at x_{i-1/2},
                          so the flux between cells i and i+1 (F[i+1/2]) is index i+1.
    Fluxes returned     : shape (N+1,), same interface indexing as cos_f.

References: Savage & Hutter (1989) J. Fluid Mech. 199; Voellmy (1955) Schweiz. Bauzeitung 73;
George & Iverson (2014) Proc. R. Soc. A 470, 20130820 (hyperbolicity); Lax & Wendroff (1960)
Commun. Pure Appl. Math. 13; Richtmyer (1963) two-step variant; van Leer (1974) J. Comput.
Phys. 14 (limiter); Rusanov (1961) USSR Comput. Math. Math. Phys. 1(2); Toro, "Riemann
Solvers and Numerical Methods for Fluid Dynamics", and LeVeque (2002), "Finite Volume
Methods for Hyperbolic Problems" (ch. 6, flux-limited high-resolution methods).
"""

import numpy as np

# Gravitational acceleration (m/s^2).
G = 9.81

# Courant number. Lax-Wendroff is stable for CFL <= 1, but the Rusanov flux it falls back
# to keeps h >= 0 only for CFL <= 0.5; 0.45 leaves a margin because dt is computed from
# cell-centre states while the flux uses interface states.
CFL = 0.45

# Dry-cell threshold (m). Below this depth the cell is treated as dry: its velocity is
# defined as 0 and hu is reset to 0. Prevents u = hu/h from blowing up at the front, where
# numerical diffusion leaves an exponentially thin tail of h. h itself is NOT zeroed,
# because deleting that tail would destroy exact mass conservation.
H_MIN = 1e-10

# Depth (m) below which an interface uses the Rusanov flux only (limiter phi = 0). Pure
# Lax-Wendroff overshoots in a thin layer and drives h negative at the dry front (tested:
# h reached -27 m and the run blew up at dx = 0.5 m without this and the limiter).
H_DRY_FIRST_ORDER = 1e-3

# Boundary conditions, passed as bc = (left, right), one string per side:
#   "wall"     : reflective, ghosts mirror the edge cells with hu negated; zero mass flux.
#   "open"     : transmissive / zero-gradient, ghosts = copies of the edge cell.
#   "periodic" : ghosts = cells on the opposite end (must be used on both sides).


def velocity(h, hu):
    """
    Depth-averaged velocity u = hu / h with a dry-cell guard.

    Equation: u = hu / h.
    Discretization: pointwise per cell; u[i] = hu[i] / h[i] if h[i] >= H_MIN else 0.
    Why: the solver stores conserved variables (h, hu); u is needed for the flux and the
        friction law. Dividing by a near-zero h at the dry front would produce huge or
        infinite u and destroy the CFL time step, so dry cells get u = 0.
    Units and shapes: h (m), hu (m^2/s), returns u (m/s); any matching shapes.
    Known weakness: the threshold is arbitrary; the thin tail below H_MIN is frozen
        (carried only by numerical diffusion, not advection).
    """
    # Dry-front guard: never divide by h < H_MIN (prevents u -> infinity at the front).
    wet = h >= H_MIN
    return np.where(wet, hu / np.where(wet, h, 1.0), 0.0)


def cfl_dt(h, hu, cos_c, dx):
    """
    Explicit time step from the CFL condition.

    Equation: wave speeds of the flux Jacobian, lambda = u +/- sqrt(g h cos(theta)).
    Discretization: dt = CFL * dx / max_i(|u[i]| + sqrt(g h[i] cos_c[i])).
    Why: an explicit scheme is stable only if no wave crosses more than a fraction of a
        cell per step. Explicit stepping was chosen over implicit because the dynamics we
        want (a moving front) has the same time scale as the waves anyway, and both the
        Lax-Wendroff and Rusanov fluxes are stable under this limit (see CFL).
    Units and shapes: h (m), hu (m^2/s), cos_c (-) shape (N,); dx (m); returns dt (s).
    Known weakness: recomputed every step, so dt shrinks when the flow is fast.
    """
    u = velocity(h, hu)
    # fastest signal speed |u| + c, with c = sqrt(g h cos(theta)) the gravity-wave speed
    smax = np.max(np.abs(u) + np.sqrt(G * h * cos_c))
    return CFL * dx / smax


def add_ghosts(h, hu, bc):
    """
    Pad the state with two ghost cells on each side to impose boundary conditions.

    Equation: none (boundary data for the fluxes at x_{1/2} and x_{N+1/2}).
    Discretization: returns arrays of shape (N+4,); indices 0, 1 and N+2, N+3 are ghosts,
        so real cell i sits at padded index i+2.
        wall:     mirror image, ghost k cells outside = edge cell k cells inside, hu negated
                  (the mass flux through the wall is then exactly 0 for both fluxes)
        open:     both ghosts = copy of the edge cell (zero gradient, lets waves leave)
        periodic: ghosts = the cells on the opposite end
    Why two: the flux at an interface needs one cell on each side; the limiter also needs
        the next cell out on each side to measure smoothness (see `flux_limiter`).
    Units and shapes: h (m), hu (m^2/s), shape (N,); bc = (left, right) strings.
    Known weakness: the "open" condition is only approximately non-reflecting.
    """
    left, right = bc

    # left ghosts, outermost first
    if left == "periodic":
        hl, hul = h[-2:], hu[-2:]
    elif left == "wall":
        hl, hul = h[1::-1], -hu[1::-1]
    else:
        hl, hul = np.full(2, h[0]), np.full(2, hu[0])

    # right ghosts, innermost first
    if right == "periodic":
        hr, hur = h[:2], hu[:2]
    elif right == "wall":
        hr, hur = h[:-3:-1], -hu[:-3:-1]
    else:
        hr, hur = np.full(2, h[-1]), np.full(2, hu[-1])

    return np.concatenate([hl, h, hr]), np.concatenate([hul, hu, hur])


def physical_flux(h, hu, cos_theta):
    """
    Physical flux F(U) of the system.

    Equation: F(U) = [hu, hu^2/h + 0.5 g cos(theta) h^2].
    Units and shapes: h (m), hu (m^2/s), cos_theta (-), matching shapes; returns
        (mass flux in m^2/s, momentum flux in m^3/s^2).
    """
    u = velocity(h, hu)
    # mass flux: (hu); momentum flux: hu * u (advection) + 0.5 g cos(theta) h^2 (pressure)
    return hu, hu * u + 0.5 * G * cos_theta * h**2


def rusanov_flux(hL, huL, hR, huR, cos_f):
    """
    Rusanov (local Lax-Friedrichs) numerical flux at every interface.

    Equation: physical flux F(U) = [hu, hu^2/h + 0.5 g cos(theta) h^2].
    Discretization (first-order finite volume, interface i+1/2 between cells i and i+1):
        a       = max(|u_L| + sqrt(g h_L cos_f), |u_R| + sqrt(g h_R cos_f))
        F[i+1/2] = 0.5 (F(U_L) + F(U_R)) - 0.5 a (U_R - U_L)
    Role: the low-order fallback of the flux-limited scheme, used where the limiter
        detects a shock or the dry front (see `hyperbolic_update`).
    Why Rusanov as the fallback: it needs only the fastest wave speed, not a full Riemann
        solve, and it is robust and positivity preserving at a dry front (h -> 0), where the
        system loses strict hyperbolicity. The rejected alternatives (Roe, HLL/HLLC, exact
        Godunov) are less diffusive but more complex and not part of the agreed scope.
    Why one cos_f for both states: with theta varying along x, using cos(theta_i) on the
        left state and cos(theta_{i+1}) on the right would make the flux inconsistent
        (F(U, U) != physical flux) at the slope break. Evaluating cos(theta) at the
        interface and using it on both sides keeps F consistent and conservative.
    Units and shapes: h (m), hu (m^2/s), cos_f (-), all shape (N+1,); returns
        (F_h in m^2/s, F_hu in m^3/s^2), shape (N+1,).
    Known weakness: the -0.5 a (U_R - U_L) term is numerical diffusion of size ~0.5 a dx;
        wherever the limiter selects this flux it smears the solution (the front, shocks).
    """
    uL = velocity(hL, huL)
    uR = velocity(hR, huR)

    # local maximum wave speed |u| + sqrt(g h cos(theta)) over both interface states
    a = np.maximum(np.abs(uL) + np.sqrt(G * hL * cos_f),
                   np.abs(uR) + np.sqrt(G * hR * cos_f))

    fhL, fhuL = physical_flux(hL, huL, cos_f)
    fhR, fhuR = physical_flux(hR, huR, cos_f)

    # central average minus the Rusanov dissipation term
    F_h = 0.5 * (fhL + fhR) - 0.5 * a * (hR - hL)
    F_hu = 0.5 * (fhuL + fhuR) - 0.5 * a * (huR - huL)
    return F_h, F_hu


def lax_wendroff_flux(hL, huL, hR, huR, cos_f, dt, dx):
    """
    Two-step (Richtmyer) Lax-Wendroff flux at every interface.

    Equation: U_t + F(U)_x = 0 with F as in `physical_flux`.
    Discretization (interface i+1/2 between cells i and i+1):
        predictor: U*[i+1/2] = 0.5 (U_L + U_R) - 0.5 dt/dx (F(U_R) - F(U_L))   (half step)
        flux:      F[i+1/2]  = F(U*[i+1/2])
    Why Lax-Wendroff: second order in space and time where the solution is smooth, so the
        rarefaction and the body of the flow are much less smeared than with Rusanov. The
        two-step Richtmyer form is the conservative version for a nonlinear system: it never
        needs the flux Jacobian, unlike the one-step form U - dt A U_x + dt^2/2 (A^2 U_x)_x.
    Units and shapes: as `rusanov_flux`, plus dt (s) and dx (m).
    Known weakness: dispersive. On its own it oscillates at shocks and goes negative at a
        dry front, so it is only used through the limiter in `hyperbolic_update`.
    """
    fhL, fhuL = physical_flux(hL, huL, cos_f)
    fhR, fhuR = physical_flux(hR, huR, cos_f)

    # predictor: half-step state at the interface
    h_half = 0.5 * (hL + hR) - 0.5 * dt / dx * (fhR - fhL)
    hu_half = 0.5 * (huL + huR) - 0.5 * dt / dx * (fhuR - fhuL)
    # Positivity safeguard: a negative predicted depth (possible at a strong shock) would
    # give a negative pressure term; clip it to 0.
    h_half = np.maximum(h_half, 0.0)

    # corrector flux: the physical flux of the half-step state
    return physical_flux(h_half, hu_half, cos_f)


def flux_limiter(hg, hug, hL, hR):
    """
    Limiter phi[i+1/2] in [0, 1] that blends Rusanov (phi = 0) and Lax-Wendroff (phi = 1).

    Discretization, for q = h and q = hu, with jumps D[i+1/2] = q[i+1] - q[i]:
        r- = D[i-1/2] / D[i+1/2],   r+ = D[i+3/2] / D[i+1/2]   (smoothness on each side)
        phi_q = min(1, psi(r-), psi(r+)),   psi(r) = (r + |r|) / (1 + |r|)   (van Leer)
        phi   = min(phi_h, phi_hu),  and phi = 0 if min(h_L, h_R) < H_DRY_FIRST_ORDER.
    Why: where the jumps on both sides are about the same (r ~ 1, smooth flow) psi ~ 1 and
        the second-order Lax-Wendroff flux is used. At an extremum (r <= 0, e.g. the bore at
        the slope break or the top of the deposit) psi = 0 and the scheme drops to Rusanov,
        which is what prevents Lax-Wendroff's oscillations. Both sides are checked so that
        no upwind direction has to be chosen for a system with waves going both ways.
        Capping phi at 1 keeps the flux a convex blend of the two.
    Units and shapes: hg, hug ghost-padded state, shape (N+4,); hL, hR interface states,
        shape (N+1,). Returns phi, shape (N+1,).
    Known weakness: this symmetric limiter is a pragmatic choice, not a proven TVD scheme
        for the system; the verification suite (positivity, mass, convergence) is the check.
        It drops to first order at smooth extrema too, which caps the observed accuracy.
    """
    phi = np.ones(hL.size)
    for q in (hg, hug):
        # D[k] = q[k+1] - q[k]; interface i+1/2 sits between padded cells i+2 and i+3
        D = np.diff(q)
        d0, dm, dp = D[1:-1], D[:-2], D[2:]
        flat = d0 == 0.0
        d0_safe = np.where(flat, 1.0, d0)
        r_minus = dm / d0_safe
        r_plus = dp / d0_safe
        psi = np.minimum((r_minus + np.abs(r_minus)) / (1.0 + np.abs(r_minus)),
                         (r_plus + np.abs(r_plus)) / (1.0 + np.abs(r_plus)))
        # a zero jump means both fluxes agree there; leave phi unconstrained
        phi = np.minimum(phi, np.where(flat, 1.0, psi))

    # Dry-front safeguard: first order next to (nearly) dry cells, where Lax-Wendroff
    # overshoots and would drive h negative.
    return np.where(np.minimum(hL, hR) < H_DRY_FIRST_ORDER, 0.0, phi)


def hyperbolic_update(h, hu, dt, dx, cos_f, bc):
    """
    One explicit conservative step of the homogeneous system U_t + F(U)_x = 0.

    Equation: h_t + (hu)_x = 0,  (hu)_t + (hu^2 + 0.5 g cos(theta) h^2)_x = 0.
    Discretization (finite volume, flux-limited Lax-Wendroff):
        F[i+1/2] = F_rusanov + phi[i+1/2] (F_LW - F_rusanov)
        U[i] <- U[i] - dt/dx * (F[i+1/2] - F[i-1/2])
    Why: Lax-Wendroff alone is second order but oscillates at shocks and goes negative at
        the dry front; Rusanov alone is robust but first order everywhere. The limiter uses
        each where it works (LeVeque 2002, ch. 6). The conservative form makes the total of
        h*dx change only through the two boundary fluxes, so mass is conserved to round-off.
    Units and shapes: h (m), hu (m^2/s) shape (N,); dt (s); dx (m); cos_f shape (N+1,).
        Returns new (h, hu) and the mass fluxes through the left and right boundaries
        (m^2/s, positive = in the +x direction) for bookkeeping of outflow.
    Known weakness: first order wherever phi = 0 (shocks, extrema, the dry front), so the
        front is still smeared and the convergence rate stays below 2 (and below 1 in L2
        near the front).
    """
    hg, hug = add_ghosts(h, hu, bc)
    # interface j (j = 0..N) lies between padded cells j+1 (left) and j+2 (right)
    hL, hR = hg[1:-2], hg[2:-1]
    huL, huR = hug[1:-2], hug[2:-1]

    F_lo_h, F_lo_hu = rusanov_flux(hL, huL, hR, huR, cos_f)
    F_hi_h, F_hi_hu = lax_wendroff_flux(hL, huL, hR, huR, cos_f, dt, dx)
    phi = flux_limiter(hg, hug, hL, hR)

    # limited flux: Rusanov plus a limited fraction of the Lax-Wendroff correction
    F_h = F_lo_h + phi * (F_hi_h - F_lo_h)
    F_hu = F_lo_hu + phi * (F_hi_hu - F_lo_hu)

    # conservative update: U[i] -= dt/dx * (F[i+1/2] - F[i-1/2])
    h_new = h - dt / dx * (F_h[1:] - F_h[:-1])
    hu_new = hu - dt / dx * (F_hu[1:] - F_hu[:-1])

    # Positivity safeguard: removes round-off negatives so sqrt(g h cos) never sees h < 0.
    # The verification suite checks that this clip never removes measurable mass.
    h_new = np.maximum(h_new, 0.0)
    return h_new, hu_new, F_h[0], F_h[-1]


def apply_sources(h, hu, dt, sin_c, cos_c, mu, xi):
    """
    Source-term substep (operator splitting): gravity, Coulomb friction, Voellmy drag.

    Equation (velocity form, h is unchanged by the sources):
        u_t = g sin(theta) - mu g cos(theta) sgn(u) - g u|u| / (xi h)
    Discretization, per cell, in this order:
        1. gravity (explicit):        u* = u + dt g sin(theta)
        2. Coulomb (projection):      u** = sgn(u*) max(|u*| - dt mu g cos(theta), 0)
        3. Voellmy drag (implicit):   u + k u|u| = u**,   k = dt g / (xi h)
              solved exactly:         |u| = 2|u**| / (1 + sqrt(1 + 4 k |u**|)),  sgn(u) = sgn(u**)
        then hu = h u.
    Why:
        - Coulomb as a projection is the stopping rule: friction can decelerate u to exactly
          0 but never reverse it. A cell at rest (u = 0) whose driving acceleration
          (gravity + pressure gradient, already in u from the flux step) is at most the
          friction stays exactly at rest. An explicit -mu g cos sgn(u) term would instead
          make u chatter around 0 every step.
        - The drag coefficient g/(xi h) blows up as h -> 0, so an explicit update would
          force tiny time steps or overshoot through u = 0. Backward Euler in u is
          unconditionally stable, never flips the sign, and its quadratic has the closed
          form above (written in the rationalized form so that xi = inf, k = 0 is exact:
          Coulomb is just the xi -> infinity limit, no separate code path).
        - The fully implicit quadratic (not a lagged linearization) keeps the steady
          Voellmy state u_t = sqrt(xi h (sin - mu cos)) an exact fixed point of the step.
    Units and shapes: h (m), hu (m^2/s), sin_c, cos_c (-) shape (N,); dt (s);
        mu = tan(phi) (-); xi (m/s^2), may be np.inf. Returns new hu (m^2/s).
    Known weakness: splitting is first order in time, so transients (not steady states)
        carry an O(dt) splitting error.
    """
    # Dry-front guard: sources act only on wet cells; dry cells keep hu = 0.
    wet = h >= H_MIN
    u = velocity(h, hu)

    # 1. gravity component along the bed: + g sin(theta)
    u = u + dt * G * sin_c

    # 2. Coulomb friction mu g cos(theta), as a projection toward u = 0.
    # Stopping rule: if |u| <= dt mu g cos(theta) the cell stops (or stays) at rest exactly.
    u = np.sign(u) * np.maximum(np.abs(u) - dt * mu * G * cos_c, 0.0)

    # 3. Voellmy drag g u|u| / (xi h), backward Euler. h_safe avoids 1/0 in dry cells
    # (their result is discarded below).
    h_safe = np.where(wet, h, 1.0)
    k = dt * G / (xi * h_safe)
    u = np.sign(u) * 2.0 * np.abs(u) / (1.0 + np.sqrt(1.0 + 4.0 * k * np.abs(u)))

    # back to the conserved variable; dry cells carry no momentum
    return np.where(wet, h * u, 0.0)


def step(h, hu, dt, dx, bed, mu, xi, bc):
    """
    Advance (h, hu) by one time step dt: hyperbolic update, then sources.

    Equation: U_t + F(U)_x = S(U), U = [h, hu] (see module docstring).
    Discretization: Lie splitting, U^{n+1} = Sources_dt( Hyperbolic_dt( U^n ) ).
    Why: splitting lets the stiff friction terms be treated (semi-)implicitly while the
        flux stays explicit. Lie splitting is first order in time; Strang splitting would
        be second order, but the limiter already drops the scheme to first order at the
        front and shocks, which dominate the error, so it was not worth the extra cost.
    Units and shapes: see `hyperbolic_update` and `apply_sources`; bed is a dict with
        "cos_c", "sin_c" (shape (N,)) and "cos_f" (shape (N+1,)). Returns
        (h, hu, left boundary mass flux, right boundary mass flux).
    Known weakness: first order in time overall (splitting), and first order in space
        wherever the limiter falls back to Rusanov.
    """
    h, hu, flux_left, flux_right = hyperbolic_update(h, hu, dt, dx, bed["cos_f"], bc)
    hu = apply_sources(h, hu, dt, bed["sin_c"], bed["cos_c"], mu, xi)
    return h, hu, flux_left, flux_right


def run_until(h, hu, t, t_end, dx, bed, mu, xi, bc):
    """
    Integrate from time t to exactly t_end.

    Discretization: repeat dt = min(cfl_dt, t_end - t); step(...).
    Why: verification compares against exact solutions at a fixed time, so the last step
        is shortened to land on t_end exactly (shortening never violates the CFL limit).
    Units and shapes: as in `step`; t, t_end (s). Returns (h, hu, t_end, mass that left
        through the left boundary, mass that left through the right boundary), the
        masses being per unit width (m^2).
    """
    out_left = 0.0
    out_right = 0.0
    while t < t_end:
        # CFL cap, shortened on the final step to hit t_end exactly
        dt = min(cfl_dt(h, hu, bed["cos_c"], dx), t_end - t)
        h, hu, fl, fr = step(h, hu, dt, dx, bed, mu, xi, bc)
        # outflow bookkeeping: positive flux at the left boundary is INTO the domain
        out_left -= dt * fl
        out_right += dt * fr
        t += dt
    return h, hu, t_end, out_left, out_right


def uniform_bed(n_cells, theta):
    """
    Bed geometry for a constant slope theta (rad): cell and interface trig values.

    Units and shapes: returns dict with "cos_c", "sin_c" (N,) and "cos_f" (N+1,).
    """
    return {
        "cos_c": np.full(n_cells, np.cos(theta)),
        "sin_c": np.full(n_cells, np.sin(theta)),
        "cos_f": np.full(n_cells + 1, np.cos(theta)),
    }
