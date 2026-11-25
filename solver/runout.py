"""
Debris-flow runout scenario on a two-reach bed, built on the finite-volume core.

This is the configuration the web app runs. It does not add any numerics of its own:
every time step is `core.step`, the same function the verification notebook tests.

Algorithm:
    1. Build the grid along the bed, s in [0, DOMAIN_LENGTH], N cells of width dx.
    2. Build the two-reach bed: slope theta_1 (user) above S_BREAK, RUNOUT_SLOPE_DEG below,
       joined by a linear ramp in theta of width BREAK_WIDTH (see `two_reach_theta`).
    3. Release: block of depth H0 and length L0 = V / (b H0) against the upstream wall,
       at rest (u = 0).
    4. Loop: dt from CFL, `core.step` with a wall upstream and open outflow downstream,
       track peak Froude number and outflow volume; every SNAP_DT record a snapshot.
    5. Stop when the front has stopped advancing (FRONT_STOP_SPEED over STOP_WINDOW) or
       at T_MAX. A thin film may still be draining down the upper reach at that point.
    6. Metrics from the solver output: deposit front position (runout), peak front speed,
       peak Froude number, max deposit thickness, mass-balance error.

Design decisions from CLAUDE.md that affect this module:
    - Width-averaged 1D: channel width b only converts depth to volume, V = sum(h b dx).
      There is no cross-channel dimension (no lateral spreading or overbank flow).
    - Two-reach bed with a fixed 5 degree runout reach starting at s = 350 m. A constant
      slope with phi < theta never produces a deposit, hence the slope break.
    - Release block at rest, wall upstream, open outflow downstream.
    - Material enters only through phi (mu = tan(phi)) and xi.
    - Fixed constants: g, H0 = 3 m, domain 600 m, CFL = 0.45 (in core).

This is a tool for building intuition about what controls runout. It is not a hazard
prediction tool: phi and xi trade off, so matching one runout distance identifies neither.
"""

import numpy as np

from . import core

# Release depth (m), fixed by the project design.
H0 = 3.0
# Length of the computational domain along the bed (m).
DOMAIN_LENGTH = 600.0
# Start of the runout reach along the bed (m) and its fixed slope (degrees).
S_BREAK = 350.0
RUNOUT_SLOPE_DEG = 5.0
# Width (m) of the linear ramp in theta centred on S_BREAK. The sharp break makes
# cos(theta) inside the pressure flux 0.5 g cos(theta) h^2 jump, which turns its
# x-derivative into a grid-dependent spike; a ramp resolved by ~20 cells at the default
# grid makes that term converge under refinement.
BREAK_WIDTH = 20.0
# Default number of cells (dx = 0.5 m on the 600 m domain). Refinement check for
# V = 1e4 m^3, b = 20 m, theta_1 = 30, phi = 20, xi = 500: runout 437 / 433 / 431 / 431 m
# at dx = 2 / 1 / 0.5 / 0.25 m, so dx = 0.5 m is within ~1 m of the converged value.
N_CELLS = 1200
# Depth (m) that defines the flow front for runout and front-speed metrics. The
# scheme drops to first order at the front and leaves a smeared, thin tail ahead of it,
# so the front is where h first drops below this value.
H_FRONT = 0.01
# Run ends when the front has advanced less than FRONT_STOP_SPEED (m/s) on average over
# the last STOP_WINDOW (s). "All velocities zero" is not usable: with phi < theta_1 a thin
# film keeps draining down the upper reach indefinitely (physically correct), and the
# numerical fluxes move mass even where u = 0 (Rusanov dissipation, and the Lax-Wendroff
# predictor, which sees the pressure gradient but not the friction holding the deposit).
# Measured for the default case: after arrest the front creeps ~0.05 m/s at dx = 1 m and
# ~0.03 m/s at dx = 0.5 m, shrinking with dx, so the creep is numerical.
FRONT_STOP_SPEED = 0.25
STOP_WINDOW = 10.0
# Hard limit on simulated time (s) in case the flow never stops.
T_MAX = 300.0
# Simulated time between snapshots (s), also the interval for front-speed differences.
SNAP_DT = 0.5


def two_reach_theta(s, theta1):
    """
    Bed slope theta(s) (rad) for the two-reach bed.

    Equation: theta = theta1 for s <= S_BREAK - BREAK_WIDTH/2, RUNOUT_SLOPE for
        s >= S_BREAK + BREAK_WIDTH/2, linear in s in between.
    Why: see BREAK_WIDTH. Evaluating the same function at cell centres (sources) and at
        interfaces (pressure flux) keeps the two consistent.
    Units and shapes: s (m), any shape; theta1 (rad). Returns theta (rad), same shape.
    Known weakness: the model has no bed-curvature (centripetal) terms, so the ramp only
        changes the slope; it does not model the extra normal stress of a concave bend.
    """
    theta2 = np.radians(RUNOUT_SLOPE_DEG)
    frac = np.clip((s - (S_BREAK - 0.5 * BREAK_WIDTH)) / BREAK_WIDTH, 0.0, 1.0)
    return theta1 + frac * (theta2 - theta1)


def bed_profile(s_faces, theta_faces):
    """
    Bed geometry for plotting: horizontal distance X(s) and elevation Z(s) (m).

    Equation: dX/ds = cos(theta), dZ/ds = -sin(theta), X(0) = Z(0) = 0.
    Discretization: trapezoid rule on the interfaces.
    Units and shapes: s_faces (m), theta_faces (rad), shape (N+1,); returns X, Z (N+1,).
    """
    ds = np.diff(s_faces)
    X = np.concatenate([[0.0], np.cumsum(ds * 0.5 * (np.cos(theta_faces[1:]) + np.cos(theta_faces[:-1])))])
    Z = -np.concatenate([[0.0], np.cumsum(ds * 0.5 * (np.sin(theta_faces[1:]) + np.sin(theta_faces[:-1])))])
    return X, Z


def front_position(h, s, dx):
    """
    Front position s_front (m): last point where h crosses H_FRONT, linearly interpolated.

    Discretization: j = last cell with h >= H_FRONT; s_front lies between centres s[j]
        and s[j+1] where the linear interpolant of h equals H_FRONT. Interpolation keeps
        s_front from jumping in steps of dx, which would make front speeds noisy.
    Units and shapes: h (m), s (m) cell centres, shape (N,); dx (m). Returns s_front (m);
        0 if no cell exceeds H_FRONT, DOMAIN_LENGTH if the front reached the outlet.
    """
    above = np.nonzero(h >= H_FRONT)[0]
    if above.size == 0:
        return 0.0
    j = above[-1]
    if j == h.size - 1:
        return DOMAIN_LENGTH
    return s[j] + (h[j] - H_FRONT) / (h[j] - h[j + 1]) * dx


def run_runout(volume, width, theta1_deg, phi_deg, xi, n_cells=N_CELLS, on_snapshot=None):
    """
    Simulate release, flow and deposition on the two-reach bed and compute metrics.

    Equation: see core (U_t + F(U)_x = S(U) with Coulomb/Voellmy friction).
    Discretization: core.step with bc = ("wall", "open"), dt from core.cfl_dt.
    Why the metrics are defined as they are:
        - runout = deposit front position s_front (m along the bed from the upstream wall)
          at the end of the run, with the front defined by H_FRONT.
        - peak front speed = max over snapshots of (s_front(t) - s_front(t - SNAP_DT)) / SNAP_DT.
        - peak Froude number = max over all steps of |u| / sqrt(g h cos(theta)) at the cell
          with the largest discharge |hu| (the flow body; at the front itself Fr -> inf).
        - mass error = (V_in_domain + V_outflow - V0) / V0, should be ~1e-15.
    Units and shapes: volume (m^3), width b (m), theta1_deg, phi_deg (degrees), xi (m/s^2).
        on_snapshot(t, h, u, s_front) is called every SNAP_DT (t in s, h in m, u in m/s,
        arrays of shape (N,)); the web worker uses it to stream results.
    Returns a dict with grid, bed geometry, front history and metrics.
    Known weakness: smearing at the front (where the limiter falls back to first-order
        Rusanov) moves the H_FRONT contour ahead of the true front,
        so runout is slightly overestimated and grid-dependent (refine n_cells to check).
        Numerical diffusion also makes an arrested deposit creep, so the stopping rule
        itself shifts runout by a few metres. If the front reaches the outlet, runout is
        reported as DOMAIN_LENGTH and only means "at least 600 m".
    """
    dx = DOMAIN_LENGTH / n_cells
    # cell interfaces s_{i-1/2} (N+1,) and centres s_i (N,)
    s_faces = np.linspace(0.0, DOMAIN_LENGTH, n_cells + 1)
    s = 0.5 * (s_faces[:-1] + s_faces[1:])

    theta1 = np.radians(theta1_deg)
    theta_c = two_reach_theta(s, theta1)
    theta_f = two_reach_theta(s_faces, theta1)
    bed = {"cos_c": np.cos(theta_c), "sin_c": np.sin(theta_c), "cos_f": np.cos(theta_f)}
    mu = np.tan(np.radians(phi_deg))

    # Release block: depth H0 over length L0 = V / (b H0). The last partly covered cell
    # gets a proportional depth so the discrete volume equals V exactly.
    L0 = volume / (width * H0)
    if L0 > S_BREAK - 0.5 * BREAK_WIDTH:
        raise ValueError(
            f"Release block is {L0:.0f} m long (V / (b h0)), longer than the upper reach "
            f"({S_BREAK - 0.5 * BREAK_WIDTH:.0f} m). Reduce the volume or widen the channel.")
    h = H0 * np.clip((L0 - s_faces[:-1]) / dx, 0.0, 1.0)
    hu = np.zeros(n_cells)
    V0 = np.sum(h) * dx * width

    t = 0.0
    v_out = 0.0
    peak_froude = 0.0
    times = [0.0]
    fronts = [front_position(h, s, dx)]
    if on_snapshot is not None:
        on_snapshot(0.0, h, np.zeros(n_cells), fronts[0])
    next_snap = SNAP_DT
    stopped = False

    while True:
        dt = core.cfl_dt(h, hu, bed["cos_c"], dx)
        # wall upstream (closed release area), open outflow at the downstream end
        h, hu, _, flux_right = core.step(h, hu, dt, dx, bed, mu, xi, ("wall", "open"))
        t += dt
        # volume that left through the outlet (m^3)
        v_out += dt * flux_right * width

        u = core.velocity(h, hu)
        # Froude number of the flow body: Fr = |u| / sqrt(g h cos(theta)) at the cell with
        # the largest discharge |hu|. The max over all cells is useless here: at a dry-bed
        # front h -> 0 with finite u, so Fr -> infinity and the value only reflects the grid.
        i = np.argmax(np.abs(hu))
        if h[i] >= H_FRONT:
            fr = abs(u[i]) / np.sqrt(core.G * h[i] * bed["cos_c"][i])
            peak_froude = max(peak_froude, fr)

        if t >= next_snap:
            times.append(t)
            fronts.append(front_position(h, s, dx))
            if on_snapshot is not None:
                on_snapshot(t, h, u, fronts[-1])
            next_snap += SNAP_DT
            # Stopping rule: front (nearly) stationary over the last STOP_WINDOW. Running
            # on would only let numerical diffusion creep the deposit (see FRONT_STOP_SPEED).
            k = round(STOP_WINDOW / SNAP_DT)
            if len(fronts) > k:
                stopped = (fronts[-1] - fronts[-1 - k]) / (times[-1] - times[-1 - k]) < FRONT_STOP_SPEED
            if stopped or t >= T_MAX:
                break

    times = np.array(times)
    fronts = np.array(fronts)
    front_speeds = np.diff(fronts) / np.diff(times)
    X, Z = bed_profile(s_faces, theta_f)

    return {
        "s": s,
        "bed_x": X,
        "bed_z": Z,
        "s_faces": s_faces,
        "times": times,
        "fronts": fronts,
        "release_length": L0,
        "t_final": t,
        "stopped": stopped,
        "runout": float(fronts[-1]),
        "peak_front_speed": float(np.max(front_speeds)),
        "peak_froude": float(peak_froude),
        "max_deposit_thickness": float(np.max(h)),
        "outflow_fraction": v_out / V0,
        "mass_error": (np.sum(h) * dx * width + v_out - V0) / V0,
    }
