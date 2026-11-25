"""
Verification suite for the solver (the default success criterion for any solver change).

    1. Ritter dry-bed dam break (flat bed, no friction): error decreases under refinement.
    2. Stoker wet-bed dam break (flat bed, no friction): a real shock, so it tests the flux
       limiter (shock position, no oscillations, L1 convergence).
    3. Voellmy flow on a uniform slope: transient and terminal velocity vs. the exact ODE solution.
    4. Mass conservation on a closed domain, and mass bookkeeping in the app's runout run.

Runs with plain Python (no pytest needed):  python tests/test_verification.py
Every test calls the same solver functions the web app uses.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from solver import core, runout  # noqa: E402


def ritter_exact(x, t, h0):
    """Ritter dry-bed dam break (dam at x = 0, reservoir x < 0): returns (h, u)."""
    c0 = np.sqrt(core.G * h0)
    xi = x / t
    h = np.where(xi < -c0, h0, np.where(xi > 2 * c0, 0.0, (2 * c0 - xi) ** 2 / (9 * core.G)))
    u = np.where((xi >= -c0) & (xi <= 2 * c0), 2.0 / 3.0 * (c0 + xi), 0.0)
    return h, u


def run_ritter(n_cells, h0=1.0, t_end=1.0, x_min=-10.0, x_max=10.0):
    """Solve the dam break with the real solver; returns x, h, hu, dx, mass before/after."""
    dx = (x_max - x_min) / n_cells
    x = x_min + (np.arange(n_cells) + 0.5) * dx
    h = np.where(x < 0.0, h0, 0.0)
    hu = np.zeros(n_cells)
    m0 = np.sum(h) * dx
    bed = core.uniform_bed(n_cells, 0.0)
    h, hu, _, out_l, out_r = core.run_until(h, hu, 0.0, t_end, dx, bed, 0.0, np.inf, ("open", "open"))
    return x, h, hu, dx, m0, np.sum(h) * dx + out_l + out_r


def test_ritter():
    errors = []
    sizes = [100, 200, 400, 800, 1600]
    for n in sizes:
        x, h, hu, dx, m0, m1 = run_ritter(n)
        h_ex, _ = ritter_exact(x, 1.0, 1.0)
        errors.append(np.sqrt(np.sum((h - h_ex) ** 2) * dx))
        assert abs(m1 - m0) / m0 < 1e-13, f"Ritter mass error {abs(m1 - m0) / m0:.2e}"
    errors = np.array(errors)
    rates = np.log2(errors[:-1] / errors[1:])
    print("  Ritter L2 errors:", np.array2string(errors, precision=3))
    print("  Ritter L2 rates: ", np.array2string(rates, precision=2))
    # Pass: error decreases monotonically with a positive observed rate. Not "rate = 1" or
    # "rate = 2": the limiter falls back to first order at the dry front and the
    # rarefaction-head kink, which dominate the error.
    assert np.all(rates > 0.4), "Ritter error is not converging"


def stoker_middle_state(hL, hR):
    """
    Middle state (h_m, u_m) and shock speed s of the wet-bed dam break (Stoker).

    The left-going rarefaction keeps u + 2c constant:  u_m = 2 (c_L - c_m).
    The right-going shock satisfies Rankine-Hugoniot:  u_m = (h_m - h_R) sqrt(g (h_m + h_R) / (2 h_m h_R)),
    with speed s = h_m u_m / (h_m - h_R). h_m is found by bisection on (h_R, h_L), so the
    test needs only NumPy (CI does not install SciPy).
    """
    cL = np.sqrt(core.G * hL)

    def mismatch(hm):
        return 2 * (cL - np.sqrt(core.G * hm)) - (hm - hR) * np.sqrt(core.G * (hm + hR) / (2 * hm * hR))

    lo, hi = hR, hL                      # mismatch(hR) > 0 > mismatch(hL)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if mismatch(mid) > 0 else (lo, mid)
    hm = 0.5 * (lo + hi)
    um = 2 * (cL - np.sqrt(core.G * hm))
    return hm, um, hm * um / (hm - hR)


def stoker_exact(x, t, hL, hR):
    """Wet-bed dam break: rarefaction, constant middle state, shock into depth hR."""
    hm, um, s = stoker_middle_state(hL, hR)
    cL, cm = np.sqrt(core.G * hL), np.sqrt(core.G * hm)
    xi = x / t
    h = np.where(xi < -cL, hL, np.where(xi < um - cm, (2 * cL - xi) ** 2 / (9 * core.G),
                                        np.where(xi < s, hm, hR)))
    u = np.where(xi < -cL, 0.0, np.where(xi < um - cm, 2.0 / 3.0 * (cL + xi),
                                         np.where(xi < s, um, 0.0)))
    return h, u


def test_stoker():
    hL, hR, t_end = 1.0, 0.1, 1.0
    hm, um, s = stoker_middle_state(hL, hR)
    # the exact solution must satisfy both Rankine-Hugoniot conditions at the shock
    assert abs(s * (hm - hR) - hm * um) < 1e-12
    assert abs(s * hm * um - (hm * um**2 + 0.5 * core.G * (hm**2 - hR**2))) < 1e-12

    L1 = []
    for n in [100, 200, 400, 800, 1600]:
        dx = 20.0 / n
        x = -10.0 + (np.arange(n) + 0.5) * dx
        h, hu = np.where(x < 0.0, hL, hR), np.zeros(n)
        m0 = np.sum(h) * dx
        h, hu, _, out_l, out_r = core.run_until(h, hu, 0.0, t_end, dx, core.uniform_bed(n, 0.0),
                                                0.0, np.inf, ("open", "open"))
        h_ex, _ = stoker_exact(x, t_end, hL, hR)
        L1.append(np.sum(np.abs(h - h_ex)) * dx)
        assert abs(np.sum(h) * dx + out_l + out_r - m0) / m0 < 1e-13

        # shock position: last crossing of the mid-jump depth, interpolated
        h_mid = 0.5 * (hm + hR)
        j = np.nonzero(h >= h_mid)[0][-1]
        x_shock = x[j] + (h[j] - h_mid) / (h[j] - h[j + 1]) * dx
        assert abs(x_shock - s * t_end) < 2 * dx, f"N={n}: shock at {x_shock:.4f}, exact {s * t_end:.4f}"
        # No oscillations. Ahead of the shock: never below h_R (a Lax-Wendroff ripple would
        # undershoot), and back to exactly h_R within 15 cells (the smeared foot is ~9 cells
        # wide on every grid). Behind it: only a tiny overshoot on the plateau.
        ahead = x > s * t_end
        assert np.min(h[ahead]) >= hR - 1e-12, "undershoot ahead of the shock"
        assert np.max(h[x > s * t_end + 15 * dx]) <= hR + 1e-12
        plateau = (x > (um - np.sqrt(core.G * hm)) * t_end + 0.3) & (x < s * t_end - 0.3)
        overshoot = np.max(h[plateau]) - hm
    print(f"  Stoker: h_m = {hm:.4f} m, shock speed {s:.4f} m/s")
    rates = np.log2(np.array(L1[:-1]) / np.array(L1[1:]))
    print("  Stoker L1 errors:", np.array2string(np.array(L1), precision=3))
    print("  Stoker L1 rates: ", np.array2string(rates, precision=2))
    print(f"  Stoker plateau overshoot at N = 1600: {overshoot:.1e} m")
    # L1, not L2, is the criterion: a captured shock limits L2 to rate <= 0.5 on average,
    # with noise from where the jump falls inside a cell.
    assert np.all(rates > 0.7), "Stoker L1 error is not converging at about first order"
    assert overshoot < 0.01 * (hm - hR), "oscillation on the plateau behind the shock"


def test_steady_voellmy():
    theta, phi, xi, h0 = np.radians(30.0), np.radians(20.0), 500.0, 2.0
    mu = np.tan(phi)
    A = core.G * (np.sin(theta) - mu * np.cos(theta))   # driving acceleration (m/s^2)
    B = core.G / (xi * h0)                               # drag coefficient (1/m)
    u_t = np.sqrt(xi * h0 * (np.sin(theta) - mu * np.cos(theta)))
    n, dx = 50, 1.0
    bed = core.uniform_bed(n, theta)
    bc = ("periodic", "periodic")

    # (a) the terminal state is an exact fixed point of one full step
    h = np.full(n, h0)
    hu = h * u_t
    dt = core.cfl_dt(h, hu, bed["cos_c"], dx)
    h1, hu1, _, _ = core.step(h, hu, dt, dx, bed, mu, xi, bc)
    assert np.max(np.abs(hu1 / h1 - u_t)) / u_t < 1e-12

    # (b) transient from rest: du/dt = A - B u^2  ->  u = u_t tanh(sqrt(A B) t)
    errs = []
    for dxk in [1.0, 0.5, 0.25]:   # smaller dx -> smaller CFL dt
        h, hu = np.full(n, h0), np.zeros(n)
        h, hu, _, _, _ = core.run_until(h, hu, 0.0, 5.0, dxk, bed, mu, xi, bc)
        errs.append(abs(hu[0] / h[0] - u_t * np.tanh(np.sqrt(A * B) * 5.0)) / u_t)
    print("  Voellmy transient rel. errors at t=5 s:", ", ".join(f"{e:.1e}" for e in errs))
    assert errs[-1] < 1e-3 and errs[0] > errs[1] > errs[2]

    # (c) long time: velocity reaches the terminal value
    h, hu = np.full(n, h0), np.zeros(n)
    h, hu, _, _, _ = core.run_until(h, hu, 0.0, 100.0, dx, bed, mu, xi, bc)
    rel = abs(hu[0] / h[0] - u_t) / u_t
    print(f"  Voellmy terminal velocity: u_t = {u_t:.4f} m/s, rel. error {rel:.1e}")
    assert rel < 1e-10


def test_mass_conservation():
    # closed domain (walls at both ends) on the app's two-reach bed
    n = 600
    dx = runout.DOMAIN_LENGTH / n
    s_faces = np.linspace(0.0, runout.DOMAIN_LENGTH, n + 1)
    s = 0.5 * (s_faces[1:] + s_faces[:-1])
    th1 = np.radians(35.0)
    bed = {"cos_c": np.cos(runout.two_reach_theta(s, th1)),
           "sin_c": np.sin(runout.two_reach_theta(s, th1)),
           "cos_f": np.cos(runout.two_reach_theta(s_faces, th1))}
    h = np.where(s < 150.0, runout.H0, 0.0)
    hu = np.zeros(n)
    m0 = np.sum(h) * dx
    h, hu, _, _, _ = core.run_until(h, hu, 0.0, 60.0, dx, bed, np.tan(np.radians(15.0)), 300.0,
                                    ("wall", "wall"))
    rel = abs(np.sum(h) * dx - m0) / m0
    print(f"  closed-domain mass error after 60 s: {rel:.1e}")
    assert rel < 1e-13

    # the app's own run (open outlet): volume in domain + outflow = initial volume
    res = runout.run_runout(10000.0, 20.0, 30.0, 20.0, 500.0)
    print(f"  app run: runout {res['runout']:.1f} m, t = {res['t_final']:.1f} s, "
          f"stopped = {res['stopped']}, mass error {res['mass_error']:.1e}")
    assert abs(res["mass_error"]) < 1e-12


if __name__ == "__main__":
    for test in (test_ritter, test_stoker, test_steady_voellmy, test_mass_conservation):
        print(test.__name__)
        test()
    print("all verification tests passed")
