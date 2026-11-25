// Web Worker: loads Pyodide + NumPy, installs the solver/ package, runs simulations.
// The main thread never runs the solver loop; it only receives messages from here.
//
// Messages to the main thread:
//   {type: "status", text}                   loading progress
//   {type: "ready"}                          Pyodide and the solver are loaded
//   {type: "geometry", data}                 grid and bed profile for the current run
//   {type: "snapshot", data}                 {t, h, u, front} every SNAP_DT of simulated time
//   {type: "done", data}                     metrics and front history
//   {type: "error", text}

// Pinned Pyodide version, loaded from the official CDN (no backend needed).
const PYODIDE_VERSION = "0.28.3";
const PYODIDE_URL = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;
// Solver files, fetched relative to this worker so the same paths work on localhost and
// under GitHub Pages' /<repo-name>/ prefix.
const SOLVER_FILES = ["__init__.py", "core.py", "runout.py"];

importScripts(PYODIDE_URL + "pyodide.js");

// Thin glue: calls the real solver and serializes results as JSON (no physics here).
const GLUE = `
import json
import numpy as np
from solver import runout

def run(params_json, post):
    p = json.loads(params_json)
    th1 = np.radians(p["theta1"])
    s_faces = np.linspace(0.0, runout.DOMAIN_LENGTH, runout.N_CELLS + 1)
    X, Z = runout.bed_profile(s_faces, runout.two_reach_theta(s_faces, th1))
    post("geometry", json.dumps({"s_faces": s_faces.tolist(), "bed_x": X.tolist(),
                                 "bed_z": Z.tolist(), "s_break": runout.S_BREAK,
                                 "domain": runout.DOMAIN_LENGTH, "h0": runout.H0}))

    def on_snapshot(t, h, u, front):
        post("snapshot", json.dumps({"t": t, "h": np.round(h, 5).tolist(),
                                     "u": np.round(u, 4).tolist(), "front": front}))

    res = runout.run_runout(p["volume"], p["width"], p["theta1"], p["phi"], p["xi"],
                            on_snapshot=on_snapshot)
    keys = ["times", "fronts", "release_length", "t_final", "stopped", "runout",
            "peak_front_speed", "peak_froude", "max_deposit_thickness",
            "outflow_fraction", "mass_error"]
    out = {k: (res[k].tolist() if hasattr(res[k], "tolist") else res[k]) for k in keys}
    post("done", json.dumps(out))
`;

async function init() {
  postMessage({ type: "status", text: "Loading Python (Pyodide " + PYODIDE_VERSION + ")…" });
  const pyodide = await loadPyodide({ indexURL: PYODIDE_URL });
  postMessage({ type: "status", text: "Loading NumPy…" });
  await pyodide.loadPackage("numpy");

  // Copy solver/*.py into Pyodide's virtual file system as an importable package.
  pyodide.FS.mkdirTree("/home/pyodide/solver");
  for (const f of SOLVER_FILES) {
    const resp = await fetch("solver/" + f);
    if (!resp.ok) throw new Error(`could not load solver/${f} (${resp.status})`);
    pyodide.FS.writeFile("/home/pyodide/solver/" + f, await resp.text());
  }
  pyodide.runPython("import sys; sys.path.insert(0, '/home/pyodide')");
  pyodide.runPython(GLUE);
  return pyodide;
}

const ready = init().then(
  (py) => { postMessage({ type: "ready" }); return py; },
  (err) => { postMessage({ type: "error", text: String(err) }); throw err; }
);

const post = (type, json) => postMessage({ type, data: JSON.parse(json) });

onmessage = async (event) => {
  const pyodide = await ready;
  try {
    pyodide.globals.get("run")(JSON.stringify(event.data), post);
  } catch (err) {
    // Python exceptions (e.g. release block too long) arrive here; keep the last line.
    const lines = String(err.message || err).trim().split("\n");
    postMessage({ type: "error", text: lines[lines.length - 1] });
  }
};
