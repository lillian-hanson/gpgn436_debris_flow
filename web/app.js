// Main thread: UI, plotting and playback only. All physics runs in worker.js (Pyodide);
// every number shown here comes from the solver's output.

const $ = (id) => document.getElementById(id);
const worker = new Worker("worker.js");

const INPUTS = {
  volume: { fmt: (v) => `${(+v).toLocaleString()} m³` },
  width: { fmt: (v) => `${v} m` },
  theta1: { fmt: (v) => `${v}°` },
  phi: { fmt: (v) => `${v}°` },
  xi: { fmt: (v) => `${v} m/s²` },
};
// Vertical exaggeration of flow depth in the side view (display only), switchable 1x / 5x.
// At 1x a 3 m flow on ~200 m of relief is only a few pixels thick; at 5x it is visible
// but the deposit looks five times taller relative to the slope than it is.
let exaggeration = 5;
// Depth below which a cell is drawn as dry in the overhead view (m); matches the
// solver's front threshold H_FRONT.
const DRY_DEPTH = 0.01;

let geometry = null;   // {s_faces, bed_x, bed_z, s_break, domain, h0}
let snaps = [];        // [{t, h, u, front}]
let result = null;     // metrics from the solver
let frame = 0;
let hMaxSeen = 1;
let playTimer = null;
let drawQueued = false;

// ---------- inputs ----------
function params() {
  const p = {};
  for (const k of Object.keys(INPUTS)) p[k] = parseFloat($(k).value);
  return p;
}
function updateInputs() {
  for (const [k, spec] of Object.entries(INPUTS)) $(k + "-out").textContent = spec.fmt($(k).value);
  const p = params();
  // Release length L0 = V / (b h0), h0 = 3 m (the solver checks the limit and reports errors)
  $("release-info").textContent = `Release block: ${(p.volume / (p.width * 3)).toFixed(0)} m long × 3 m deep`;
}
for (const k of Object.keys(INPUTS)) $(k).addEventListener("input", updateInputs);
updateInputs();
for (const radio of document.querySelectorAll('input[name="exag"]')) {
  radio.addEventListener("change", (e) => { exaggeration = +e.target.value; queueDraw(); });
}

$("run").addEventListener("click", () => {
  stopPlay();
  snaps = []; result = null; geometry = null; frame = 0; hMaxSeen = 1;
  clearMetrics();
  $("run").disabled = true;
  setStatus("Running…");
  worker.postMessage(params());
});

// ---------- worker messages ----------
worker.onmessage = ({ data: msg }) => {
  if (msg.type === "status") setStatus(msg.text);
  else if (msg.type === "ready") {
    $("run").disabled = false; $("run").textContent = "Run simulation";
    setStatus("Ready. Set inputs and run.");
  } else if (msg.type === "geometry") geometry = msg.data;
  else if (msg.type === "snapshot") {
    snaps.push(msg.data);
    hMaxSeen = Math.max(hMaxSeen, ...msg.data.h);
    frame = snaps.length - 1;
    setStatus(`Running… t = ${msg.data.t.toFixed(1)} s, front at ${msg.data.front.toFixed(0)} m`);
    queueDraw();
  } else if (msg.type === "done") {
    result = msg.data;
    showMetrics(result);
    const slider = $("time");
    slider.max = snaps.length - 1; slider.value = frame; slider.disabled = false;
    $("play").disabled = false;
    $("run").disabled = false;
    setStatus(result.stopped ? `Done: front stopped at t = ${result.t_final.toFixed(1)} s.`
                             : `Stopped at the time limit (${result.t_final.toFixed(0)} s); front still moving.`);
    queueDraw();
  } else if (msg.type === "error") {
    setStatus(msg.text, true);
    $("run").disabled = false;
    if ($("run").textContent === "Loading…") $("run").textContent = "Run simulation";
  }
};

function setStatus(text, isError = false) {
  $("status").textContent = text;
  $("status").classList.toggle("error", isError);
}

// ---------- metrics ----------
const METRIC_IDS = ["m-runout", "m-speed", "m-froude", "m-thick", "m-mass", "m-time"];
function clearMetrics() { for (const id of METRIC_IDS) $(id).textContent = "–"; }
function showMetrics(r) {
  $("m-runout").textContent = r.runout >= geometry.domain
    ? `≥ ${geometry.domain} m (left domain)` : `${r.runout.toFixed(0)} m`;
  $("m-speed").textContent = `${r.peak_front_speed.toFixed(1)} m/s`;
  $("m-froude").textContent = r.peak_froude.toFixed(2);
  $("m-thick").textContent = `${r.max_deposit_thickness.toFixed(2)} m`;
  $("m-mass").textContent = Math.abs(r.mass_error).toExponential(1);
  $("m-time").textContent = `${r.t_final.toFixed(0)} s · ${(100 * r.outflow_fraction).toFixed(0)}% out`;
}

// ---------- playback ----------
$("time").addEventListener("input", (e) => { stopPlay(); frame = +e.target.value; queueDraw(); });
$("play").addEventListener("click", () => {
  if (playTimer) return stopPlay();
  if (frame >= snaps.length - 1) frame = 0;
  $("play").textContent = "Pause";
  playTimer = setInterval(() => {
    frame = Math.min(frame + 1, snaps.length - 1);
    $("time").value = frame;
    queueDraw();
    if (frame === snaps.length - 1) stopPlay();
  }, 50);
});
function stopPlay() {
  clearInterval(playTimer); playTimer = null;
  $("play").textContent = "Play";
}

// ---------- drawing helpers ----------
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

function setupCanvas(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth, h = canvas.clientHeight;
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
    canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
  }
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  ctx.font = "11px system-ui, sans-serif";
  return { ctx, w, h };
}

function niceTicks(min, max, count) {
  const span = max - min || 1;
  const step0 = span / count, mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= step0);
  const ticks = [];
  for (let v = Math.ceil(min / step) * step; v <= max + 1e-9; v += step) ticks.push(+v.toFixed(10));
  return ticks;
}

// Axes with recessive grid; returns scale functions for the plot box.
function axes(ctx, w, h, xr, yr, xlabel, ylabel) {
  const m = { l: 48, r: 12, t: 8, b: 32 };
  const X = (v) => m.l + (v - xr[0]) / (xr[1] - xr[0]) * (w - m.l - m.r);
  const Y = (v) => h - m.b - (v - yr[0]) / (yr[1] - yr[0]) * (h - m.t - m.b);
  ctx.strokeStyle = css("--grid"); ctx.lineWidth = 1; ctx.fillStyle = css("--text-muted");
  ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const v of niceTicks(xr[0], xr[1], Math.max(3, Math.floor(w / 90)))) {
    ctx.beginPath(); ctx.moveTo(X(v), m.t); ctx.lineTo(X(v), h - m.b); ctx.stroke();
    ctx.fillText(v, X(v), h - m.b + 4);
  }
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (const v of niceTicks(yr[0], yr[1], 4)) {
    ctx.beginPath(); ctx.moveTo(m.l, Y(v)); ctx.lineTo(w - m.r, Y(v)); ctx.stroke();
    ctx.fillText(v, m.l - 6, Y(v));
  }
  ctx.fillStyle = css("--text-secondary"); ctx.textAlign = "center"; ctx.textBaseline = "bottom";
  ctx.fillText(xlabel, (m.l + w - m.r) / 2, h);
  ctx.save(); ctx.translate(11, (m.t + h - m.b) / 2); ctx.rotate(-Math.PI / 2);
  ctx.textBaseline = "middle"; ctx.fillText(ylabel, 0, 0); ctx.restore();
  return { X, Y, m };
}

const centres = (faces) => faces.slice(0, -1).map((v, i) => 0.5 * (v + faces[i + 1]));

// ---------- plots ----------
function queueDraw() {
  if (drawQueued) return;
  drawQueued = true;
  requestAnimationFrame(() => { drawQueued = false; drawAll(); });
}

function drawAll() {
  $("time-out").textContent = snaps.length ? `${snaps[frame].t.toFixed(1)} s` : "0.0 s";
  for (const id of ["side", "overhead", "deposit", "front"]) setupCanvas($(id));
  if (!geometry || !snaps.length) return;
  drawSide(); drawOverhead(); drawDeposit(); drawFront();
}

function drawSide() {
  const { ctx, w, h } = setupCanvas($("side"));
  const snap = snaps[frame];
  const xc = centres(geometry.bed_x), zc = centres(geometry.bed_z);
  const zmin = Math.min(...geometry.bed_z) - 5;
  const zmax = geometry.bed_z[0] + exaggeration * hMaxSeen + 5;
  const { X, Y } = axes(ctx, w, h, [0, geometry.bed_x.at(-1)], [zmin, zmax],
                        "horizontal distance (m)", "elevation (m)");
  // ground
  ctx.fillStyle = css("--bed-fill");
  ctx.beginPath(); ctx.moveTo(X(0), Y(zmin));
  geometry.bed_x.forEach((x, i) => ctx.lineTo(X(x), Y(geometry.bed_z[i])));
  ctx.lineTo(X(geometry.bed_x.at(-1)), Y(zmin)); ctx.closePath(); ctx.fill();
  // flow: bed + exaggerated depth
  ctx.fillStyle = css("--flow-fill"); ctx.strokeStyle = css("--flow"); ctx.lineWidth = 2;
  ctx.beginPath();
  xc.forEach((x, i) => ctx.lineTo(X(x), Y(zc[i] + exaggeration * snap.h[i])));
  for (let i = xc.length - 1; i >= 0; i--) ctx.lineTo(X(xc[i]), Y(zc[i]));
  ctx.closePath(); ctx.fill();
  ctx.beginPath();
  xc.forEach((x, i) => snap.h[i] >= DRY_DEPTH ? ctx.lineTo(X(x), Y(zc[i] + exaggeration * snap.h[i])) : ctx.moveTo(X(x), Y(zc[i])));
  ctx.stroke();
  // bed line
  ctx.strokeStyle = css("--bed"); ctx.lineWidth = 1.5; ctx.beginPath();
  geometry.bed_x.forEach((x, i) => ctx.lineTo(X(x), Y(geometry.bed_z[i]))); ctx.stroke();
  // slope break label
  const ib = geometry.s_faces.findIndex((s) => s >= geometry.s_break);
  ctx.fillStyle = css("--text-muted"); ctx.textAlign = "left"; ctx.textBaseline = "bottom";
  ctx.fillText("slope break (s = 350 m)", X(geometry.bed_x[ib]) + 4, Y(geometry.bed_z[ib]) - 4);
}

function drawOverhead() {
  const { ctx, w, h } = setupCanvas($("overhead"));
  const snap = snaps[frame];
  const L = geometry.domain, ml = 48, mr = 12;
  const X = (s) => ml + s / L * (w - ml - mr);
  // strip height proportional to b (50 m slider maximum = 56 px)
  const bh = Math.max(4, params().width / 50 * 56), y0 = (h - 22 - bh) / 2;
  const flow = css("--flow");
  ctx.fillStyle = css("--dry"); ctx.fillRect(X(0), y0, X(L) - X(0), bh);
  geometry.s_faces.slice(0, -1).forEach((s, i) => {
    if (snap.h[i] < DRY_DEPTH) return;
    // single-hue sequential encoding: opacity of the flow colour ∝ depth
    ctx.globalAlpha = 0.15 + 0.85 * Math.min(1, snap.h[i] / hMaxSeen);
    ctx.fillStyle = flow;
    ctx.fillRect(X(s), y0, X(geometry.s_faces[i + 1]) - X(s) + 0.5, bh);
  });
  ctx.globalAlpha = 1;
  ctx.strokeStyle = css("--text-muted"); ctx.setLineDash([3, 3]);
  ctx.beginPath(); ctx.moveTo(X(geometry.s_break), y0 - 4); ctx.lineTo(X(geometry.s_break), y0 + bh + 4); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = css("--text-muted"); ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const s of niceTicks(0, L, 6)) ctx.fillText(s, X(s), h - 18);
  ctx.textAlign = "left"; ctx.textBaseline = "middle";
  ctx.fillText(`b = ${params().width} m`, 2, y0 + bh / 2);
  ctx.textAlign = "right"; ctx.textBaseline = "top";
  ctx.fillText(`distance along bed s (m) · darkest = ${hMaxSeen.toFixed(1)} m depth`, w - mr, 0);
}

let depositScale = null, frontScale = null;

function drawDeposit() {
  const { ctx, w, h } = setupCanvas($("deposit"));
  const snap = snaps[frame], s = centres(geometry.s_faces);
  const sc = axes(ctx, w, h, [0, geometry.domain], [0, Math.ceil(hMaxSeen)], "distance along bed s (m)", "depth h (m)");
  depositScale = sc;
  // vertical stretch of this chart vs. true shape: pixels per metre of depth / pixels per metre along s
  $("stretch").textContent = `×${Math.round((sc.Y(0) - sc.Y(1)) / (sc.X(1) - sc.X(0)))}`;
  ctx.fillStyle = css("--flow-fill"); ctx.strokeStyle = css("--flow"); ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(sc.X(0), sc.Y(0));
  s.forEach((x, i) => ctx.lineTo(sc.X(x), sc.Y(snap.h[i])));
  ctx.lineTo(sc.X(geometry.domain), sc.Y(0)); ctx.closePath(); ctx.fill();
  ctx.beginPath(); s.forEach((x, i) => ctx.lineTo(sc.X(x), sc.Y(snap.h[i]))); ctx.stroke();
}

function drawFront() {
  const { ctx, w, h } = setupCanvas($("front"));
  const tmax = Math.max(10, snaps.at(-1).t);
  const sc = axes(ctx, w, h, [0, tmax], [0, geometry.domain], "time (s)", "front position s (m)");
  frontScale = sc;
  ctx.strokeStyle = css("--flow"); ctx.lineWidth = 2; ctx.beginPath();
  snaps.forEach((p) => ctx.lineTo(sc.X(p.t), sc.Y(p.front))); ctx.stroke();
  // current-time marker, with a surface ring so it stays visible on the line
  const p = snaps[frame];
  ctx.fillStyle = css("--accent"); ctx.strokeStyle = css("--surface-1"); ctx.lineWidth = 2;
  ctx.beginPath(); ctx.arc(sc.X(p.t), sc.Y(p.front), 5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
}

// ---------- hover tooltips ----------
function hover(canvasId, read) {
  const c = $(canvasId), tip = $("tooltip");
  c.addEventListener("mousemove", (e) => {
    const r = c.getBoundingClientRect();
    const text = snaps.length && geometry ? read(e.clientX - r.left) : null;
    if (!text) { tip.hidden = true; return; }
    tip.textContent = text; tip.hidden = false;
    tip.style.left = `${e.clientX + 12}px`; tip.style.top = `${e.clientY + 12}px`;
  });
  c.addEventListener("mouseleave", () => { tip.hidden = true; });
}
const invX = (sc, px, x0, x1) => x0 + (px - sc.X(x0)) / (sc.X(x1) - sc.X(x0)) * (x1 - x0);

hover("deposit", (px) => {
  if (!depositScale) return null;
  const s = invX(depositScale, px, 0, geometry.domain);
  const i = Math.floor(s / geometry.domain * snaps[frame].h.length);
  if (i < 0 || i >= snaps[frame].h.length) return null;
  return `s = ${s.toFixed(0)} m · h = ${snaps[frame].h[i].toFixed(2)} m · u = ${snaps[frame].u[i].toFixed(2)} m/s`;
});
hover("front", (px) => {
  if (!frontScale) return null;
  const t = invX(frontScale, px, 0, Math.max(10, snaps.at(-1).t));
  let best = 0;
  snaps.forEach((p, i) => { if (Math.abs(p.t - t) < Math.abs(snaps[best].t - t)) best = i; });
  return `t = ${snaps[best].t.toFixed(1)} s · front at ${snaps[best].front.toFixed(1)} m`;
});

window.addEventListener("resize", queueDraw);
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", queueDraw);
