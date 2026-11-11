// Web universe (task UI1): builds settings for the engine, runs them in a worker, animates the result.
"use strict";

const G_N = 6.67430e-11;
const K = 1 / (4 * Math.PI * 8.8541878188e-12);
const $ = (id) => document.getElementById(id);

let PRESETS = {};
let lastSettings = null, lastResult = null, compareResult = null;
let frame = 0, playing = false, nextId = 1;
const pending = new Map();

// ---------- scenes: pairs in their centre-of-mass frame ----------
function pair(a, b, r, vrel) {
  const ma = PRESETS[a].m, mb = PRESETS[b].m, M = ma + mb;
  return [
    { preset: a, x: [-mb / M * r, 0], v: [0, -mb / M * vrel], fixed: false },
    { preset: b, x: [ma / M * r, 0], v: [0, ma / M * vrel], fixed: false },
  ];
}
const strength = (a, b) => G_N * PRESETS[a].m * PRESETS[b].m - K * PRESETS[a].q * PRESETS[b].q;
const reduced = (a, b) => PRESETS[a].m * PRESETS[b].m / (PRESETS[a].m + PRESETS[b].m);
function circular(a, b, r) {
  const k = strength(a, b), mu = reduced(a, b);
  return { v: Math.sqrt(k / (mu * r)), period: 2 * Math.PI * Math.sqrt(mu * r ** 3 / k) };
}

const SCENES = {
  "hydrogen": {
    title: "Hydrogen-like pair (demo T1.4)",
    note: "A proton and an electron on the circle of L0. Switch to L1a to see the charge share at work.",
    build() { const c = circular("proton", "electron", 5.29177e-11);
      return { bodies: pair("proton", "electron", 5.29177e-11, c.v), t_end: 4.5 * c.period, eta: 1e-3, law: "L1a" }; },
  },
  "earth-moon": {
    title: "Earth and Moon, one month",
    note: "Two neutral bodies: L0 gives Newton's orbit, and L1a adds nothing (no charge, no Composer force).",
    build() { const c = circular("Earth", "moon", 3.844e8);
      return { bodies: pair("Earth", "moon", 3.844e8, c.v), t_end: c.period, eta: 1e-3, law: "L0" }; },
  },
  "protons": {
    title: "Two protons released at rest",
    note: "Coulomb repulsion minus a tiny gravity. Test T0.1c compares this with a closed form.",
    build() { const k = -strength("proton", "proton"), mu = reduced("proton", "proton");
      return { bodies: pair("proton", "proton", 1e-10, 0), t_end: 3 * Math.sqrt(mu * 1e-30 / k), eta: 1e-3, law: "L0" }; },
  },
  "fall": {
    title: "Earth and Moon falling from rest",
    note: "Released at rest, they fall towards each other. Stops before contact.",
    build() { const gm = strength("Earth", "moon") / reduced("Earth", "moon");
      const t = 0.9 * (Math.PI / 2) * Math.sqrt((3.844e8) ** 3 / (2 * gm));
      return { bodies: pair("Earth", "moon", 3.844e8, 0), t_end: t, eta: 1e-3, law: "L0" }; },
  },
};

// ---------- worker ----------
const worker = new Worker("worker.js");
worker.onmessage = ({ data }) => {
  if (data.type === "status") setStatus(data.text);
  else if (data.type === "ready") {
    PRESETS = data.presets; $("version").textContent = data.version;
    setStatus("Engine ready", "ok"); $("run").disabled = false; initScenes();
  } else if (data.type === "error") setStatus("Could not start the engine: " + data.text, "err");
  else if (data.type === "result" || data.type === "failed") {
    const p = pending.get(data.id); pending.delete(data.id);
    data.type === "result" ? p.resolve(data.result) : p.reject(new Error(data.text));
  }
};
function runInWorker(settings) {
  const id = nextId++;
  return new Promise((resolve, reject) => { pending.set(id, { resolve, reject }); worker.postMessage({ id, settings }); });
}
function setStatus(text, kind = "") { const s = $("status"); s.textContent = text; s.className = "status " + kind; }

// ---------- editors ----------
function initScenes() {
  const sel = $("scene");
  for (const [key, s] of Object.entries(SCENES)) sel.add(new Option(s.title, key));
  sel.add(new Option("Custom (edit below)", "custom"));
  sel.onchange = () => loadScene(sel.value);
  loadScene("hydrogen");
}
function loadScene(key) {
  if (key === "custom") { $("scene-note").textContent = "Your own universe: up to 16 bodies, any law, any angle."; return; }
  const s = SCENES[key], built = s.build();
  $("scene-note").textContent = s.note;
  $("t-end").value = built.t_end.toPrecision(6);
  $("eta").value = built.eta;
  document.querySelector(`input[name=law][value=${built.law}]`).checked = true;
  lawChanged();
  renderBodies(built.bodies);
}
function bodyColor(q) { return q > 0 ? "var(--pos)" : q < 0 ? "var(--neg)" : "var(--neutral)"; }
function renderBodies(bodies) {
  const box = $("bodies"); box.innerHTML = "";
  bodies.forEach((b) => box.appendChild(bodyRow(b)));
}
const short = (v) => Number(Number(v).toPrecision(7));  // display precision of positions and speeds
function bodyRow(b) {
  b = { ...b, x: b.x.map(short), v: b.v.map(short) };
  const div = document.createElement("div"); div.className = "body";
  const names = Object.keys(PRESETS);
  const m = b.preset ? PRESETS[b.preset].m : b.m, q = b.preset ? PRESETS[b.preset].q : b.q;
  div.innerHTML = `
    <div class="body-head"><span class="dot"></span>
      <select class="kind">${names.map((n) => `<option ${n === b.preset ? "selected" : ""}>${n}</option>`).join("")}
        <option value="" ${b.preset ? "" : "selected"}>custom m, q</option></select>
      <button class="ghost icon del" title="remove">✕</button></div>
    <div class="grid4">
      <label>m (kg)<input class="m" type="number" step="any" value="${m}"></label>
      <label>q (C)<input class="q" type="number" step="any" value="${q}"></label>
      <label>x (m)<input class="x" type="number" step="any" value="${b.x[0]}"></label>
      <label>y (m)<input class="y" type="number" step="any" value="${b.x[1]}"></label>
      <label>vx (m/s)<input class="vx" type="number" step="any" value="${b.v[0]}"></label>
      <label>vy (m/s)<input class="vy" type="number" step="any" value="${b.v[1]}"></label>
      <label class="check wide"><input class="fixed" type="checkbox" ${b.fixed ? "checked" : ""}> fixed in place</label>
    </div>`;
  const kind = div.querySelector(".kind"), mi = div.querySelector(".m"), qi = div.querySelector(".q");
  const sync = () => { const custom = !kind.value; mi.disabled = qi.disabled = !custom;
    if (!custom) { mi.value = PRESETS[kind.value].m; qi.value = PRESETS[kind.value].q; }
    div.querySelector(".dot").style.background = bodyColor(Number(qi.value)); $("scene").value = $("scene").value; };
  kind.onchange = () => { sync(); $("scene").value = "custom"; };
  [mi, qi].forEach((i) => (i.oninput = () => (div.querySelector(".dot").style.background = bodyColor(Number(qi.value)))));
  div.querySelectorAll("input").forEach((i) => i.addEventListener("change", () => ($("scene").value = "custom")));
  div.querySelector(".del").onclick = () => { div.remove(); $("scene").value = "custom"; };
  sync();
  return div;
}
$("add-body").onclick = () => {
  if ($("bodies").children.length >= 16) return;
  $("bodies").appendChild(bodyRow({ preset: "proton", x: [0, 0], v: [0, 0], fixed: false }));
  $("scene").value = "custom";
};
function lawChanged() {
  const l1a = document.querySelector("input[name=law]:checked").value === "L1a";
  $("l1a-options").hidden = !l1a;
}
document.querySelectorAll("input[name=law]").forEach((r) => (r.onchange = lawChanged));

// ---------- settings ----------
function readSettings(lawOverride) {
  const law = lawOverride || document.querySelector("input[name=law]:checked").value;
  const bodies = [...$("bodies").children].map((div, i) => {
    const kind = div.querySelector(".kind").value, num = (c) => Number(div.querySelector(c).value);
    const b = { x: [num(".x"), num(".y")], v: [num(".vx"), num(".vy")] };
    if (kind) b.preset = kind; else { b.name = `body ${i + 1}`; b.m = num(".m"); b.q = num(".q"); }
    if (div.querySelector(".fixed").checked) b.fixed = true;
    return b;
  });
  if (bodies.length < 1) throw new Error("add at least one body");
  const lawSpec = law === "L0" ? { id: "L0" }
    : { id: "L1a", lambda: Number($("lambda").value), sign: $("sign").value, terms: $("terms").value };
  const numerics = { t_end: Number($("t-end").value), eta: Number($("eta").value), record_every: 1 };
  if (!$("soft").checked) numerics.softening = 0.0;
  return {
    scenario: "web-" + $("scene").value, universe: { phi_bar_deg: Number($("phi").value), dimension: 2, units: "SI" },
    law: lawSpec, bodies, numerics, kind: law === "L0" ? "mathematics" : "finding", preregistered: false,
  };
}

// ---------- run ----------
$("run").onclick = async () => {
  let settings;
  try { settings = readSettings(); } catch (e) { setStatus(e.message, "err"); return; }
  const compare = settings.law.id === "L1a" && $("compare").checked;
  $("run").disabled = true; stop();
  const started = performance.now();
  const tick = setInterval(() => setStatus(`Running the engine… ${((performance.now() - started) / 1000).toFixed(0)} s`), 500);
  try {
    lastResult = await runInWorker(settings);
    compareResult = compare ? await runInWorker(readSettings("L0")) : null;
    lastSettings = settings;
    setStatus(`Done in ${((performance.now() - started) / 1000).toFixed(1)} s: ${lastResult.steps.toLocaleString()} steps`, "ok");
    showResult();
  } catch (e) { setStatus("The engine refused: " + e.message, "err"); }
  finally { clearInterval(tick); $("run").disabled = false; }
};

function showResult() {
  $("badge").hidden = lastResult.law !== "L1a";
  $("dl-result").disabled = false;
  const n = lastResult.t.length - 1;
  Object.assign($("scrub"), { max: n, value: 0, disabled: false }); $("play").disabled = false;
  frame = 0; computeView(); renderReadout(); renderLegend(); play();
}

// ---------- drawing ----------
let view = null;
function computeView() {
  const pts = [lastResult, compareResult].filter(Boolean).flatMap((r) => r.paths.flat());
  let [x0, x1, y0, y1] = [Infinity, -Infinity, Infinity, -Infinity];
  for (const [x, y] of pts) { x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y); }
  const span = Math.max(x1 - x0, y1 - y0) || 1;
  view = { cx: (x0 + x1) / 2, cy: (y0 + y1) / 2, span: span * 1.15 };
}
function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
function draw() {
  const c = $("view"), ctx = c.getContext("2d"), W = c.width, H = c.height;
  ctx.clearRect(0, 0, W, H);
  if (!lastResult || !view) return;
  const s = Math.min(W, H) / view.span;
  const P = ([x, y]) => [W / 2 + (x - view.cx) * s, H / 2 - (y - view.cy) * s];
  const masses = lastSettings.bodies.map((b) => (b.preset ? PRESETS[b.preset].m : b.m));
  const charges = lastSettings.bodies.map((b) => (b.preset ? PRESETS[b.preset].q : b.q));
  const radius = (m) => { const lo = Math.log10(Math.min(...masses)), hi = Math.log10(Math.max(...masses));
    return hi === lo ? 7 : 4 + 8 * (Math.log10(m) - lo) / (hi - lo); };
  const layers = compareResult ? [[compareResult, true], [lastResult, false]] : [[lastResult, false]];
  for (const [res, ghost] of layers) {
    const f = Math.min(frame, res.t.length - 1);
    res.paths.forEach((path, i) => {
      const col = ghost ? css("--l0") : charges[i] > 0 ? css("--pos") : charges[i] < 0 ? css("--neg") : css("--neutral");
      ctx.strokeStyle = col; ctx.lineWidth = ghost ? 1.2 : 1.6; ctx.setLineDash(ghost ? [5, 4] : []);
      ctx.beginPath(); path.slice(0, f + 1).forEach((p, k) => (k ? ctx.lineTo(...P(p)) : ctx.moveTo(...P(p)))); ctx.stroke();
      ctx.setLineDash([]); ctx.fillStyle = col; ctx.beginPath(); ctx.arc(...P(path[f]), radius(masses[i]) * (ghost ? 0.7 : 1), 0, 2 * Math.PI); ctx.fill();
    });
  }
  // scale bar
  const bar = niceLength(view.span / 5), px = bar * s;
  ctx.strokeStyle = css("--muted"); ctx.fillStyle = css("--muted"); ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(20, H - 22); ctx.lineTo(20 + px, H - 22); ctx.stroke();
  ctx.font = "12px system-ui"; ctx.fillText(fmt(bar, "m"), 20, H - 30);
  $("clock").textContent = "t = " + fmt(lastResult.t[Math.min(frame, lastResult.t.length - 1)], "s");
  $("scrub").value = frame;
}
function niceLength(x) { const p = 10 ** Math.floor(Math.log10(x)); return [1, 2, 5, 10].map((k) => k * p).find((v) => v >= x) || x; }
function fmt(v, unit) { return (Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-2 && v !== 0) ? v.toExponential(3) : v.toPrecision(4)) + " " + unit; }

function play() { playing = true; $("play").textContent = "❚❚"; requestAnimationFrame(loop); }
function stop() { playing = false; $("play").textContent = "▶"; }
function loop() {
  if (!playing || !lastResult) return;
  frame += Math.max(1, Math.round(Number($("speed").value) / 5));
  if (frame >= lastResult.t.length - 1) { frame = lastResult.t.length - 1; draw(); stop(); return; }
  draw(); requestAnimationFrame(loop);
}
$("play").onclick = () => { if (playing) stop(); else { if (frame >= lastResult.t.length - 1) frame = 0; play(); } };
$("scrub").oninput = () => { stop(); frame = Number($("scrub").value); draw(); };

// ---------- readout ----------
function renderLegend() {
  const items = lastSettings.bodies.map((b, i) => {
    const q = b.preset ? PRESETS[b.preset].q : b.q;
    return `<span><i class="dot" style="background:${bodyColor(q)}"></i>${b.preset || b.name}</span>`;
  });
  if (compareResult) items.push(`<span><i class="dot" style="background:var(--l0)"></i>dashed: the same start under L0</span>`);
  $("legend").innerHTML = items.join("");
}
function renderReadout() {
  const r = lastResult, d = r.drift, e = (x) => (x == null ? "n/a" : Number(x).toExponential(2));
  const cells = [
    ["law", r.law + (r.law === "L1a" ? " (hypothesis)" : " (reference)")],
    ["steps", r.steps.toLocaleString()],
    ["energy (Re) change at end", e(d.energy_re)],
    ["largest energy excursion", e(r.energy_excursion)],
    ["momentum drift", e(d.momentum)],
    ["mass / charge drift", `${e(d.mass)} / ${e(d.charge)}`],
  ];
  let html = cells.map(([k, v]) => `<div><small>${k}</small><b>${v}</b></div>`).join("");
  if (r.law === "L1a") html += `<div class="wide"><small>Note</small> Under L1a the L0 energy is not a conserved quantity, so its change is not an error. The Composer energy needs an approved derivation (not in 0.1.0).</div>`;
  html += `<div><small>engine</small><b>cosmicmass ${r.version}</b></div>`;
  $("readout").innerHTML = html;
}

// ---------- downloads ----------
function download(name, obj) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([JSON.stringify(obj, null, 2)], { type: "application/json" }));
  a.download = name; a.click(); URL.revokeObjectURL(a.href);
}
$("dl-settings").onclick = () => { try { download("settings.json", readSettings()); } catch (e) { setStatus(e.message, "err"); } };
$("dl-result").onclick = () => lastResult && download("result.json", { settings: lastSettings, stamp: lastResult.stamp,
  steps: lastResult.steps, drift: lastResult.drift, energy_excursion: lastResult.energy_excursion });

window.addEventListener("resize", draw);
