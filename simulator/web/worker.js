// Web worker: loads Pyodide and the unchanged engine `cosmicmass`, then runs scenes (task UI1).
importScripts("https://cdn.jsdelivr.net/pyodide/v0.26.4/full/pyodide.js");

let pyodide = null;

async function boot() {
  postMessage({ type: "status", text: "Loading Python…" });
  pyodide = await loadPyodide();
  postMessage({ type: "status", text: "Loading NumPy, mpmath, SymPy…" });
  await pyodide.loadPackage(["numpy", "mpmath", "sympy"]);
  postMessage({ type: "status", text: "Loading the engine…" });
  const manifest = await (await fetch("engine-files.json")).json();
  for (const file of manifest.files) {
    const response = await fetch(manifest.root + file);
    if (!response.ok) throw new Error(`cannot load ${file} (${response.status})`);
    const path = "/home/pyodide/engine/" + file;
    pyodide.FS.mkdirTree(path.substring(0, path.lastIndexOf("/")));
    pyodide.FS.writeFile(path, await response.text());
  }
  pyodide.FS.writeFile("/home/pyodide/bridge.py", await (await fetch("bridge.py")).text());
  pyodide.runPython("import sys; sys.path[:0] = ['/home/pyodide/engine', '/home/pyodide']; import bridge");
  const presets = JSON.parse(pyodide.runPython("bridge.presets()"));
  const version = pyodide.runPython("import cosmicmass; cosmicmass.__version__");
  postMessage({ type: "ready", presets, version });
}

const ready = boot().catch((err) => postMessage({ type: "error", text: String(err) }));

onmessage = async (event) => {
  await ready;
  const { id, settings } = event.data;
  try {
    pyodide.globals.set("settings_json", JSON.stringify(settings));
    const out = pyodide.runPython("bridge.run_scene(settings_json)");
    postMessage({ type: "result", id, result: JSON.parse(out) });
  } catch (err) {
    const text = String(err).split("\n").filter(Boolean).slice(-1)[0] || String(err);
    postMessage({ type: "failed", id, text });
  }
};
