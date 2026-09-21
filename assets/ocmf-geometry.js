/* ──────────────────────────────────────────────────────────────────────────────
   CentrivaLab O-CMF cartridge geometry and ΔP–rpm model
   Single source of truth shared by gamma.html, speed.html and simulator.html.
   Exposes a global `OCMF` in the browser (module.exports in node). Plain script, no build step.

   Pressure at the membrane:  ΔP = ½ ρ ω² (r_m² − r_s²)
     r_m = r_max − δ_m                     (δ_m intrinsic to the cartridge)
     r_s = r_max − δ_s(V)                  (free-surface radius; depends on the feed volume V)
   Over a run from V0 to V0 − Vp at constant speed, the mean of ΔP(V0) and ΔP(V0 − Vp) equals
   ½ ρ ω² [r_m² − ½(r_s0² + r_s1²)], so one "effective" term serves both conversions.
   ────────────────────────────────────────────────────────────────────────────── */
(function (root) {
const OCMF = (() => {
  "use strict";

  const CARTRIDGES = {
    "94-25": {
      name: "CentrivaLab 94/25",
      deltaM: 0.0268,                                            // m, r_max − r_m
      // Measured free-surface curve δ_s(V) = A + B·V + C·V² (mm, V in mL), fitted to 9 CAD
      // measurements spanning 7.4–22.8 mL. Quadratic because the equipotential is a cylinder,
      // not a plane. Cartridge-intrinsic: valid for any 25° rotor. Do not extrapolate.
      ds: { A: 18.7197, B: 2.03995, C: -0.022293 },
      vol: { min: 7.4, max: 22.8, nom: 22.5, vpNom: 0.5 },     // mL
      rotorSelect: true,
      rotors: [
        { id: "sigma12159", label: "Sigma 12159 — 6×94 mL, 25°, r_max 98 mm, 15 500 rpm",              rmax: 0.098, maxRpm: 15500 },
        { id: "fiberlite",  label: "Thermo FiberLite F15-6×100y — 6×100 mL, 25°, r_max 98 mm, 15 000 rpm", rmax: 0.098, maxRpm: 15000 },
        { id: "sigma12155", label: "Sigma 12155 — 4×94 mL, 26°, r_max 91 mm, 20 000 rpm",              rmax: 0.091, maxRpm: 20000 },
        { id: "fx6100",     label: "Beckman FX6100/JA-10.100 — 6×100 mL, 25°, r_max 98 mm, 10 200 rpm", rmax: 0.098, maxRpm: 10200 },
        { id: "vf694",      label: "Beckman VF 6.94 — 6×94 mL, 25°, r_max 106 mm, 10 000 rpm",         rmax: 0.106, maxRpm: 10000 }
      ]
    },
    "50-34": {
      name: "CentrivaLab 50/34",
      deltaM: 0.0318,                                            // m  (r_max − r_m = 0.108 − 0.0762)
      deltaS: 0.0367, alpha: 1.60e-3,                            // m, m/mL: δ_s(V) = δ_s* + α·V (validated vs spec table)
      vol: { min: 11.0, max: 13.0, nom: 13.0, vpNom: 0.5 },
      rotorSelect: false,
      rotors: [ { id: "ja20", label: "Beckman JA-20 — r_max 108 mm, 20 000 rpm", rmax: 0.1080, maxRpm: 20000 } ]
    }
  };

  const rotor = (cart, id) => CARTRIDGES[cart].rotors.find(r => r.id === id) || CARTRIDGES[cart].rotors[0];

  // δ_s(V) in metres
  function deltaS(cart, V) {
    const c = CARTRIDGES[cart];
    return c.ds ? (c.ds.A + c.ds.B * V + c.ds.C * V * V) / 1000 : c.deltaS + c.alpha * V;
  }
  // Inverse of δ_s for the 94/25 (mm in). The physical branch is the root inside the tube.
  function v94FromDs(ds_mm) {
    const { A, B, C } = CARTRIDGES["94-25"].ds;
    const disc = B * B - 4 * C * (A - ds_mm);
    return disc < 0 ? null : (-B + Math.sqrt(disc)) / (2 * C);
  }
  const rm = (cart, rmax) => rmax - CARTRIDGES[cart].deltaM;
  const rs = (cart, rmax, V) => rmax - deltaS(cart, V);

  // Geometry of a run from V0 down to V0 − Vp (Vp = 0 → single operating point)
  function geomRun(cart, rmax, V0, Vp = 0) {
    const r_m = rm(cart, rmax), r_s0 = rs(cart, rmax, V0), r_s1 = rs(cart, rmax, V0 - Vp);
    return { r_m, r_s0, r_s1,
             sq0: r_m * r_m - r_s0 * r_s0, sq1: r_m * r_m - r_s1 * r_s1,
             sqEff: r_m * r_m - 0.5 * (r_s0 * r_s0 + r_s1 * r_s1) };
  }
  const omega = rpm => rpm * 2 * Math.PI / 60;
  const dpFromRpm = (rpm, rho, sq) => 0.5 * rho * omega(rpm) ** 2 * sq / 1e5;          // bar
  const rpmFromDp = (dP, rho, sq) => sq > 0 ? Math.sqrt(2 * dP * 1e5 / (rho * sq)) * 60 / (2 * Math.PI) : 0;
  const amOf = (rpm, r_m) => omega(rpm) ** 2 * r_m;                                      // m/s² at the membrane

  // Water density, kg/m³, T in °C (0–50 °C, ±0.05 kg/m³)
  const waterRho = t => 999.84 + 0.0680 * t - 0.00909 * t * t + 0.0000976 * t ** 3 - 1.7e-6 * t ** 4;

  // Fill a <select> with the rotors of a cartridge (option value = rotor id)
  function fillRotorSelect(sel, cart, keepId) {
    const prev = keepId || sel.value;
    sel.innerHTML = "";
    for (const r of CARTRIDGES[cart].rotors) sel.add(new Option(r.label, r.id));
    if (CARTRIDGES[cart].rotors.some(r => r.id === prev)) sel.value = prev;
  }

  return { CARTRIDGES, rotor, deltaS, v94FromDs, rm, rs, geomRun, omega, dpFromRpm, rpmFromDp, amOf, waterRho, fillRotorSelect, G: 9.81 };
})();
if (typeof module !== "undefined" && module.exports) module.exports = OCMF; else root.OCMF = OCMF;
})(typeof window !== "undefined" ? window : globalThis);
