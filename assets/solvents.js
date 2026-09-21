/* ──────────────────────────────────────────────────────────────────────────────
   CentrivaLab solvent property table
   Shared by planner.html (and any future OSN-capable tool). Exposes a global `SOLVENTS`
   in the browser and module.exports in node. Plain script, no build step.

   Each solvent: density at 25 °C (kg/m³) with a linear thermal expansion term,
   and dynamic viscosity as an Andrade fit ln μ = A + B/T (μ in mPa·s, T in K)
   through the 20 °C and 40 °C handbook values. Valid roughly 10–50 °C.
   `dense: true` flags solvents denser than most organic solutes (DCM, CHCl₃, DMSO):
   the retained layer is then lighter than the bulk and the O-CMF buoyancy model does not apply.
   Sources: CRC Handbook / DIPPR compilations, rounded to 3 significant figures.
   ────────────────────────────────────────────────────────────────────────────── */
(function (root) {
  "use strict";

  // μ20, μ40 in mPa·s; rho25 in kg/m³; drho in kg/m³/K (≈ −ρ·α_p)
  const RAW = [
    { id: "water",        name: "Water",               M: 18.0,  rho25: 997,  drho: -0.26, mu20: 1.002, mu40: 0.653 },
    { id: "methanol",     name: "Methanol",            M: 32.0,  rho25: 787,  drho: -0.94, mu20: 0.594, mu40: 0.456 },
    { id: "ethanol",      name: "Ethanol",             M: 46.1,  rho25: 785,  drho: -0.85, mu20: 1.194, mu40: 0.834 },
    { id: "ipa",          name: "2-Propanol (IPA)",    M: 60.1,  rho25: 781,  drho: -0.84, mu20: 2.40,  mu40: 1.38  },
    { id: "butanol",      name: "1-Butanol",           M: 74.1,  rho25: 806,  drho: -0.77, mu20: 2.95,  mu40: 1.78  },
    { id: "acetone",      name: "Acetone",             M: 58.1,  rho25: 784,  drho: -1.13, mu20: 0.324, mu40: 0.275 },
    { id: "mek",          name: "Methyl ethyl ketone", M: 72.1,  rho25: 800,  drho: -1.02, mu20: 0.423, mu40: 0.348 },
    { id: "ethylacetate", name: "Ethyl acetate",       M: 88.1,  rho25: 894,  drho: -1.17, mu20: 0.452, mu40: 0.362 },
    { id: "thf",          name: "THF",                 M: 72.1,  rho25: 883,  drho: -1.07, mu20: 0.482, mu40: 0.397 },
    { id: "acetonitrile", name: "Acetonitrile",        M: 41.1,  rho25: 777,  drho: -1.06, mu20: 0.369, mu40: 0.301 },
    { id: "dmf",          name: "DMF",                 M: 73.1,  rho25: 944,  drho: -0.92, mu20: 0.924, mu40: 0.664 },
    { id: "toluene",      name: "Toluene",             M: 92.1,  rho25: 862,  drho: -0.93, mu20: 0.590, mu40: 0.470 },
    { id: "xylene",       name: "Xylene (mixed)",      M: 106.2, rho25: 861,  drho: -0.87, mu20: 0.650, mu40: 0.520 },
    { id: "heptane",      name: "n-Heptane",           M: 100.2, rho25: 680,  drho: -0.85, mu20: 0.418, mu40: 0.341 },
    { id: "hexane",       name: "n-Hexane",            M: 86.2,  rho25: 655,  drho: -0.90, mu20: 0.310, mu40: 0.260 },
    { id: "cyclohexane",  name: "Cyclohexane",         M: 84.2,  rho25: 774,  drho: -0.94, mu20: 0.980, mu40: 0.713 },
    { id: "dcm",          name: "Dichloromethane",     M: 84.9,  rho25: 1317, drho: -1.80, mu20: 0.437, mu40: 0.370, dense: true },
    { id: "chloroform",   name: "Chloroform",          M: 119.4, rho25: 1480, drho: -1.90, mu20: 0.563, mu40: 0.464, dense: true },
    { id: "dmso",         name: "DMSO",                M: 78.1,  rho25: 1096, drho: -0.98, mu20: 2.24,  mu40: 1.48,  dense: true }
  ];

  const T20 = 293.15, T40 = 313.15;
  const SOLVENTS = {};
  for (const s of RAW) {
    const B = Math.log(s.mu20 / s.mu40) / (1 / T20 - 1 / T40);
    const A = Math.log(s.mu20) - B / T20;
    SOLVENTS[s.id] = Object.assign({}, s, { andrade: { A, B } });
  }

  // Properties at temperature t (°C): rho [kg/m³], mu [Pa·s], nu [m²/s]
  function props(id, t) {
    const s = SOLVENTS[id];
    const T = t + 273.15;
    const rho = s.rho25 + s.drho * (t - 25);
    const mu = Math.exp(s.andrade.A + s.andrade.B / T) * 1e-3;
    return { rho, mu, nu: mu / rho, M: s.M, name: s.name, dense: !!s.dense };
  }

  const api = { SOLVENTS, props, list: RAW.map(s => ({ id: s.id, name: s.name })) };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.SOLV = api;
})(typeof window !== "undefined" ? window : globalThis);
