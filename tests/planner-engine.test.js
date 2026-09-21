// node tests/planner-engine.test.js — sanity checks of the Test Planner engine
const SOLV = require("../assets/solvents.js");
const P = require("../assets/planner-engine.js");
const OCMF = require("../assets/ocmf-geometry.js");
const assert = (c, m) => { if (!c) { console.error("FAIL:", m); process.exitCode = 1; } else console.log("ok  ", m); };
const f = (x, d = 3) => Number(x).toFixed(d);

// ── solvent table ──
const tol = SOLV.props("toluene", 25);
console.log("toluene 25 °C:", f(tol.rho, 0), "kg/m³", f(tol.mu * 1e3, 3), "mPa·s");
assert(Math.abs(tol.mu * 1e3 - 0.56) < 0.02, "toluene viscosity at 25 °C ≈ 0.56 mPa·s");
const w = SOLV.props("water", 25);
assert(Math.abs(w.mu * 1e3 - 0.89) < 0.02, "water viscosity at 25 °C ≈ 0.89 mPa·s");

// ── solute in toluene: 500 Da, ρ_s 1300 ──
const su = P.soluteProps(tol, 500, 1300, 25);
console.log("D =", su.D.toExponential(2), "beta =", su.beta.toExponential(2));
assert(su.D > 5e-10 && su.D < 1e-9, "Stokes–Einstein D of a 500 Da solute in toluene ~ 7e-10");
assert(Math.abs(su.beta - 3.84e-4) < 0.1e-4, "beta ≈ 3.8e-4 (g/L)^-1");

// ── spiral: k at u = 0.15 m/s, 28 mil, D = 8e-10 (reference numbers of the worked example) ──
const t = P.swK(0.15, P.SPACERS["28mil"], 6.81e-7, 8e-10);
console.log("S&M u=0.15: Re", f(t.Re, 0), "k", f(t.k * P.LMH, 0), "LMH");
assert(Math.abs(t.k * P.LMH - 114 * (0.95 / 1.0) ** (0.875 - 1)) < 15, "k ≈ 110–120 LMH for d_h≈1 mm");

// ── module simulation: 4040, 28 mil, Q 2 m³/h, Y 10 %, 20 bar, C0 0.1 g/L, R 1 ──
const inp = { module: P.MODULES["4040"], spacer: P.SPACERS["28mil"], nElements: 1, Q_m3h: 2, dPin_bar: 20, T: 25,
  solvent: tol, solute: { D: 8e-10, beta: 3.84e-4, M: 500, i: 1 }, C0: 0.1, Rint: 1.0, mode: "recovery", Y: 0.10 };
const m = P.simulateModule(inp);
console.log("module: uIn", f(m.uIn, 3), "uOut", f(m.uOut, 3), "Jmean", f(m.Jmean_lmh, 1), "LMH", "Lp", f(m.Lp_lmhbar, 2),
  "cmIn", f(m.cmIn), "cmOut", f(m.cmOut), "cmMeanA", f(m.cmMeanA), "Γin", f(m.GammaIn, 2), "Γout", f(m.GammaOut, 2),
  "Δp_ch", f(m.dpChannel / 1e5, 2), "bar", "Y", f(m.Y, 4), "Robs", f(m.Robs, 3));
assert(Math.abs(m.Y - 0.10) < 1e-3, "recovery solved to 10 %");
assert(m.Jmean_lmh > 20 && m.Jmean_lmh < 30, "J = Y·Q/A = 0.1·2000/7.9 ≈ 25 LMH");
assert(m.uOut < m.uIn && m.cmOut > m.cmIn, "velocity falls and c_m rises along the module");
assert(m.cbOut > 0.1 / 0.9 * 0.99, "retentate ≈ C0/(1−Y) for R=1");
assert(m.flags.length === 0, "no flags for a standard case: " + m.flags.join(" | "));

// permeance mode round trip
const m2 = P.simulateModule(Object.assign({}, inp, { mode: "permeance", Lp_lmhbar: m.Lp_lmhbar }));
assert(Math.abs(m2.Y - m.Y) < 1e-6, "permeance mode reproduces the recovery");

// ── CentrivaLab plan for the module average ──
const rot = OCMF.rotor("94-25", "sigma12159");
const plan = P.planCentrivaLab({ cart: "94-25", rmax: rot.rmax, maxRpm: rot.maxRpm, V0: 22.5, Vp: 0.5, T: 25,
  solvent: tol, solute: { D: 8e-10, beta: 3.84e-4 }, Rint: 1.0, target: { dPh: m.dPmean, J: m.Jmean, cm: m.cmMeanA } });
console.log("plan: rpm", f(plan.rpm, 0), "a/g", f(plan.ag, 0), "C0cl", f(plan.C0, 4), "g/L", "Γcl", f(plan.Gamma, 2), "cm", f(plan.cm),
  "k", f(plan.k * P.LMH, 0), "t_run", f(plan.tRun_min, 1), "min", plan.flags);
assert(Math.abs(plan.cm - m.cmMeanA) / m.cmMeanA < 1e-6, "tube c_m matches the module mean c_m");
assert(plan.rpm > 11000 && plan.rpm < 12500, "≈ 11 700 rpm for 20 bar in the 94/25 + Sigma 12159");
assert(!plan.over, "within rotor limit");

// ── worked example of the note: spiral u=0.15, C0 0.1, J 50 LMH → CentrivaLab C0 ≈ 0.085 g/L at 20 bar ──
const kS = P.swK(0.15, { dh: 1.0e-3, eps: 0.89, h: 0.71e-3 }, 6.81e-7, 8e-10).k;
const cmT = 0.1 * (1 + P.gammaFilm(50 / P.LMH / kS, 1));
const plan2 = P.planCentrivaLab({ cart: "94-25", rmax: rot.rmax, maxRpm: rot.maxRpm, V0: 22.5, Vp: 0.0001, T: 25,
  solvent: { rho: 867, nu: 6.81e-7 }, solute: { D: 8e-10, beta: 3.84e-4 }, Rint: 1.0, target: { dPh: 20e5, J: 50 / P.LMH, cm: cmT } });
console.log("worked example: cm target", f(cmT), "→ C0cl", f(plan2.C0, 4), "Γcl", f(plan2.Gamma, 2));
assert(Math.abs(plan2.C0 - 0.085) < 0.005, "C0,cl ≈ 0.085 g/L as in the worked example");

// ── β ≤ 0 is refused ──
const bad = P.planCentrivaLab({ cart: "94-25", rmax: rot.rmax, maxRpm: rot.maxRpm, V0: 22.5, Vp: 0.5, T: 25,
  solvent: SOLV.props("dcm", 25), solute: { D: 8e-10, beta: -1e-4 }, Rint: 1, target: { dPh: 20e5, J: 1e-5, cm: 0.2 } });
assert(bad.error === "beta", "negative beta returns an error");

// ── 6 × 8040 at 50 % recovery: wide c_m range, three-point programme ──
const big = P.simulateModule(Object.assign({}, inp, { module: P.MODULES["8040"], nElements: 6, Q_m3h: 10, Y: 0.5 }));
console.log("6×8040 Y=50 %: cmIn", f(big.cmIn), "cmOut", f(big.cmOut), "ratio", f(big.cmOut / big.cmIn, 2), "Δp", f(big.dpChannel / 1e5, 2), "bar", big.flags);
assert(big.cmOut / big.cmIn > 1.5, "c_m ratio > 1.5 triggers the multi-point programme");
