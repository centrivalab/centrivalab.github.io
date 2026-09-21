// node tests/batch-engine.test.js — sanity checks of the Batch Planner engine
const SOLV = require("../assets/solvents.js");
const P = require("../assets/planner-engine.js");
const B = require("../assets/batch-engine.js");
const OCMF = require("../assets/ocmf-geometry.js");
const assert = (c, m) => { if (!c) { console.error("FAIL:", m); process.exitCode = 1; } else console.log("ok  ", m); };
const f = (x, d = 3) => Number(x).toFixed(d);

const tol = SOLV.props("toluene", 25);
const spP = { M: 500, rhoS: 1300, i: 1 };
const spI = { kind: "solute", M: 150, rhoS: 1100, i: 1 };
const cfg = { module: P.MODULES["4040"], spacer: P.SPACERS["28mil"], nElements: 1, Q_m3h: 2, dPin_bar: 20, T: 25,
  Lp: 1.5 / P.LMH / 1e5, RP: 0.98, RI: 0.2, sys: { solvent: tol, oldSolvent: null, kmu: 0 }, spP, spI, piFactor: 1 };

// ── module state at a tank state ──
const m = B.moduleState(cfg, { V: 0.1, cP: 10, cI: 1, phi: 0 });
console.log("module @10 g/L: J", f(m.J_lmh, 1), "LMH cmP", f(m.cmP, 2), "cpP", f(m.cpP, 3), "cmI", f(m.cmI, 3), "cpI", f(m.cpI, 3));
assert(m.J_lmh > 0 && m.cmP > 10 && m.cpP < 1, "two-solute module state is physical");
assert(m.cmI > 1 && m.cpI > 0.5 * m.cmI, "low-rejection impurity: little polarisation, permeates");

// ── concentration stage, R=1: c = c0·CF, yield 1 ──
const cfgR1 = Object.assign({}, cfg, { RP: 1, spI: null });
const r1 = B.simulateBatch({ V0_L: 100, cP0: 5, cI0: 0, stages: [{ type: "conc", cf: 5 }], cfg: cfgR1 });
console.log("conc ×5, R=1: t", f(r1.tTotal / 3600, 2), "h, cP", f(r1.cPfinal, 2), "yield", f(r1.yieldP, 4));
assert(Math.abs(r1.cPfinal - 25) < 0.05, "c_final = c0·CF for R = 1");
assert(Math.abs(r1.yieldP - 1) < 1e-6, "no product loss for R = 1");
assert(Math.abs(r1.Vfinal - 0.02) < 1e-5, "V_final = V0/CF");

// ── diafiltration: impurity decays as exp(−(1−R_I)·DV) (R_I = 0.2 → factor e^−0.8·DV) ──
const r2 = B.simulateBatch({ V0_L: 100, cP0: 5, cI0: 1, stages: [{ type: "df", dv: 3 }], cfg: Object.assign({}, cfg, { RP: 1 }) });
const expI = 1 * Math.exp(-(1 - 0.2 * 1.0) * 3);   // ignoring polarisation of I (small): (1−R_obs,I)·DV
console.log("DF 3 DV: cI", f(r2.cIfinal, 4), "expected ≈", f(expI, 4), "cP", f(r2.cPfinal, 3), "V", f(r2.Vfinal, 4));
assert(Math.abs(r2.cIfinal / expI - 1) < 0.15, "impurity wash-out ≈ exp(−(1−R_I)·DV)");
assert(Math.abs(r2.Vfinal - 0.1) < 1e-6 && Math.abs(r2.cPfinal - 5) < 1e-6, "constant volume, product retained (R=1)");

// ── sequence with losses: yield < 1 when R_P < 1 ──
const r3 = B.simulateBatch({ V0_L: 100, cP0: 5, cI0: 1, stages: [{ type: "conc", cf: 4 }, { type: "df", dv: 4 }, { type: "conc", cf: 2 }], cfg, cSat: 60 });
console.log("seq: t", f(r3.tTotal / 3600, 2), "h, cP", f(r3.cPfinal, 2), "cI", f(r3.cIfinal, 4), "yield", f(r3.yieldP, 3), "cm max", f(r3.cmPmax, 1), "sat?", r3.satExceeded, "stages", r3.stages.length);
assert(r3.yieldP < 1 && r3.yieldP > 0.8, "product loss with R_P = 0.98 is a few %");
assert(r3.stages.length === 3 && r3.stages.every(s => !s.incomplete), "three stages completed");
assert(Math.abs(r3.Vfinal - 0.0125) < 1e-5, "V_final = V0/(4·2)");
// mass balance: product in tank + lost = initial
const massOk = Math.abs(r3.Vfinal * r3.cPfinal + (1 - r3.yieldP) * r3.mP0 - r3.mP0) / r3.mP0 < 1e-6;
assert(massOk, "product mass balance closes");

// ── solvent exchange: old solvent fraction decays as exp(−DV) ──
const cfgX = Object.assign({}, cfg, { spI: { kind: "solvent" }, sys: { solvent: SOLV.props("ethanol", 25), oldSolvent: tol, kmu: 0 } });
const r4 = B.simulateBatch({ V0_L: 50, cP0: 5, cI0: 0, phi0: 1, stages: [{ type: "df", dv: 3 }], cfg: cfgX });
console.log("solvent swap 3 DV: phi", f(r4.phiFinal, 4), "expected", f(Math.exp(-3), 4));
assert(Math.abs(r4.phiFinal - Math.exp(-3)) < 0.01, "old solvent fraction = e^−DV");

// ── sampling and tube plan ──
const states = B.sampleStates(r3, 6);
console.log("sampled:", states.map(s => `${s.stageType}${s.stageIndex} cP=${f(s.cP, 1)} cI=${f(s.cI, 3)}`).join(" | "));
assert(states.length === 6, "six tube states");
assert(new Set(states.map(s => s.stageIndex)).size === 3, "all three stages sampled");
const rot = OCMF.rotor("94-25", "sigma12159");
const plan = B.planTubes({ cart: "94-25", rmax: rot.rmax, maxRpm: rot.maxRpm, V0: 22.5, Vp: 0.5, cfg, states });
console.log("tubes: rpm", f(plan.rpm, 0), "ΔP spread", f(100 * plan.dPspread, 2), "% pair imbalance", f(plan.pairImbalance, 3), "g", plan.flags);
plan.tubes.forEach(tb => console.log(`  pos ${tb.position} cP_batch ${f(tb.state.cP, 2)} → C0P ${f(tb.C0, 3)} C0I ${f(tb.C0I, 4)} cm ${f(tb.cm, 2)} (target ${f(tb.wall.cmP, 2)}) ΔP ${f(tb.dP_bar, 2)} m ${f(tb.massG, 2)} g`));
assert(plan.tubes.every(tb => Math.abs(tb.cm - tb.wall.cmP) / tb.wall.cmP < 1e-5), "every tube reproduces the module wall concentration");
assert(plan.tubes.every(tb => tb.position >= 1 && tb.position <= 6) && new Set(plan.tubes.map(t => t.position)).size === 6, "six distinct rotor positions");
assert(plan.tubes.every(tb => !tb.hasI || Math.abs(tb.tracer.cm - tb.wall.cmI) < 1e-9), "impurity wall concentration reproduced");

// ── measured response closes the loop ──
const meas = plan.tubes.map(tb => ({ Vp_mL: 0.5, t_min: tb.tRun_min, cpP: tb.cp, cpI: tb.hasI ? tb.tracer.cp : undefined }));
const mr = B.measuredResponse(plan.tubes, meas, cfg);
console.log("measured Lp", mr.pts.map(q => f(q.Lp_lmhbar, 2)).join(","), "Rint", mr.pts.map(q => f(q.Rint, 3)).join(","));
assert(mr.pts.every(q => Math.abs(q.Lp_lmhbar - 1.5) < 0.05 && Math.abs(q.Rint - 0.98) < 0.005), "entering the predicted results recovers Lp and R_int");
const r3b = B.simulateBatch({ V0_L: 100, cP0: 5, cI0: 1, stages: [{ type: "conc", cf: 4 }, { type: "df", dv: 4 }, { type: "conc", cf: 2 }], cfg, lpOf: mr.lpOf, rpOf: mr.rpOf });
assert(Math.abs(r3b.tTotal / r3.tTotal - 1) < 0.02, "re-integration with recovered response reproduces the batch time");

// ── J·c curve and DF optimum (Ng, Lundblad & Mitra) ──
const jc = B.jcCurve(cfg, 1, 200, 40);
console.log("J·c optimum at c ≈", f(jc.cOpt, 1), "g/L; J·c max", f(jc.JcMax, 0));
assert(jc.cOpt > 1 && jc.cOpt < 200, "J·c has an interior maximum with osmotic pressure");
const scan = B.scanDfConcentration({ V0_L: 100, cP0: 5, cI0: 1, cfg }, 4, 50, 5, 50, 12);
console.log("time scan: optimum c_DF ≈", f(scan.cOpt, 1), "g/L, t_min", f(scan.tMin / 3600, 2), "h; J·c optimum", f(jc.cOpt, 1));
assert(scan.cOpt > 0, "time scan finds an optimum");

// ── C·J programme states and analytical checks ──
const cj = B.cjStates(5, 100, 6, true, 0);
assert(cj.length === 6 && cj.filter(s => s.phi === 1).length === 3, "3 × 2 solvents for the exchange C·J programme");
assert(Math.abs(B.analyticLoss(0.98, 8, 4) - (1 - Math.exp(-0.02 * (Math.log(8) + 4)))) < 1e-12, "analytical product loss formula");
console.log("analytic loss R=0.98, VCF 8, N 4:", f(100 * B.analyticLoss(0.98, 8, 4), 2), "% ; simulated:", f(100 * (1 - r3.yieldP), 2), "%");
assert(Math.abs(B.analyticLoss(0.98, 8, 4) - (1 - r3.yieldP)) < 0.05, "simulated loss close to the constant-R analytical value");
