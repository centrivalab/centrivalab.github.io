/* ──────────────────────────────────────────────────────────────────────────────
   CentrivaLab Test Planner — calculation engine (no DOM)
   Used by planner.html; runnable in node for tests (see tests/planner-engine.test.js).

   Question answered: a membrane will run in a spiral-wound element (or a vessel of N elements)
   at a given feed flow, recovery and pressure. Which CentrivaLab run (rotor speed, feed
   concentration) exposes the membrane to the same wall conditions — the same flux J_v and the
   same concentration at the membrane surface c_m — so that the permeate measured in the tube
   is the permeate the module would give?

   Spiral-wound side: 1-D plug-flow model along the feed channel.
     mass transfer   Schock & Miquel (1987): Sh = 0.065 Re^0.875 Sc^0.25, Re = u d_h/ν (u interstitial)
     channel Δp      λ = 6.23 Re^−0.3,  dp/dx = λ ρ u² / (2 d_h)
     film model      c_m = c_p + (c_b − c_p) e^{J_v/k},  c_p = (1 − R_int) c_m
     osmotic         van't Hoff, π = i c R T / M   (adequate for dilute OSN feeds; use the
                     Pitzer-based CP Calculator for concentrated aqueous electrolytes)
   O-CMF side (Geraldes et al. 2026):  k = 0.28 U* 𝒢^{1/4} Sc^{−1/4},
     U* = (a_m D)^{1/3},  𝒢 = R_obs β c_b J_v / U*,  a_m = ω² r_m.

   Units inside the engine are SI except where stated: concentrations in kg/m³ (= g/L),
   pressures in Pa, flux in m/s. Public inputs/outputs use bar, LMH, g/L, m³/h for convenience.
   ────────────────────────────────────────────────────────────────────────────── */
(function (root, factory) {
  if (typeof module !== "undefined" && module.exports) module.exports = factory(require("./ocmf-geometry.js"));
  else root.PLANNER = factory(root.OCMF);
})(typeof window !== "undefined" ? window : globalThis, function (OCMF) {
  "use strict";

  const RGAS = 8.314, KB = 1.380649e-23, NA = 6.02214076e23;
  const K2 = 0.28;                          // O-CMF correlation constant (Eq. 26 of the source paper)
  const A_CARTRIDGE = 1.45e-4;              // m², effective membrane area of both CentrivaLab cartridges
  const LMH = 3.6e6;                        // m/s → L m⁻² h⁻¹

  // ── Reference spiral-wound elements (typical 40" constructions; use "custom" for a datasheet) ──
  // area: active membrane area [m²]; W: feed path length of one element [m].
  // The feed-channel cross-section follows from area, spacer thickness h and W:
  //   A_ch = n_leaves · L_leaf · h = (area / 2W) · h
  const MODULES = {
    "2540": { label: "2.5″ × 40″  (≈ 2.5 m²)", area: 2.5,  W: 0.96, qTyp: 0.6 },
    "4040": { label: "4″ × 40″  (≈ 7.9 m²)",   area: 7.9,  W: 0.96, qTyp: 2.0 },
    "8040": { label: "8″ × 40″  (≈ 38 m²)",    area: 38.0, W: 0.96, qTyp: 10.0 },
    "custom": { label: "Custom element", area: 7.9, W: 0.96, qTyp: 2.0 }
  };

  // ── Feed spacers. d_h and ε from Schock & Miquel (1987) Table 1 for the FT30 (28 mil) and Desal UF
  //    (46 mil) nets; the 31 and 34 mil entries are scaled from the FT30 net with the same filament
  //    specific surface (d_h = 4ε / (2/h + (1−ε) S_v), S_v ≈ 8 450 m⁻¹). ──
  const SPACERS = {
    "28mil": { label: "28 mil (0.71 mm) — standard RO/NF (FT30-type)", h: 0.71e-3, dh: 0.95e-3, eps: 0.89 },
    "31mil": { label: "31 mil (0.79 mm)",                              h: 0.79e-3, dh: 1.03e-3, eps: 0.89 },
    "34mil": { label: "34 mil (0.86 mm) — fouling-resistant / OSN",   h: 0.86e-3, dh: 1.09e-3, eps: 0.89 },
    "46mil": { label: "46 mil (1.17 mm) — wide UF/OSN net (Desal UF)", h: 1.17e-3, dh: 1.39e-3, eps: 0.81 },
    "custom": { label: "Custom spacer", h: 0.71e-3, dh: 0.95e-3, eps: 0.89 }
  };

  // ── Solute transport properties in a solvent {rho, mu, nu} ──
  // D: Stokes–Einstein with the molecular radius from M and the solute density (r = (3M/4πρ_s N_A)^{1/3}).
  // β = (∂ρ/∂c)/ρ for a volume-additive mixture: (1 − ρ_solv/ρ_solute)/ρ_solv  [m³/kg = (g/L)⁻¹].
  function soluteProps(sol, M, rhoSolute, T) {
    const r = Math.cbrt(3 * (M / 1000) / (4 * Math.PI * rhoSolute * NA));
    const D = KB * (T + 273.15) / (6 * Math.PI * sol.mu * r);
    const beta = (1 - sol.rho / rhoSolute) / sol.rho;
    return { D, beta, r };
  }

  // Film model: Γ = R_int (E − 1) / (R_int + (1 − R_int) E), E = e^{φ}  (c_m = c_b (1 + Γ))
  function gammaFilm(phi, Rint) { const E = Math.exp(Math.min(phi, 20)); return Rint * (E - 1) / (Rint + (1 - Rint) * E); }
  const piVH = (c, i, M, T) => i * (c / (M / 1000)) * RGAS * (T + 273.15);   // Pa, c in kg/m³

  // Spiral-wound mass transfer and friction at interstitial velocity u
  function swK(u, sp, nu, D) {
    const Re = u * sp.dh / nu, Sc = nu / D;
    return { Re, Sc, k: 0.065 * Math.pow(Re, 0.875) * Math.pow(Sc, 0.25) * D / sp.dh, lambda: 6.23 * Math.pow(Re, -0.3) };
  }

  // Local wall state for given bulk c_b, mass transfer k, permeance Lp [m/s/Pa], hydraulic TMP [Pa]
  function wallState(cb, k, Lp, dPh, Rint, osm) {
    let cm = cb, J = 0;
    for (let it = 0; it < 200; it++) {
      const cp = (1 - Rint) * cm;
      J = Lp * (dPh - (piVH(cm, osm.i, osm.M, osm.T) - piVH(cp, osm.i, osm.M, osm.T)));
      if (J <= 0) return { J: 0, cm: cb, cp: (1 - Rint) * cb, Gamma: 0, phi: 0, blocked: true };
      const Gamma = gammaFilm(J / k, Rint);
      const cmNew = cb * (1 + Gamma);
      const err = Math.abs(cmNew - cm) / (cb + 1e-12);
      cm = cm + 0.5 * (cmNew - cm);
      if (err < 1e-8) return { J, cm, cp: (1 - Rint) * cm, Gamma, phi: J / k };
    }
    return { J, cm, cp: (1 - Rint) * cm, Gamma: cm / cb - 1, phi: J / k, noconv: true };
  }

  /* simulateModule
     inp = { module: {area, W}, spacer: {h, dh, eps}, nElements, Q_m3h, dPin_bar, T,
             solvent: {rho, nu}, solute: {D, beta, M, i}, C0, Rint,
             mode: "recovery" | "permeance", Y (fraction) | Lp_lmhbar, cells }               */
  function simulateModule(inp) {
    const N = inp.cells || 60 * inp.nElements;
    const Atot = inp.module.area * inp.nElements, Ltot = inp.module.W * inp.nElements;
    const Ach = inp.module.area * inp.spacer.h / (2 * inp.module.W);
    const Q0 = inp.Q_m3h / 3600, dA = Atot / N, dx = Ltot / N;
    const osm = { i: inp.solute.i || 1, M: inp.solute.M, T: inp.T };
    const { rho, nu } = inp.solvent, D = inp.solute.D;

    function run(Lp) {
      let Q = Q0, cb = inp.C0, p = inp.dPin_bar * 1e5, Vp = 0, mp = 0;
      let sumCmA = 0, sumCmJ = 0, sumJ = 0, sumDP = 0;
      const prof = [];
      for (let n = 0; n < N; n++) {
        const u = Q / (inp.spacer.eps * Ach);
        const t = swK(u, inp.spacer, nu, D);
        const w = wallState(cb, t.k, Lp, p, inp.Rint, osm);
        prof.push({ x: (n + 0.5) * dx, u, Re: t.Re, k: t.k, cb, cm: w.cm, cp: w.cp, J: w.J, Gamma: w.Gamma, p, dPh: p });
        sumCmA += w.cm * dA; sumCmJ += w.cm * w.J * dA; sumJ += w.J * dA; sumDP += p * dA;
        const Qn = Q - w.J * dA;
        cb = (Q * cb - w.J * dA * w.cp) / Qn;
        Vp += w.J * dA; mp += w.J * dA * w.cp; Q = Qn;
        p -= t.lambda * rho * u * u / (2 * inp.spacer.dh) * dx;
        if (p <= 0 || Q <= 0) break;
      }
      const first = prof[0], last = prof[prof.length - 1];
      return {
        Lp, prof, Y: Vp / Q0, Sc: nu / D,
        Jmean: sumJ / Atot, dPmean: sumDP / Atot, dPin: inp.dPin_bar * 1e5, dPout: last.p, dpChannel: inp.dPin_bar * 1e5 - last.p,
        cmMeanA: sumCmA / Atot, cmMeanJ: sumJ > 0 ? sumCmJ / sumJ : first.cm, cmIn: first.cm, cmOut: last.cm,
        cbOut: cb, cpMix: Vp > 0 ? mp / Vp : 0, uIn: first.u, uOut: last.u, ReIn: first.Re, ReOut: last.Re,
        kIn: first.k, kOut: last.k, GammaIn: first.Gamma, GammaOut: last.Gamma, Ach, Atot, Ltot
      };
    }

    let res;
    if (inp.mode === "permeance") res = run(inp.Lp_lmhbar / LMH / 1e5);
    else {
      // find the permeance that yields the requested recovery (Y is monotonic in Lp)
      let lo = Math.log(1e-3 / LMH / 1e5), hi = Math.log(1e3 / LMH / 1e5);
      for (let it = 0; it < 60; it++) {
        const mid = 0.5 * (lo + hi); res = run(Math.exp(mid));
        if (res.Y < inp.Y) lo = mid; else hi = mid;
      }
      res = run(Math.exp(0.5 * (lo + hi)));
    }
    res.Robs = 1 - res.cpMix / inp.C0;
    res.Lp_lmhbar = res.Lp * LMH * 1e5;
    res.Jmean_lmh = res.Jmean * LMH;
    // representative points for a multi-point programme: inlet, area-median, outlet
    const pick = (i, label) => { const q = res.prof[i]; return { label, x: q.x, dPh: q.dPh, J: q.J, cm: q.cm, cb: q.cb, Gamma: q.Gamma, u: q.u }; };
    res.points = [pick(0, "Inlet"), pick(Math.floor(res.prof.length / 2), "Mid-length"), pick(res.prof.length - 1, "Outlet")];
    res.flags = [];
    if (res.ReOut < 50 || res.ReIn > 1500) res.flags.push("Re outside the range of the Schock & Miquel data (≈ 50–1500).");
    if (res.dpChannel > 1e5 * inp.nElements) res.flags.push("Feed-channel pressure drop exceeds ≈ 1 bar per element.");
    if (res.prof.some(q => q.J <= 0)) res.flags.push("Osmotic pressure blocks permeation in part of the module.");
    const piOut = piVH(res.cmOut, osm.i, osm.M, osm.T);
    if (piOut > 0.05 * res.dPmean) res.flags.push(`Osmotic pressure at the wall reaches ${(piOut / 1e5).toFixed(1)} bar (van't Hoff) — no longer negligible.`);
    return res;
  }

  /* planCentrivaLab
     Finds the CentrivaLab feed concentration and speed that reproduce the wall state (J, c_m) of a target.
     inp = { cart, rmax, maxRpm, V0, Vp, T, solvent: {rho, nu}, solute: {D, beta}, Rint,
             target: {dPh [Pa], J [m/s], cm [kg/m³]} }                                             */
  function planCentrivaLab(inp) {
    const g = OCMF.geomRun(inp.cart, inp.rmax, inp.V0, inp.Vp);
    const rhoFeed = inp.solvent.rho;                    // dilute: the β·c term is < 0.1 %
    const dP_bar = inp.target.dPh / 1e5;
    const rpm = OCMF.rpmFromDp(dP_bar, rhoFeed, g.sqEff);
    const a_m = OCMF.amOf(rpm, g.r_m);
    const { D, beta } = inp.solute, nu = inp.solvent.nu, J = inp.target.J, Rint = inp.Rint;
    const y = inp.Vp / inp.V0;
    const Ustar = Math.cbrt(a_m * D), Sc = nu / D;

    // wall state of a tube run starting at feed concentration C0 (run-mean bulk concentration)
    function tube(C0) {
      let Robs = Rint, out = null;
      for (let it = 0; it < 100; it++) {
        const cbMean = 0.5 * C0 * (1 + (1 - y * (1 - Robs)) / (1 - y));
        const G = Math.max(Robs, 1e-6) * beta * cbMean * J / Ustar;
        const k = K2 * Ustar * Math.pow(G, 0.25) * Math.pow(Sc, -0.25);
        const Gamma = gammaFilm(J / k, Rint);
        const Rn = 1 - (1 - Rint) * (1 + Gamma);
        out = { C0, cbMean, G, k, phi: J / k, Gamma, Robs: Rn, cm: cbMean * (1 + Gamma), cp: (1 - Rint) * cbMean * (1 + Gamma) };
        if (Math.abs(Rn - Robs) < 1e-10) break;
        Robs = 0.5 * (Robs + Rn);
      }
      return out;
    }
    if (!(beta > 0)) return { rpm, a_m, error: "beta" };
    // c_m is monotonic in C0: bisection in log C0
    let lo = Math.log(1e-5), hi = Math.log(1e4), t = null;
    for (let it = 0; it < 80; it++) {
      const mid = 0.5 * (lo + hi); t = tube(Math.exp(mid));
      if (t.cm < inp.target.cm) lo = mid; else hi = mid;
    }
    t = tube(Math.exp(0.5 * (lo + hi)));
    const res = Object.assign({
      rpm, a_m, ag: a_m / OCMF.G, dP_bar, over: rpm > inp.maxRpm, maxRpm: inp.maxRpm,
      dP0: OCMF.dpFromRpm(rpm, rhoFeed, g.sq0), dP1: OCMF.dpFromRpm(rpm, rhoFeed, g.sq1),
      Ustar, Sc, J_lmh: J * LMH, tRun_min: inp.Vp * 1e-6 / (J * A_CARTRIDGE) / 60, y, geom: g
    }, t);
    res.flags = [];
    if (res.over) res.flags.push(`Required speed ${Math.round(rpm)} rpm exceeds the rotor limit (${inp.maxRpm} rpm).`);
    if (res.phi > 1.5) res.flags.push("J_v/k > 1.5 in the tube: beyond the validated range of the O-CMF correlation.");
    if (res.G < 7e-9 || res.G > 5e-6) res.flags.push("Buoyancy parameter G outside the range where the O-CMF correlation was validated (aqueous salts): treat the cartridge prediction as an extrapolation.");
    if (res.Sc < 700 || res.Sc > 2200) res.flags.push("Sc outside the validated range 700–2200 (aqueous salts): treat the O-CMF prediction as an extrapolation.");
    return res;
  }

  return { MODULES, SPACERS, soluteProps, gammaFilm, piVH, swK, wallState, simulateModule, planCentrivaLab, A_CARTRIDGE, LMH, K2 };
});
