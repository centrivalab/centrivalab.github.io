/* ──────────────────────────────────────────────────────────────────────────────
   CentrivaLab Batch Planner — calculation engine (no DOM)
   Used by batch.html; runnable in node (tests/batch-engine.test.js). Builds on planner-engine.js.

   A batch membrane operation — concentration, diafiltration (impurity wash-out or solvent
   exchange), or a sequence of such stages — is integrated for a spiral-wound module fed from
   a recirculated tank. Six CentrivaLab tubes are then placed along the batch trajectory, each
   reproducing the wall state (flux and membrane-surface concentrations of product and second
   species) of one batch instant at a shared rotor speed. With measured tube results the batch
   is re-integrated using the membrane's measured response instead of the assumed one.

   Species: P (product, retained) and optionally I — a solute impurity (own rejection) or the
   old solvent in a solvent exchange (rejection 0, no polarisation, mixture properties by
   volume-fraction mixing: ρ linear, ln μ linear).
   Tank state: V [m³], cP, cI [kg/m³], phi (volume fraction of old solvent, solvent exchange only).
   Module: 1-D plug flow along the feed channel as in planner-engine.js, two solutes.
   Cartridge: O-CMF buoyancy model (product drives buoyancy; I treated as a tracer).
   ────────────────────────────────────────────────────────────────────────────── */
(function (root, factory) {
  if (typeof module !== "undefined" && module.exports) module.exports = factory(require("./planner-engine.js"), require("./ocmf-geometry.js"));
  else root.BATCH = factory(root.PLANNER, root.OCMF);
})(typeof window !== "undefined" ? window : globalThis, function (P, OCMF) {
  "use strict";
  const LMH = P.LMH;

  // ── solvent mixture and feed properties at a tank state ──
  // sys = { solvent: {rho, mu}, oldSolvent: {rho, mu} | null, kmu (relative viscosity per kg/m³ of P), T }
  function feedProps(sys, st) {
    let rho = sys.solvent.rho, mu = sys.solvent.mu;
    if (sys.oldSolvent && st.phi > 0) {
      rho = st.phi * sys.oldSolvent.rho + (1 - st.phi) * sys.solvent.rho;
      mu = Math.exp(st.phi * Math.log(sys.oldSolvent.mu) + (1 - st.phi) * Math.log(sys.solvent.mu));
    }
    mu *= 1 + (sys.kmu || 0) * st.cP;
    return { rho, mu, nu: mu / rho };
  }
  // species transport in that feed: D (Stokes–Einstein at the feed viscosity, or the user value scaled by μ), β
  function speciesProps(sp, feed, T) {
    const r = Math.cbrt(3 * (sp.M / 1000) / (4 * Math.PI * sp.rhoS * 6.02214076e23));
    const D = sp.D ? sp.D * (sp.muRef / feed.mu) : 1.380649e-23 * (T + 273.15) / (6 * Math.PI * feed.mu * r);
    const beta = (1 - feed.rho / sp.rhoS) / feed.rho;
    return { D, beta };
  }
  const piOf = (c, sp, T, piFactor) => (sp ? P.piVH(c, sp.i || 1, sp.M, T) * (piFactor || 1) : 0);

  // ── module at a tank state: 1-D along the channel, two solutes, permeance Lp [m/s/Pa] ──
  // cfg = { module, spacer, nElements, Q_m3h, dPin_bar, T, Lp, RP, RI, sys, spP, spI, piFactor, cells }
  function moduleState(cfg, st) {
    const N = cfg.cells || 30 * cfg.nElements;
    const Atot = cfg.module.area * cfg.nElements, Ltot = cfg.module.W * cfg.nElements;
    const Ach = cfg.module.area * cfg.spacer.h / (2 * cfg.module.W);
    const feed = feedProps(cfg.sys, st);
    const pP = speciesProps(cfg.spP, feed, cfg.T);
    const hasI = cfg.spI && cfg.spI.kind === "solute" && st.cI > 0;
    const pI = hasI ? speciesProps(cfg.spI, feed, cfg.T) : null;
    const Q0 = cfg.Q_m3h / 3600, dA = Atot / N, dx = Ltot / N;
    let Q = Q0, cbP = st.cP, cbI = st.cI, p = cfg.dPin_bar * 1e5;
    let sumJ = 0, sumJcpP = 0, sumJcpI = 0, sumCmP = 0, sumCmI = 0, sumDP = 0, cmPmax = 0, blocked = false;
    for (let n = 0; n < N; n++) {
      const u = Q / (cfg.spacer.eps * Ach);
      const tP = P.swK(u, cfg.spacer, feed.nu, pP.D);
      const kI = hasI ? P.swK(u, cfg.spacer, feed.nu, pI.D).k : 0;
      // local wall state, both solutes coupled through the osmotic pressure
      let cmP = cbP, cmI = cbI, J = 0, GP = 0, GI = 0;
      for (let it = 0; it < 80; it++) {
        const cpP = (1 - cfg.RP) * cmP, cpI = hasI ? (1 - cfg.RI) * cmI : 0;
        const dpi = piOf(cmP, cfg.spP, cfg.T, cfg.piFactor) - piOf(cpP, cfg.spP, cfg.T, cfg.piFactor)
                  + (hasI ? piOf(cmI, cfg.spI, cfg.T, cfg.piFactor) - piOf(cpI, cfg.spI, cfg.T, cfg.piFactor) : 0);
        J = cfg.Lp * (p - dpi);
        if (J <= 0) { J = 0; blocked = true; break; }
        GP = P.gammaFilm(J / tP.k, cfg.RP); GI = hasI ? P.gammaFilm(J / kI, cfg.RI) : 0;
        const nP = cbP * (1 + GP), nI = cbI * (1 + GI);
        const err = Math.abs(nP - cmP) / (cbP + 1e-12) + (hasI ? Math.abs(nI - cmI) / (cbI + 1e-12) : 0);
        cmP += 0.5 * (nP - cmP); cmI += 0.5 * (nI - cmI);
        if (err < 1e-7) break;
      }
      const cpP = (1 - cfg.RP) * cmP, cpI = hasI ? (1 - cfg.RI) * cmI : (st.cI > 0 ? cbI : 0);   // old solvent: permeates freely
      sumJ += J * dA; sumJcpP += J * cpP * dA; sumJcpI += J * cpI * dA; sumCmP += cmP * dA; sumCmI += cmI * dA; sumDP += p * dA;
      cmPmax = Math.max(cmPmax, cmP);
      const Qn = Q - J * dA;
      cbP = (Q * cbP - J * dA * cpP) / Qn; cbI = (Q * cbI - J * dA * cpI) / Qn; Q = Qn;
      p -= tP.lambda * feed.rho * u * u / (2 * cfg.spacer.dh) * dx;
      if (p <= 0) { blocked = true; break; }
    }
    const Jm = sumJ / Atot;
    return { J: Jm, J_lmh: Jm * LMH, cpP: sumJ > 0 ? sumJcpP / sumJ : 0, cpI: sumJ > 0 ? sumJcpI / sumJ : 0,
             cmP: sumCmP / Atot, cmI: sumCmI / Atot, cmPmax, dPh: sumDP / Atot, Atot, feed, pP, pI, blocked, Ypass: sumJ / Q0 };
  }

  /* simulateBatch
     inp = { V0_L, cP0, cI0, phi0, stages: [{type:"conc", cf | cFinal | VfinalL} | {type:"df", dv | cIfinal | phiFinal}],
             cfg (moduleState cfg), cSat (kg/m³, optional), lpOf?, rpOf? (measured response, functions of cmP) }
     Returns the trajectory and per-stage summary.                                                        */
  function simulateBatch(inp) {
    const cfg = Object.assign({}, inp.cfg);
    const st = { V: inp.V0_L * 1e-3, cP: inp.cP0, cI: inp.cI0 || 0, phi: inp.phi0 || 0 };
    const V0 = st.V, mP0 = st.V * st.cP, mI0 = st.V * st.cI;
    let t = 0, lostP = 0, Vperm = 0, Vadded = 0;
    const traj = [], stages = [];
    const maxSteps = 4000;
    for (let si = 0; si < inp.stages.length; si++) {
      const sg = inp.stages[si];
      const s0 = { t, V: st.V, cP: st.cP, cI: st.cI, phi: st.phi, Vadded: 0, dv: 0 };
      let steps = 0, VaddedStage = 0, lostPstage = 0;
      const done = () => {
        if (sg.type === "conc") {
          if (sg.cf) return st.V <= s0.V / sg.cf;
          if (sg.cFinal) return st.cP >= sg.cFinal;
          if (sg.VfinalL) return st.V <= sg.VfinalL * 1e-3;
        } else {
          if (sg.dv) return VaddedStage / s0.V >= sg.dv;
          if (sg.cIfinal) return (cfg.spI && cfg.spI.kind === "solvent" ? st.phi * cfg.sys.oldSolvent.rho : st.cI) <= sg.cIfinal;
          if (sg.phiFinal) return st.phi <= sg.phiFinal;
        }
        return true;
      };
      while (!done() && steps < maxSteps) {
        if (inp.lpOf) cfg.Lp = inp.lpOf(st.cP);
        if (inp.rpOf) cfg.RP = inp.rpOf(st.cP);
        const m = moduleState(cfg, st);
        if (!(m.J > 0)) { stages.push({ type: sg.type, blocked: true, tStart: s0.t, tEnd: t }); return finish("blocked"); }
        // step size: 1/400 of the initial batch volume of permeate, or the stage remainder
        let dV = V0 / 400;
        if (sg.type === "conc") {
          const Vend = sg.cf ? s0.V / sg.cf : sg.VfinalL ? sg.VfinalL * 1e-3 : 0;
          if (Vend > 0) dV = Math.min(dV, Math.max(st.V - Vend, 1e-9));
        } else if (sg.dv) dV = Math.min(dV, Math.max(sg.dv * s0.V - VaddedStage, 1e-9));
        const dt = dV / (m.J * m.Atot);
        const solvOut = st.phi;                                     // old-solvent fraction in the permeate
        traj.push({ t, stage: si, V: st.V, cP: st.cP, cI: st.cI, phi: st.phi, J_lmh: m.J_lmh, cmP: m.cmP, cmI: m.cmI, cmPmax: m.cmPmax,
                    cpP: m.cpP, cpI: m.cpI, dPh: m.dPh, dv: VaddedStage / s0.V, yieldP: 1 - lostP / mP0, feed: m.feed });
        const mP = st.V * st.cP - m.cpP * dV, mI = st.V * st.cI - m.cpI * dV;
        lostP += m.cpP * dV; lostPstage += m.cpP * dV; Vperm += dV;
        if (sg.type === "conc") { st.V -= dV; }
        else { VaddedStage += dV; Vadded += dV; st.phi = cfg.spI && cfg.spI.kind === "solvent" ? st.phi * (1 - dV / st.V) : 0; }
        st.cP = mP / st.V; st.cI = cfg.spI && cfg.spI.kind === "solvent" ? st.phi * cfg.sys.oldSolvent.rho : mI / st.V;
        t += dt; steps++;
      }
      stages.push({ type: sg.type, tStart: s0.t, tEnd: t, Vstart: s0.V, Vend: st.V, cPstart: s0.cP, cPend: st.cP, cIstart: s0.cI, cIend: st.cI,
                    dv: VaddedStage / s0.V, lostP: lostPstage, incomplete: steps >= maxSteps });
    }
    return finish("ok");

    function finish(status) {
      const m = moduleState(cfg, st);
      traj.push({ t, stage: inp.stages.length - 1, V: st.V, cP: st.cP, cI: st.cI, phi: st.phi, J_lmh: m.J_lmh, cmP: m.cmP, cmI: m.cmI, cmPmax: m.cmPmax,
                  cpP: m.cpP, cpI: m.cpI, dPh: m.dPh, dv: stages.length ? stages[stages.length - 1].dv : 0, yieldP: 1 - lostP / mP0, feed: m.feed, end: true });
      const cmMaxAll = Math.max(...traj.map(q => q.cmPmax));
      return { status, traj, stages, tTotal: t, yieldP: 1 - lostP / mP0, Vfinal: st.V, cPfinal: st.cP, cIfinal: st.cI, phiFinal: st.phi,
               Vperm, Vadded, mP0, mI0, cmPmax: cmMaxAll, satExceeded: inp.cSat ? cmMaxAll > inp.cSat : false,
               tSat: inp.cSat ? (traj.find(q => q.cmPmax > inp.cSat) || {}).t : undefined };
    }
  }

  // ── choose n sampling instants along the trajectory (log-spaced in c_P within concentration
  //    stages, even in diavolumes within diafiltration stages, tubes shared by stage duration) ──
  function sampleStates(res, n, stagesIn) {
    const S = res.stages.length, out = [];
    const durs = res.stages.map(s => s.tEnd - s.tStart), tT = durs.reduce((a, b) => a + b, 0) || 1;
    let alloc = res.stages.map(d => Math.max(1, Math.round(n * (d.tEnd - d.tStart) / tT)));
    while (alloc.reduce((a, b) => a + b, 0) > n) { const i = alloc.indexOf(Math.max(...alloc)); if (alloc[i] > 1) alloc[i]--; else break; }
    while (alloc.reduce((a, b) => a + b, 0) < n) { const i = durs.indexOf(Math.max(...durs.map((d, j) => d / alloc[j]))); alloc[i]++; }
    for (let si = 0; si < S; si++) {
      const pts = res.traj.filter(q => q.stage === si || (q.end && si === S - 1));
      const k = alloc[si], sg = res.stages[si], first = si === 0;
      for (let j = 0; j < k; j++) {
        const fr = k === 1 ? 1 : (first ? j / (k - 1) : (j + 1) / k);       // include the initial state only in the first stage
        let target;
        if (sg.type === "conc") target = pts.reduce((best, q) => Math.abs(Math.log(q.cP) - (Math.log(sg.cPstart) + fr * Math.log(sg.cPend / sg.cPstart))) < Math.abs(Math.log(best.cP) - (Math.log(sg.cPstart) + fr * Math.log(sg.cPend / sg.cPstart))) ? q : best, pts[0]);
        else target = pts.reduce((best, q) => Math.abs(q.dv - fr * sg.dv) < Math.abs(best.dv - fr * sg.dv) ? q : best, pts[0]);
        out.push(Object.assign({ stageIndex: si, stageType: sg.type, frac: fr }, target));
      }
    }
    return out;
  }

  /* planTubes — cartridge runs reproducing the wall state at each sampled batch state, one shared speed.
     tc = { cart, rmax, maxRpm, V0, Vp, cfg (moduleState cfg), states: [tank states], replicate: 1|2 }   */
  function planTubes(tc) {
    const cfg = tc.cfg, T = cfg.T;
    const walls = tc.states.map(st => ({ st, m: moduleState(cfg, st) }));
    // shared speed from the mean target ΔP at the mean feed density
    const dPmean = walls.reduce((s, w) => s + w.m.dPh, 0) / walls.length;
    const rhoMean = walls.reduce((s, w) => s + w.m.feed.rho * (1 + w.m.pP.beta * w.st.cP), 0) / walls.length;
    const g = OCMF.geomRun(tc.cart, tc.rmax, tc.V0, tc.Vp);
    const rpm = OCMF.rpmFromDp(dPmean / 1e5, rhoMean, g.sqEff);
    const tubes = walls.map((w, i) => {
      const { m, st } = w;
      const hasI = cfg.spI && cfg.spI.kind === "solute" && st.cI > 0;
      const plan = P.planCentrivaLab({
        cart: tc.cart, rmax: tc.rmax, maxRpm: tc.maxRpm, V0: tc.V0, Vp: tc.Vp, T, rpm,
        solvent: { rho: m.feed.rho, nu: m.feed.nu }, solute: { D: m.pP.D, beta: m.pP.beta }, Rint: cfg.RP,
        target: { dPh: m.dPh, J: m.J, cm: m.cmP },
        tracer: hasI ? { D: m.pI.D, Rint: cfg.RI, cm: m.cmI } : null
      });
      const C0I = hasI ? plan.tracer.C0 : (cfg.spI && cfg.spI.kind === "solvent" ? st.cI : 0);
      const massG = tc.V0 * 1e-3 * (plan.rhoFeed || m.feed.rho);       // g, feed mass in the tube
      return Object.assign({ index: i, state: st, wall: m, C0I, phi: st.phi, massG, hasI, error: plan.error }, plan);
    });
    // rotor layout: adjacent concentrations in opposite positions (1↔4, 2↔5, 3↔6)
    const order = [...tubes.keys()].sort((a, b) => tubes[a].massG - tubes[b].massG);
    const positions = new Array(tubes.length);
    const slots = [[1, 4], [2, 5], [3, 6]];
    for (let k = 0; k < order.length; k += 2) {
      const pair = slots[Math.floor(k / 2) % 3];
      positions[order[k]] = pair[0]; if (order[k + 1] !== undefined) positions[order[k + 1]] = pair[1];
    }
    tubes.forEach((tb, i) => { tb.position = positions[i]; });
    const dPs = tubes.map(tb => tb.dP_bar), dPspread = (Math.max(...dPs) - Math.min(...dPs)) / (dPmean / 1e5);
    const masses = tubes.map(tb => tb.massG);
    const pairImbalance = Math.max(...slots.map(([a, b]) => {
      const ta = tubes.find(tb => tb.position === a), tb2 = tubes.find(tb => tb.position === b);
      return ta && tb2 ? Math.abs(ta.massG - tb2.massG) : 0;
    }));
    const flags = [];
    if (rpm > tc.maxRpm) flags.push(`Required speed ${Math.round(rpm)} rpm exceeds the rotor limit (${tc.maxRpm} rpm).`);
    if (dPspread > 0.03) flags.push(`ΔP differs by ${(100 * dPspread).toFixed(1)} % between tubes at the shared speed because the feeds have different densities.`);
    if (pairImbalance > 0.1) flags.push(`Opposite tubes differ by up to ${pairImbalance.toFixed(2)} g: check the imbalance tolerance of your centrifuge or trim the fill volumes.`);
    if (tubes.some(tb => tb.error === "beta")) flags.push("Product lighter than the feed solvent in at least one tube (β ≤ 0): the buoyancy model does not apply.");
    return { rpm, dPmean_bar: dPmean / 1e5, rhoMean, tubes, dPspread, pairImbalance, flags, geom: g };
  }

  /* measuredResponse — from tube results, the membrane response at each tube's wall concentration.
     meas[i] = { Vp_mL, t_min, cpP, cpI? }; tubes from planTubes; cfg for osmotic terms.
     Returns per-tube {cm, Lp, Rint, RintI} and interpolating functions lpOf(cm), rpOf(cm) (log-linear, clamped). */
  function measuredResponse(tubes, meas, cfg) {
    const pts = [];
    tubes.forEach((tb, i) => {
      const r = meas[i]; if (!r || !(r.Vp_mL > 0) || !(r.t_min > 0)) return;
      const J = r.Vp_mL * 1e-6 / (P.A_CARTRIDGE * r.t_min * 60);
      // wall concentration at the measured flux (same tube CP model, tube's C0 and speed)
      let Robs = cfg.RP, Gamma = 0, cm = tb.cbMean;
      for (let it = 0; it < 100; it++) {
        const G = Math.max(Robs, 1e-6) * tb.wall.pP.beta * tb.cbMean * J / tb.Ustar;
        const k = P.K2 * tb.Ustar * Math.pow(G, 0.25) * Math.pow(tb.Sc, -0.25);
        const Rint = r.cpP !== undefined && r.cpP >= 0 ? Math.max(0, Math.min(1, 1 - r.cpP / cm)) : cfg.RP;
        Gamma = P.gammaFilm(J / k, Rint); cm = tb.cbMean * (1 + Gamma);
        const Rn = 1 - (1 - Rint) * (1 + Gamma);
        if (Math.abs(Rn - Robs) < 1e-9) break; Robs = 0.5 * (Robs + Rn);
      }
      const Rint = r.cpP !== undefined && r.cpP >= 0 ? Math.max(0, Math.min(1, 1 - r.cpP / cm)) : cfg.RP;
      const cp = (1 - Rint) * cm;
      const dpi = piOf(cm, cfg.spP, cfg.T, cfg.piFactor) - piOf(cp, cfg.spP, cfg.T, cfg.piFactor);
      const Lp = J / Math.max(tb.dP_bar * 1e5 - dpi, 1);
      const RintI = tb.hasI && r.cpI !== undefined && r.cpI >= 0 ? Math.max(0, Math.min(1, 1 - r.cpI / (tb.tracer.cm))) : undefined;
      pts.push({ index: i, J_lmh: J * LMH, cm, Rint, Lp, Lp_lmhbar: Lp * LMH * 1e5, RintI, Gamma });
    });
    pts.sort((a, b) => a.cm - b.cm);
    const interp = key => x => {
      if (!pts.length) return undefined;
      if (pts.length === 1 || x <= pts[0].cm) return pts[0][key];
      if (x >= pts[pts.length - 1].cm) return pts[pts.length - 1][key];
      for (let i = 1; i < pts.length; i++) if (x <= pts[i].cm) {
        const a = pts[i - 1], b = pts[i], f = (Math.log(x) - Math.log(a.cm)) / (Math.log(b.cm) - Math.log(a.cm));
        return a[key] + f * (b[key] - a[key]);
      }
    };
    return { pts, lpOf: interp("Lp"), rpOf: interp("Rint") };
  }


  /* jcCurve — flux and J·c versus product concentration at the module operating point.
     The diafiltration time for N diavolumes is N·m_P/(A·J(c)·c) (Ng, Lundblad & Mitra, Sep. Sci. Technol. 11 (1976) 499):
     it is shortest at the concentration where J·c is largest. Under gel polarisation that is c_lim/e; here J(c) comes
     from the module model (osmotic pressure, viscosity, polarisation) or from measured tube results.        */
  function jcCurve(cfg, cMin, cMax, n, opts) {
    const pts = [];
    for (let i = 0; i < n; i++) {
      const c = Math.exp(Math.log(cMin) + (Math.log(cMax) - Math.log(cMin)) * i / (n - 1));
      const c2 = Object.assign({}, cfg);
      if (opts && opts.lpOf) c2.Lp = opts.lpOf(c);
      if (opts && opts.rpOf) c2.RP = opts.rpOf(c);
      const m = moduleState(c2, { V: 1, cP: c, cI: opts && opts.cI ? opts.cI : 0, phi: opts && opts.phi ? opts.phi : 0 });
      pts.push({ c, J_lmh: m.J_lmh, Jc: m.J_lmh * c, cmP: m.cmP, blocked: m.blocked });
    }
    const valid = pts.filter(q => q.J_lmh > 0);
    const best = valid.reduce((b, q) => (q.Jc > b.Jc ? q : b), valid[0] || pts[0]);
    const atEdge = best && (best === valid[valid.length - 1]);
    return { pts, cOpt: best ? best.c : undefined, JcMax: best ? best.Jc : 0, atEdge };
  }

  /* scanDfConcentration — total time and product yield of the sequence
     [concentrate to c_DF] → [diafilter N diavolumes] → [concentrate to c_final] as a function of c_DF.
     General form of the Ng et al. criterion: includes osmotic and viscosity effects and product loss. */
  function scanDfConcentration(base, dv, cFinal, cMin, cMax, n) {
    const rows = [];
    for (let i = 0; i < n; i++) {
      const cDF = Math.exp(Math.log(cMin) + (Math.log(cMax) - Math.log(cMin)) * i / (n - 1));
      const stages = [];
      if (cDF > base.cP0 * 1.001) stages.push({ type: "conc", cFinal: cDF });
      stages.push({ type: "df", dv });
      if (cFinal > cDF * 1.001) stages.push({ type: "conc", cFinal });
      const r = simulateBatch(Object.assign({}, base, { stages }));
      rows.push({ cDF, tTotal: r.tTotal, yieldP: r.yieldP, status: r.status, cmPmax: r.cmPmax });
    }
    const ok = rows.filter(r => r.status === "ok");
    const best = ok.reduce((b, r) => (r.tTotal < b.tTotal ? r : b), ok[0]);
    return { rows, cOpt: best ? best.cDF : undefined, tMin: best ? best.tTotal : undefined };
  }


  /* cjStates — tube states for a C·J (diafiltration-optimisation) programme: concentrations log-spaced between
     cMin and cMax; in a solvent exchange, half the tubes in the initial solvent and half in the final one
     (the TFF practice of measuring the C·J curve in both buffers).                                  */
  function cjStates(cMin, cMax, n, exchange, cI) {
    const out = [];
    const per = exchange ? Math.max(1, Math.floor(n / 2)) : n;
    const cs = [...Array(per).keys()].map(i => Math.exp(Math.log(cMin) + (Math.log(cMax) - Math.log(cMin)) * (per === 1 ? 1 : i / (per - 1))));
    for (const phi of exchange ? [1, 0] : [0]) for (const c of cs) out.push({ V: 1, cP: c, cI: cI || 0, phi, stageType: "cj", stageIndex: phi });
    return out.slice(0, n);
  }
  // Analytical checks (Millipore TFF brief): product loss and contaminant remaining for constant R
  const analyticLoss = (R, VCF, N) => 1 - Math.exp((R - 1) * (Math.log(VCF) + N));
  const analyticRemaining = (R, N) => Math.exp((R - 1) * N);

  return { feedProps, speciesProps, moduleState, simulateBatch, sampleStates, planTubes, measuredResponse, jcCurve, scanDfConcentration, cjStates, analyticLoss, analyticRemaining, LMH };
});
