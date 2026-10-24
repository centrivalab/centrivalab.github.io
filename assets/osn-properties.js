/* OSN property layer for the CentrivaLab OSN simulator.
   Solvents, solvent + water mixtures, and solute property estimates (density/β, viscosity, diffusivity, osmotic pressure).
   Every number here is documented, with its source and its validation status, in guide-osn-properties.html.
   No DOM access: usable from Node for testing. */
(function (root) {
  const R = 8.314462618, NA = 6.02214076e23, KB = 1.380649e-23;

  // ρ in kg/m³ at 25 °C (Choupo/ChemSep v8.3; acetone, ethyl acetate, DMF: median of NIST ThermoML values);
  // μ(T) [Pa·s]: kind "eq101" ln μ = A + B/T + C ln T + D T^E ; "andrade" ln μ = A + B/T ; "vogel" ln μ = A + B/(T − C).
  // φ = Wilke–Chang association factor. src: where the viscosity law comes from.
  const SOLVENTS = {
    water:        { label: "Water (reference)",       M: 18.015, rho: 997.0,  phi: 2.6, mu: { kind: "vogel",   A: -10.5366, B: 541.4, C: 144.0 }, src: "Choupo (Vogel)", status: "ok" },
    acetonitrile: { label: "Acetonitrile",            M: 41.052, rho: 776.6,  phi: 1.0, mu: { kind: "eq101",   A: -31.531, B: 1522.1, C: 3.3306, D: -5.9061e-6, E: 2 }, src: "Choupo / ChemSep Eq.101", status: "ok" },
    methanol:     { label: "Methanol",                M: 32.042, rho: 787.3,  phi: 1.9, mu: { kind: "andrade", A: -11.5712, B: 1207.4 }, src: "Choupo (Andrade)", status: "ok" },
    ethanol:      { label: "Ethanol",                 M: 46.069, rho: 787.5,  phi: 1.5, mu: { kind: "andrade", A: -12.6467, B: 1736.8 }, src: "Choupo (Andrade)", status: "ok" },
    isopropanol:  { label: "Isopropanol",             M: 60.096, rho: 782.2,  phi: 1.0, mu: { kind: "eq101",   A: -7.4407, B: 2259.7, C: -1.1149, D: 2.963e-7, E: 2 }, src: "Choupo / ChemSep Eq.101", status: "ok" },
    acetone:      { label: "Acetone",                 M: 58.080, rho: 784.6,  phi: 1.0, mu: { kind: "andrade", A: -10.970, B: 861.2 },  src: "fit to 31 NIST points (not ChemSep)", status: "fit" },
    ethylacetate: { label: "Ethyl acetate",           M: 88.106, rho: 894.4,  phi: 1.0, mu: { kind: "andrade", A: -11.080, B: 989.5 },  src: "fit to 27 NIST points (not ChemSep)", status: "fit" },
    toluene:      { label: "Toluene",                 M: 92.141, rho: 861.9,  phi: 1.0, mu: { kind: "andrade", A: -11.0586, B: 1062.3 }, src: "Choupo (Andrade)", status: "ok" },
    thf:          { label: "Tetrahydrofuran",         M: 72.107, rho: 881.4,  phi: 1.0, mu: { kind: "eq101",   A: -10.047, B: 896.61, C: -0.11864, D: 2.5625e-7, E: 2 }, src: "Choupo / ChemSep Eq.101", status: "ok" },
    dmf:          { label: "N,N-Dimethylformamide",   M: 73.095, rho: 943.9,  phi: 1.0, mu: { kind: "andrade", A: -10.719, B: 1072.7 }, src: "fit to 146 NIST points (not in Choupo)", status: "fit" },
    dmso:         { label: "Dimethyl sulfoxide",      M: 78.133, rho: 1096.1, phi: 1.0, mu: { kind: "eq101",   A: 33.492, B: 223.29, C: -7.3227, D: 1.4294e-5, E: 2 }, src: "Choupo / ChemSep Eq.101", status: "ok" }
  };

  // Organic + water viscosity: ln μ = x1 ln μ1 + x2 ln μ2 + x1 x2 Σ_k (a_k + b_k/T)(x1 − x2)^k , x1 = mole fraction of the organic.
  // Fitted to NIST ThermoML binary data (ambient pressure); Tmin/Tmax = data range; nArt = number of source articles.
  const MIX = {
    acetonitrile: { c: [-0.150761, 42.1188, 4.3087, -1705.99, 6.08107, -1563.05], Tmin: 298.15, Tmax: 318.15, nArt: 1, medErr: 1.4, note: "single source article — limited" },
    methanol:     { c: [-4.29177, 2031.61, 6.56434, -2469.03, -7.54496, 2707.81], Tmin: 273.15, Tmax: 323.15, nArt: 3, medErr: 1.6 },
    ethanol:      { c: [-5.06711, 2285.02, 11.3803, -4429.45, -6.90746, 3144.01], Tmin: 268.1,  Tmax: 328.15, nArt: 6, medErr: 3.7 },
    dmf:          { c: [-5.03251, 2438.85, 7.06585, -3279.48, 2.34879, -110.681],  Tmin: 293.15, Tmax: 353.15, nArt: 4, medErr: 1.7 }
  };

  function muPure(key, T) {
    const m = SOLVENTS[key].mu;
    if (m.kind === "eq101")   return Math.exp(m.A + m.B / T + m.C * Math.log(T) + m.D * Math.pow(T, m.E));
    if (m.kind === "andrade") return Math.exp(m.A + m.B / T);
    return Math.exp(m.A + m.B / (T - m.C));
  }

  /* Solvent (pure or organic + water). spec = { solvent, T (K), xOrg (mole fraction of the organic, 0..1; ignored if pure) , mixWater (bool) }.
     Returns ρ [kg/m³], μ [Pa·s], ν [m²/s], M [g/mol], phi, warnings[], tags[]. */
  function solventProps(spec) {
    const T = spec.T, s = SOLVENTS[spec.solvent], warn = [], tags = [];
    if (!spec.mixWater || spec.solvent === "water" || !MIX[spec.solvent]) {
      const mu = muPure(spec.solvent, T);
      if (s.status === "fit") tags.push(`viscosity law fitted by CentrivaLab to NIST data (${s.src})`);
      if (T < 283 || T > 343) warn.push("temperature outside 10–70 °C: viscosity law is extrapolated");
      return { rho: s.rho, mu, nu: mu / s.rho, M: s.M, phi: s.phi, warn, tags };
    }
    const x1 = Math.min(1, Math.max(0, spec.xOrg)), x2 = 1 - x1, w = SOLVENTS.water, mx = MIX[spec.solvent];
    const mu1 = muPure(spec.solvent, T), mu2 = muPure("water", T), d = x1 - x2, b = x1 * x2;
    let corr = 0;
    for (let k = 0; k < 3; k++) corr += (mx.c[2 * k] + mx.c[2 * k + 1] / T) * Math.pow(d, k);
    const mu = Math.exp(x1 * Math.log(mu1) + x2 * Math.log(mu2) + b * corr);
    const Mm = x1 * s.M + x2 * w.M, wm = x1 * s.M / Mm;                      // mass fraction of the organic
    const rho = 1 / (wm / s.rho + (1 - wm) / w.rho);                         // ideal volume additivity (≤ 1–3 % vs NIST)
    if (T < mx.Tmin - 2 || T > mx.Tmax + 2) warn.push(`temperature outside the ${mx.Tmin.toFixed(0)}–${mx.Tmax.toFixed(0)} K range of the mixture-viscosity data: extrapolated`);
    tags.push(`mixture viscosity: Redlich–Kister fit to ${mx.nArt} NIST article${mx.nArt > 1 ? "s" : ""}, median error ${mx.medErr}%${mx.note ? " — " + mx.note : ""}`);
    return { rho, mu, nu: mu / rho, M: Mm, phi: x1 * s.phi + x2 * w.phi, warn, tags };
  }

  /* Solute. spec = { cls: "small"|"globular"|"peg"|"custom", M (g/mol), vbar (cm³/g|null), D0 (m²/s|null), etaInt (mL/g|null), kH, A2 (mol·m³/kg²) } */
  function pegVbar(M) { const n = Math.max(1, (M - 18.015) / 44.053); return (37.0 * n + 18.0) / M; }   // cm³/g (see guide §3)
  function soluteProps(spec, solv, T) {
    const warn = [], tags = [], M = spec.M, cls = spec.cls;
    let vbar = spec.vbar, vtag = "user value";
    if (!(vbar > 0)) {
      if (cls === "peg") { vbar = pegVbar(M); vtag = "group model, validated in water only"; }
      else { vbar = 0.73; vtag = "default 0.73 cm³/g — estimate";
             tags.push(cls === "small" ? "v̄ default 0.73 cm³/g (organic molecules span roughly 0.65–0.9 cm³/g): estimate — enter a measured value if you have one"
                                       : "v̄ default 0.73 cm³/g (typical protein value; short peptides in water 0.58–0.74): estimate"); }
    }
    const beta = (1 - solv.rho * vbar * 1e-3) / solv.rho;                       // m³/kg
    // Hydrodynamic radius [m]
    let D0 = spec.D0, Dtag = "user value";
    if (!(D0 > 0)) {
      if (cls === "small") {
        const Vb = 1.35 * M * vbar;                                              // cm³/mol, calibrated on the two anthraquinones (guide §6)
        D0 = 7.4e-8 * Math.sqrt(solv.phi * solv.M) * T / (solv.mu * 1e3 * Math.pow(Vb, 0.6)) * 1e-4;
        Dtag = "Wilke–Chang"; tags.push("D: Wilke–Chang with V_b ≈ 1.35·M·v̄; validated ±6 % (small aromatics) to −21 % (API-sized dyes) in methanol/ethanol only");
        if (M > 600) warn.push("Wilke–Chang is not meant for M > ~600 g/mol: use the globular or polymer option");
      } else {
        const Rh = cls === "peg" ? 0.0145e-9 * Math.pow(M, 0.571)
                                 : 1.25 * Math.cbrt(3 * M * vbar * 1e-6 / (4 * Math.PI * NA));    // compact sphere × 1.25 hydration
        D0 = KB * T / (6 * Math.PI * solv.mu * Rh); Dtag = "Stokes–Einstein";
        tags.push(cls === "peg" ? "D: Stokes–Einstein with R_h = 0.0145·M^0.571 nm (water power law; constants not confirmed in the original; not validated in organic solvents)"
                                : "D: Stokes–Einstein, compact sphere with 25 % hydration shell (checked on lysozyme in water only)");
      }
    }
    let eta = spec.etaInt;                                                        // mL/g
    if (eta === null || eta === undefined || !(eta >= 0)) {
      if (cls === "peg") { eta = 0.0125 * Math.pow(M, 0.78); tags.push("[η] = 0.0125·M^0.78 mL/g: constants not confirmed in the original source, but they reproduce the measured viscosity of PEG 35 kDa in water within 1 %; assumed also in the organic solvent"); }
      else { eta = 2.5 * vbar * 1.95; tags.push("[η] from Einstein sphere with hydration (≈ 3.6 mL/g): estimate"); }
    }
    const kH = spec.kH > 0 ? spec.kH : (cls === "peg" ? 0.4 : 1.0);
    const etaM3 = eta * 1e-3;                                                     // m³/kg
    const nu0 = solv.nu;
    const relVisc = c => 1 + etaM3 * c + kH * etaM3 * etaM3 * c * c;
    const A2 = spec.A2 || 0, macro = (cls === "peg" || cls === "globular");
    // D(c): for macromolecules the collective (mutual) diffusivity rises with concentration through the thermodynamic factor while the
    // friction rises with the viscosity; for PEG 35 kDa in water (measured, Pereira 2026) the net effect is D/D0 ≈ 1 + 0.15 [η] c.
    // The same form is used for every macromolecule (calibrated on that single case); small molecules keep a constant D.
    const Dc = c => macro ? D0 * (1 + 0.15 * etaM3 * c) : D0;
    // Osmotic pressure [Pa], c in kg/m³:
    //  PEG coil: semi-dilute scaling  π = (RT/M) c [1 + 0.73 (c[η])^1.3]  (exponent from des Cloizeaux; 0.73 calibrated on PEG 35 kDa/water)
    //  globular: Carnahan–Starling hard spheres with the hydrated volume fraction φ = 1.95 v̄ c
    //  small:    ideal. A user A2 (mol·m³/kg²) adds the virial term in all cases.
    const piFun = c => {
      const base = c * R * T / ((spec.Mn > 0 ? spec.Mn : M) * 1e-3);
      let z = 1;
      if (cls === "peg") z = 1 + 0.73 * Math.pow(etaM3 * c, 1.3);
      else if (cls === "globular") { const f = Math.min(0.62, 1.95 * vbar * 1e-3 * c); z = (1 + f + f * f - f * f * f) / Math.pow(1 - f, 3); }
      return base * z + R * T * A2 * c * c;
    };
    if (cls === "peg") tags.push("π: semi-dilute scaling calibrated on PEG 35 kDa in water; the measured law of that batch implies a lower effective M_n (≈ 12 kDa), so π at low c is under-estimated");
    if (cls === "globular") tags.push("π: Carnahan–Starling hard spheres with a hydration shell; no charge or association effects");
    return {
      vbar, vtag, beta, D0, Dtag, eta, kH, M, A2,
      nu: c => nu0 * relVisc(c),
      D: Dc,
      pi: piFun,
      warn, tags
    };
  }

  const api = { SOLVENTS, MIX, muPure, solventProps, soluteProps, pegVbar };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.OSNP = api;
})(typeof window !== "undefined" ? window : this);
