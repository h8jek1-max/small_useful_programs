/*
 * klima_agro.js – výpočetní vrstva agroklimatické stránky (bez DOM, bez knihoven; prohlížeč i Node).
 *
 * Vstup: R = Float32Array(NY*36*NF) z KlimaDekoder.dekodujBunku, roc = KlimaDekoder.rocni(...), jih = bool.
 *   const A  = KlimaAgro.rocniUkazatele(K, R, roc, jih, maskaDekad);   // {id: Float32Array(NY)}  podle PARAMETRY
 *   const Ku = KlimaAgro.klimaRoku(K, R, roc, jih);                    // podklad pro potenciál (po letech)
 *   const P  = KlimaAgro.potencial(druhy, Ku, rokOd, rokDo, opt);       // {podil, tvrdy, vahaCelkem, n, proDruh}
 * Poznámky k přesnosti: součty počtů za celý rok, roční TG a roční TXx jsou z databáze přesné; hodnoty po dekádách
 * a odvozené veličiny (úhrn srážek = RR1 × SDII, vegetační sezóna, stupňodny) jsou přibližné.
 */
(function (root) {
  'use strict';

  const DEK_DNU = 365 / 36;                    // průměrná délka dekády ve dnech
  const COUNTS = ['fd', 'tr', 'rr1', 'r20', 'su20', 'su25', 'su30', 'su35'];

  // ---------- katalog parametrů (pořadí = pořadí skupin v bublině) ----------
/*  const SKUPINY = [
    { id: 'teplota', nazev: 'Teploty' },
    { id: 'dny', nazev: 'Počty dní' },
    { id: 'sezona', nazev: 'Sezóna' },
    { id: 'mraz', nazev: 'Termíny mrazů' },
    { id: 'srazky', nazev: 'Srážky' },
    { id: 'potencial', nazev: 'Klimatický zemědělský potenciál' }
  ];*/
  const SKUPINY = [
    { id: 'teplota', nazev: 'Vyber všechno z teplot' },
    { id: 'dny', nazev: 'Vyber všechno z vlastností dní' },
    { id: 'sezona', nazev: 'Vyber všechno z délky sezón' },
    { id: 'mraz', nazev: 'Vyber všechno z termínů mrazu' },
    { id: 'srazky', nazev: 'Vyber všechno z parametrů srážek' },
    { id: 'potencial', nazev: 'Klimatický zemědělský potenciál' }
  ];
  const PARAMETRY = [
    { id: 'tnn', skupina: 'teplota', nazev: 'Nejnižší teplota', jedn: '°C', agreg: 'min', presne: false, barva: '#1565c0', riziko: 'pod' },
    { id: 'tn', skupina: 'teplota', nazev: 'Průměrná min. teplota', jedn: '°C', agreg: 'prumer', presne: false, barva: '#0288d1' },
    { id: 'tg', skupina: 'teplota', nazev: 'Průměrná teplota', jedn: '°C', agreg: 'prumer', presne: true, barva: '#ef6c00' },
    { id: 'tx', skupina: 'teplota', nazev: 'Průměrná max. teplota', jedn: '°C', agreg: 'prumer', presne: false, barva: '#e53935' },
    { id: 'txx', skupina: 'teplota', nazev: 'Nejvyšší teplota', jedn: '°C', agreg: 'max', presne: true, barva: '#b71c1c', riziko: 'nad' },
    { id: 'fd', skupina: 'dny', nazev: 'Mrazové dny', jedn: 'dny', agreg: 'soucet', presne: true, barva: '#4fc3f7' },
    { id: 'su20', skupina: 'dny', nazev: 'Dny nad 20 °C', jedn: 'dny', agreg: 'soucet', presne: true, barva: '#fdd835' },
    { id: 'su25', skupina: 'dny', nazev: 'Dny nad 25 °C', jedn: 'dny', agreg: 'soucet', presne: true, barva: '#fb8c00' },
    { id: 'su30', skupina: 'dny', nazev: 'Dny nad 30 °C', jedn: 'dny', agreg: 'soucet', presne: true, barva: '#e53935', riziko: 'nad' },
    { id: 'su35', skupina: 'dny', nazev: 'Dny nad 35 °C', jedn: 'dny', agreg: 'soucet', presne: true, barva: '#7f1d1d', riziko: 'nad' },
    { id: 'tr', skupina: 'dny', nazev: 'Noci nad 20 °C', jedn: 'dny', agreg: 'soucet', presne: true, barva: '#8e24aa' },
    { id: 'rr1', skupina: 'dny', nazev: 'Deštivé dny', jedn: 'dny', agreg: 'soucet', presne: true, barva: '#7cb342' },
    { id: 'r20', skupina: 'dny', nazev: 'Dny s deštěm nad 20 mm', jedn: 'dny', agreg: 'soucet', presne: true, barva: '#1b5e20' },
    { id: 'bezmrazi', skupina: 'sezona', nazev: 'Bezmrazové období', jedn: 'dny', agreg: 'zvl', presne: true, barva: '#26a69a', riziko: 'pod' },
    { id: 'vegSezona', skupina: 'sezona', nazev: 'Vegetační sezóna (průměr od 5 °C)', jedn: 'dny', agreg: 'zvl', presne: false, barva: '#43a047', riziko: 'pod' },
    { id: 'jarniMraz', skupina: 'mraz', nazev: 'Poslední jarní mráz', jedn: 'číslo dekády', agreg: 'zvl', presne: true, barva: '#7e57c2' },
    { id: 'podzimniMraz', skupina: 'mraz', nazev: 'První podzimní mráz', jedn: 'číslo dekády', agreg: 'zvl', presne: true, barva: '#5c6bc0' },
    { id: 'sdii', skupina: 'srazky', nazev: 'Intenzita srážek', jedn: 'mm/den', agreg: 'vazeny', presne: false, barva: '#66bb6a', riziko: 'nad' },
    { id: 'srazky', skupina: 'srazky', nazev: 'Úhrn srážek', jedn: 'mm', agreg: 'soucet', presne: false, barva: '#388e3c', riziko: 'pod' }
  ];
  const PARAM = {}; PARAMETRY.forEach(p => { PARAM[p.id] = p; });

  function maskaVse() { return new Uint8Array(36).fill(1); }
  function jeVse(m) { for (let i = 0; i < 36; i++) if (!m[i]) return false; return true; }

  // ---------- roční ukazatele z dekódované buňky ----------
  // maska: Uint8Array(36) kalendářních dekád (0 = 1.–10. 1.), výchozí všechny. Výsledek: Float32Array(NY) pro každé id.
  function rocniUkazatele(K, R, roc, jih, maska, y0, y1) {
    const NY = K.NY, NF = K.NF, IX = K.IX;
    y0 = y0 === undefined ? 0 : y0; y1 = y1 === undefined ? NY - 1 : y1;
    maska = maska || maskaVse();
    const vse = jeVse(maska), nm = maska.reduce((a, b) => a + b, 0);
    const O = {}; PARAMETRY.forEach(p => { O[p.id] = new Float32Array(NY).fill(NaN); });
    for (let y = y0; y <= y1; y++) {
      const o = (y - y0) * 36 * NF;
      const v = (d, f) => R[o + d * NF + IX[f]];
      let s = { tg: 0, tn: 0, tx: 0 }, tnn = Infinity, txx = -Infinity, soucet = {}, wsdii = 0, sr = 0, srAll = 0, vegDek = 0;
      COUNTS.forEach(c => { soucet[c] = 0; });
      for (let d = 0; d < 36; d++) {
        const tg = v(d, 'tg');
        if (tg >= 5) vegDek++;
        srAll += v(d, 'rr1') * v(d, 'sdii');
        if (!maska[d]) continue;
        s.tg += tg; s.tn += v(d, 'tn'); s.tx += v(d, 'tx');
        tnn = Math.min(tnn, v(d, 'tnn')); txx = Math.max(txx, v(d, 'txx'));
        COUNTS.forEach(c => { soucet[c] += v(d, c); });
        wsdii += v(d, 'rr1') * v(d, 'sdii'); sr += v(d, 'rr1');
      }
      O.tg[y] = vse ? roc.tg[y] : s.tg / Math.max(nm, 1);
      O.tn[y] = s.tn / Math.max(nm, 1); O.tx[y] = s.tx / Math.max(nm, 1);
      O.tnn[y] = tnn; O.txx[y] = vse ? roc.txx[y] : txx;
      COUNTS.forEach(c => { O[c][y] = vse ? roc.soucty[c][y] : soucet[c]; });
      O.sdii[y] = sr > 0 ? wsdii / sr : 0;
      O.srazky[y] = wsdii;
      // mráz: ls/fa z databáze (číslo dekády 1..18 v „posunutém“ roce, na jihu od července)
      const ls = roc.posledniJarniMraz[y], fa = roc.prvniPodzimniMraz[y];
      const lsi = ls - 1, fai = fa > 0 ? fa + 17 : 36;
      O.bezmrazi[y] = (fai - lsi - 1) * DEK_DNU;
      const kal = i => ((i + (jih ? 18 : 0)) % 36) + 1;            // posunutý index -> číslo kalendářní dekády 1..36
      O.jarniMraz[y] = ls > 0 ? kal(lsi) : NaN;
      O.podzimniMraz[y] = fa > 0 ? kal(fai) : NaN;
      O.vegSezona[y] = vegDek * DEK_DNU;
    }
    return O;
  }


  // ---------- rychlý výpočet pro mapu: jen potřebné veličiny, bez vynucování ročních součtů ----------
  // Přesné roční hodnoty se při výběru celého roku berou přímo z uložených ročních dat (bez dekódování modelu).
  const PRESNE_ROCNI = ['tg', 'txx', 'fd', 'tr', 'rr1', 'r20', 'su20', 'su25', 'su30', 'su35', 'jarniMraz', 'podzimniMraz', 'bezmrazi'];
  function jeJenRocni(id, maska) { return jeVse(maska) && PRESNE_ROCNI.indexOf(id) >= 0; }
  function featuryRychle(K, ids, potencialni, maska) {
    const s = new Set();
    if (potencialni) ['tg', 'tnn', 'rr1', 'sdii'].forEach(f => s.add(f));
    ids.forEach(id => {
      if (jeJenRocni(id, maska)) return;
      if (['tg', 'tn', 'tx', 'tnn', 'txx', 'fd', 'tr', 'rr1', 'r20', 'su20', 'su25', 'su30', 'su35'].indexOf(id) >= 0) s.add(id);
      else if (id === 'sdii' || id === 'srazky') { s.add('rr1'); s.add('sdii'); }
      else if (id === 'vegSezona') s.add('tg');
      else if (id === 'bezmrazi' || id === 'jarniMraz' || id === 'podzimniMraz') { }
    });
    return Array.from(s).map(f => K.IX[f]);
  }
  function potrebujeModel(K, ids, potencialni, maska) { return potencialni || ids.some(id => !jeJenRocni(id, maska)); }
  function potrebujeTnn(ids, potencialni) { return potencialni || ids.indexOf('tnn') >= 0; }

  // hodnota parametru pro každý rok y0..y1 (Float32Array(NY)); P z K.dekodujRychle (jen potřebné veličiny), roc z K.rocni(tile,k,y0,y1)
  function ukazatelRychle(K, P, roc, jih, maska, y0, y1, id) {
    const NY = K.NY, NF = K.NF, IX = K.IX, o = new Float32Array(NY).fill(NaN), vse = jeVse(maska), nm = maska.reduce((a, b) => a + b, 0);
    for (let y = y0; y <= y1; y++) {
      const b = (y - y0) * 36 * NF, v = (d, f) => P[b + d * NF + IX[f]];
      if (jeJenRocni(id, maska)) {
        if (id === 'tg') o[y] = roc.tg[y]; else if (id === 'txx') o[y] = roc.txx[y];
        else if (id === 'bezmrazi' || id === 'jarniMraz' || id === 'podzimniMraz') {
          const ls = roc.posledniJarniMraz[y], fa = roc.prvniPodzimniMraz[y], lsi = ls - 1, fai = fa > 0 ? fa + 17 : 36;
          const kal = i => ((i + (jih ? 18 : 0)) % 36) + 1;
          o[y] = id === 'bezmrazi' ? (fai - lsi - 1) * DEK_DNU : id === 'jarniMraz' ? (ls > 0 ? kal(lsi) : NaN) : (fa > 0 ? kal(fai) : NaN);
        } else o[y] = roc.soucty[id][y];
        continue;
      }
      let s = 0, mn = Infinity, mx = -Infinity, w = 0, sr = 0, dv = 0;
      for (let d = 0; d < 36; d++) {
        if (id === 'vegSezona') { if (v(d, 'tg') >= 5) dv++; continue; }
        if (id === 'bezmrazi' || id === 'jarniMraz' || id === 'podzimniMraz') continue;
        if (!maska[d]) continue;
        if (id === 'sdii' || id === 'srazky') { const r = Math.max(v(d, 'rr1'), 0); w += r; sr += r * v(d, 'sdii'); continue; }
        const x = id === 'tnn' || id === 'txx' || id === 'tg' || id === 'tn' || id === 'tx' ? v(d, id) : Math.max(v(d, id), 0);
        s += x; mn = Math.min(mn, x); mx = Math.max(mx, x);
      }
      if (id === 'vegSezona') o[y] = dv * DEK_DNU;
      else if (id === 'bezmrazi' || id === 'jarniMraz' || id === 'podzimniMraz') {
        const ls = roc.posledniJarniMraz[y], fa = roc.prvniPodzimniMraz[y], lsi = ls - 1, fai = fa > 0 ? fa + 17 : 36, kal = i => ((i + (jih ? 18 : 0)) % 36) + 1;
        o[y] = id === 'bezmrazi' ? (fai - lsi - 1) * DEK_DNU : id === 'jarniMraz' ? (ls > 0 ? kal(lsi) : NaN) : (fa > 0 ? kal(fai) : NaN);
      }
      else if (id === 'sdii') o[y] = w > 0 ? sr / w : 0;
      else if (id === 'srazky') o[y] = (vse && w > 0 && roc.soucty.rr1[y] > 0) ? sr * roc.soucty.rr1[y] / w : sr;
      else if (id === 'tnn') o[y] = mn; else if (id === 'txx') o[y] = mx;
      else if (id === 'tg' || id === 'tn' || id === 'tx') o[y] = s / Math.max(nm, 1);
      else o[y] = s;                                         // počty dní při výběru části roku
    }
    return o;
  }

  // ---------- statistiky nad řadou ----------
  function okno(serie, i0, i1) { const o = []; for (let i = Math.max(0, i0); i <= Math.min(serie.length - 1, i1); i++) if (serie[i] === serie[i]) o.push(serie[i]); return o; }
  function prumer(a) { return a.length ? a.reduce((x, y) => x + y, 0) / a.length : NaN; }
  function percentil(a, p) {                       // p v 0..100, lineární interpolace
    if (!a.length) return NaN;
    const s = Array.from(a).sort((x, y) => x - y), x = (p / 100) * (s.length - 1), i = Math.floor(x), f = x - i;
    return i + 1 < s.length ? s[i] + (s[i + 1] - s[i]) * f : s[i];
  }
  function riziko(a, mez, smer) {                  // podíl let, kdy hodnota je pod / nad mezí (0..1)
    if (!a.length) return NaN;
    return a.filter(v => smer === 'pod' ? v < mez : v > mez).length / a.length;
  }
  function trend(serie) {                          // sklon na 10 let (jednotky řady) a koeficient R²
    const n = [], y = [];
    serie.forEach((v, i) => { if (v === v) { n.push(i); y.push(v); } });
    const m = n.length; if (m < 3) return { sklon10: NaN, r2: NaN };
    const mx = prumer(n), my = prumer(y); let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < m; i++) { sxy += (n[i] - mx) * (y[i] - my); sxx += (n[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
    return { sklon10: sxx ? 10 * sxy / sxx : NaN, r2: sxx && syy ? (sxy * sxy) / (sxx * syy) : NaN };
  }
  function klouzavy(serie, n) {                    // centrovaný klouzavý průměr délky n
    const o = new Float32Array(serie.length).fill(NaN), h = Math.floor(n / 2);
    for (let i = 0; i < serie.length; i++) { const a = okno(serie, i - h, i + h); if (a.length >= Math.ceil(n / 2)) o[i] = prumer(a); }
    return o;
  }

  // ---------- podklad pro potenciál: roční klima ----------
  // USDA: dolní hranice zóny 1..13 (průměrné roční minimum, °C)
  const ZONA_SPODNI = [null, -51.1, -45.6, -40, -34.4, -28.9, -23.3, -17.8, -12.2, -6.7, -1.1, 4.4, 10, 15.6];
  function zonaZTeploty(t) { let z = 1; for (let i = 1; i <= 13; i++) if (t >= ZONA_SPODNI[i]) z = i; return z; }

  function klimaRoku(K, R, roc, jih, y0, y1) {
    const NY = K.NY, NF = K.NF, IX = K.IX;
    y0 = y0 === undefined ? 0 : y0; y1 = y1 === undefined ? NY - 1 : y1;
    const o = { NY, tnn: new Float32Array(NY), srazky: new Float32Array(NY), bezmrazi: new Float32Array(NY),
      tgSezona: new Float32Array(NY), gdd5: new Float32Array(NY), gdd10: new Float32Array(NY), gdd5Sezona: new Float32Array(NY), gdd10Sezona: new Float32Array(NY) };
    for (let y = y0; y <= y1; y++) {
      const b = (y - y0) * 36 * NF;
      const ls = roc.posledniJarniMraz[y], fa = roc.prvniPodzimniMraz[y], lsi = ls - 1, fai = fa > 0 ? fa + 17 : 36;
      let tnn = Infinity, sr = 0, ts = 0, ns = 0, g5 = 0, g10 = 0, g5s = 0, g10s = 0;
      for (let d = 0; d < 36; d++) {
        const tg = R[b + d * NF + IX.tg];
        tnn = Math.min(tnn, R[b + d * NF + IX.tnn]);
        sr += R[b + d * NF + IX.rr1] * R[b + d * NF + IX.sdii];
        const a5 = Math.max(tg - 5, 0) * DEK_DNU, a10 = Math.max(tg - 10, 0) * DEK_DNU;
        g5 += a5; g10 += a10;
        const s = (d - (jih ? 18 : 0) + 36) % 36;                       // posunutý index dekády
        if (s > lsi && s < fai) { ts += tg; ns++; g5s += a5; g10s += a10; }
      }
      o.tnn[y] = tnn; o.srazky[y] = sr; o.bezmrazi[y] = (fai - lsi - 1) * DEK_DNU;
      o.tgSezona[y] = ns ? ts / ns : NaN; o.gdd5[y] = g5; o.gdd10[y] = g10; o.gdd5Sezona[y] = g5s; o.gdd10Sezona[y] = g10s;
    }
    return o;
  }

  // ---------- druhy a potenciál ----------
  function nactiDruhy(json) {                      // {pole, radky} -> pole objektů
    return json.radky.map(r => { const d = {}; json.pole.forEach((k, i) => { d[k] = r[i]; }); return d; });
  }
  const sat = x => Math.max(0, Math.min(1, x));
  const rampaNahoru = (x, nula, jedna) => nula === jedna ? (x >= jedna ? 1 : 0) : sat((x - nula) / (jedna - nula));
  function lichobeznik(x, amin, omin, omax, amax) {  // 0 mimo absolutní meze, 1 v optimu, lineárně mezi
    if (x <= amin || x >= amax) return 0;
    if (omin !== null && omax !== null) {
      if (x >= omin && x <= omax) return 1;
      return x < omin ? sat((x - amin) / Math.max(omin - amin, 1e-9)) : sat((amax - x) / Math.max(amax - omax, 1e-9));
    }
    return 1;
  }
  const ne = v => v === null || v === undefined || v !== v;

  // způsobilost druhu v jednom roce/období (0..1): nejslabší z použitelných kritérií (zákon minima)
  function vhodnost(dr, k) {
    const f = [];
    const jednol = dr.jednoleta === 1;
    const baze5 = (dr.oblast || '').indexOf('borealni') >= 0;           // BEDD_C: u chladných oblastí od 5 °C, jinak od 10 °C
    const zona = zonaZTeploty(k.tnn);
    if (!jednol) {
      if (!ne(dr.usda_a)) f.push(rampaNahoru(k.tnn, ZONA_SPODNI[dr.usda_a] - 2, ZONA_SPODNI[dr.usda_a] + 2));
      if (!ne(dr.letal_klid) && dr.ecocrop_shoda !== 'ne') f.push(rampaNahoru(k.tnn, dr.letal_klid - 2, dr.letal_klid + 2));
    } else {
      const potreba = !ne(dr.cyklus_min) ? dr.cyklus_min : (!ne(dr.veg_min) ? dr.veg_min : null);
      if (potreba !== null) f.push(rampaNahoru(k.bezmrazi, potreba * 0.6, potreba));
    }
    if (!ne(dr.usda_b)) f.push(1 - sat((zona - dr.usda_b) / 2));      // příliš teplo pro druh (horní zóna)
    if (!ne(dr.t_abs_min) && !ne(dr.t_abs_max) && !ne(k.tgSezona)) {
      f.push(lichobeznik(k.tgSezona, dr.t_abs_min, dr.t_opt_min ?? null, dr.t_opt_max ?? null, dr.t_abs_max));
    }
    if (!ne(dr.bedd) && dr.bedd > 0) {
      const g = baze5 ? (jednol ? k.gdd5Sezona : k.gdd5) : (jednol ? k.gdd10Sezona : k.gdd10);
      f.push(sat(g / dr.bedd));
    }
    const rmin = !ne(dr.r_abs_min) ? dr.r_abs_min : (!ne(dr.min_srazky) ? dr.min_srazky : null);
    if (rmin !== null) f.push(rampaNahoru(k.srazky, rmin * 0.75, rmin));
    if (!ne(dr.r_abs_max)) f.push(1 - rampaNahoru(k.srazky, dr.r_abs_max, dr.r_abs_max * 1.25));
    return f.length ? Math.min.apply(null, f) : NaN;
  }

  // potenciál v okně let [rokOd, rokDo] (indexy 0..NY-1): klima = průměr přes okno; vaha = kcal na 100 g
  function potencial(druhy, Ku, i0, i1, opt) {
    opt = opt || {};
    const prah = opt.prah ?? 0.5;
    const ag = a => prumer(okno(a, i0, i1));
    const k = { tnn: ag(Ku.tnn), srazky: ag(Ku.srazky), bezmrazi: ag(Ku.bezmrazi), tgSezona: ag(Ku.tgSezona),
      gdd5: ag(Ku.gdd5), gdd10: ag(Ku.gdd10), gdd5Sezona: ag(Ku.gdd5Sezona), gdd10Sezona: ag(Ku.gdd10Sezona) };
    let sw = 0, ss = 0, st = 0, n = 0; const proDruh = opt.detail ? [] : null;
    for (const dr of druhy) {
      const w = dr.kcal > 0 ? dr.kcal : 0; if (!w) continue;
      const s = vhodnost(dr, k); if (s !== s) continue;
      sw += w; ss += w * s; st += w * (s >= prah ? 1 : 0); n++;
      if (proDruh) proDruh.push({ klic: dr.klic, cesky: dr.cesky, s });
    }
    return { podil: sw ? ss / sw : NaN, tvrdy: sw ? st / sw : NaN, vahaCelkem: sw, n, klima: k, proDruh };
  }

  const API = { PARAMETRY, PARAM, SKUPINY, maskaVse, rocniUkazatele, okno, prumer, percentil, riziko, trend, klouzavy,
    klimaRoku, nactiDruhy, vhodnost, potencial, zonaZTeploty, DEK_DNU, featuryRychle, potrebujeModel, potrebujeTnn, ukazatelRychle, jeJenRocni };
  if (typeof module !== 'undefined' && module.exports) module.exports = API; else root.KlimaAgro = API;
})(typeof window !== 'undefined' ? window : globalThis);
