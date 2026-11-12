/*
 * klima_dekoder.js – dekodér dlaždic databáze klimatických indikátorů (bez knihoven, prohlížeč i Node).
 *
 * Použití v prohlížeči (stejně jako u staré databáze):
 *   const reg  = await (await fetch('registrace_klima.json')).json();
 *   const kody = await (await fetch(reg.soubor_kody)).arrayBuffer();      // jednou, 1,9 MB
 *   const K    = new KlimaDekoder(reg, kody);
 *   // dlaždice: POST {"tile_ids":[...]} na get_klima_post.php -> binární stream
 *   const bloky = K.parsujStream(await odpoved.arrayBuffer());            // pole Uint8Array|null v pořadí dotazu
 *   const tile  = await K.rozbalGzip(bloky[0]);                           // jedna dlaždice
 *   const R     = K.dekodujBunku(tile, indexBunkyVDlazdici, jeJih);       // Float32Array(NY*36*NF)
 *       // hodnota: R[(rok*36 + dekada)*NF + K.IX.tnn], rok = 0 pro K.reg.rok_od
 *   const roc   = K.rocni(tile, indexBunkyVDlazdici);                     // přesné roční hodnoty bez dekódování
 *
 * Pořadí buněk v dlaždici: řádek po řádku podle masky pevniny (stejně jako u staré databáze).
 * Výsledek dekódování je stejný jako výsledek Pythonu (dekodovat v stavba_db.py), až na zanedbatelné rozdíly
 * ve volbě mezi dekádami se stejným pořadím při dělení ročních součtů o jednotky.
 */
(function (root) {
  'use strict';

  const fr = Math.fround;

  function f16(h) {
    const s = (h & 0x8000) ? -1 : 1, e = (h >> 10) & 0x1f, m = h & 0x3ff;
    if (e === 0) return s * m * 5.960464477539063e-8;
    if (e === 31) return m ? NaN : s * Infinity;
    return s * (1 + m / 1024) * Math.pow(2, e - 15);
  }

  function rint(x) {                       // zaokrouhlení na nejbližší, 0,5 na sudé (jako numpy)
    const f = Math.floor(x), d = x - f;
    if (d < 0.5) return f;
    if (d > 0.5) return f + 1;
    return (f % 2 === 0) ? f : f + 1;
  }

  function interp(x, xp, fp) {             // np.interp: mimo rozsah krajní hodnota
    const n = xp.length;
    if (x <= xp[0]) return fp[0];
    if (x >= xp[n - 1]) return fp[n - 1];
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xp[mid] <= x) lo = mid; else hi = mid; }
    return fp[lo] + (fp[lo + 1] - fp[lo]) / (xp[lo + 1] - xp[lo]) * (x - xp[lo]);
  }

  // Rozdělí roční součet S do 36 dekád podle vah Q (viz skaluj v odhad_vq2.py). Všechna pole délky 36.
  function skaluj(Q, S, W, baz, cap, Bf, out) {
    const w = new Float64Array(36), x = new Float64Array(36), hlava = new Float64Array(36), vol = new Float64Array(36);
    const frac = new Float64Array(36), fl = new Float64Array(36);
    let sb = 0; for (let i = 0; i < 36; i++) sb += baz[i];
    const zb = Math.max(S - sb, 0);
    let sw = 0; for (let i = 0; i < 36; i++) { w[i] = Q[i] * W[i]; sw += w[i]; }
    if (sw <= 0) { sw = 0; for (let i = 0; i < 36; i++) { w[i] = Bf[i] * W[i]; sw += w[i]; } }
    if (sw <= 0) { sw = 0; for (let i = 0; i < 36; i++) { w[i] = W[i]; sw += w[i]; } }
    for (let i = 0; i < 36; i++) x[i] = sw > 0 ? w[i] / Math.max(sw, 1e-12) * zb : 0;
    for (let i = 0; i < 36; i++) hlava[i] = Math.max(cap[i] - baz[i], 0);
    for (let it = 0; it < 3; it++) {
      let sx = 0, sv = 0;
      for (let i = 0; i < 36; i++) { x[i] = Math.min(x[i], hlava[i]); vol[i] = Math.max(hlava[i] - x[i], 0) * W[i]; sx += x[i]; sv += vol[i]; }
      const manko = Math.min(zb - sx, sv);
      for (let i = 0; i < 36; i++) x[i] += sv > 0 ? vol[i] / Math.max(sv, 1e-12) * manko : 0;
    }
    let sfl = 0;
    for (let i = 0; i < 36; i++) { x[i] = Math.min(x[i], hlava[i]); fl[i] = Math.floor(x[i]); sfl += fl[i]; }
    const rem = Math.max(rint(zb - sfl), 0);
    for (let i = 0; i < 36; i++) frac[i] = (fl[i] < hlava[i] && W[i]) ? x[i] - fl[i] : -1;
    const poradi = Array.from({ length: 36 }, (_, i) => i).sort((a, b) => (frac[b] - frac[a]) || (a - b));
    const add = new Float64Array(36);
    for (let r = 0; r < 36 && r < rem; r++) { const i = poradi[r]; if (frac[i] >= 0) add[i] = 1; }
    for (let i = 0; i < 36; i++) out[i] = baz[i] + fl[i] + add[i];
  }

  class KlimaDekoder {
    constructor(reg, kodyBuf) {
      this.reg = reg;
      const M = this.M = reg.dekoder;
      this.NF = M.NF; this.NY = M.NY; this.NT = M.NT; this.nper = M.nper;
      this.IX = {}; M.FEAT.forEach((f, i) => { this.IX[f] = i; });
      // pokřivení (warp): uzly x (hranice) a y (kumulativní sklony)
      this.wx = []; this.wy = []; this.wp = [];
      M.FEAT.forEach((f) => {
        const h = M.warp[f].hranice, s = M.warp[f].sklony, y = [h[0]];
        for (let k = 0; k < s.length; k++) y.push(y[k] + s[k] * (h[k + 1] - h[k]));
        this.wx.push(h); this.wy.push(y); this.wp.push(M.priorita[f]);
      });
      // kódové knihy: C1 (K1,NF), C2 (K1,K2,NF), CA (K,42) jako float16 little-endian
      const dv = new DataView(kodyBuf);
      const n1 = M.K1 * M.NF, n2 = M.K1 * M.K2 * M.NF, n3 = M.K_ANOM * 42;
      const nacti = (start, n) => { const a = new Float32Array(n); for (let i = 0; i < n; i++) a[i] = f16(dv.getUint16((start + i) * 2, true)); return a; };
      this.C1 = nacti(0, n1); this.C2 = nacti(n1, n2); this.CA = nacti(n1 + n2, n3);
      this.yp = new Uint8Array(M.NY);                 // rok -> období (základ pro 3 období)
      { let y = 0; for (let i = 0; i < M.nper; i++) for (let k = 0; k < M.ny[i]; k++) this.yp[y++] = i; }
      this.qc = (M.strop_tnn === null || M.strop_tnn === undefined) ? null : Math.floor(M.strop_tnn / M.krok_tnn + 1e-9);
    }

    // ---------- přenos ----------
    parsujStream(buf) {                                // uint32 BE délka + gzip blok, v pořadí dotazu
      const u8 = new Uint8Array(buf), dv = new DataView(buf), out = []; let pos = 0;
      while (pos + 4 <= u8.length) {
        const len = dv.getUint32(pos, false); pos += 4;
        out.push(len > 0 ? u8.subarray(pos, pos + len) : null); pos += len;
      }
      return out;
    }

    // Dlaždice se stahuje po částech: 'rocni' (přesné roční hodnoty), 'model' (základ, index, výjimky SDII), 'tnn' (TNn po dekádách).
    // Každá část se rozbalí do téhož objektu dlaždice: tile = await K.rozbalGzip(gz, 'rocni', tile)
    async rozbalGzip(gz, cast, tile) {                  // prohlížeč: DecompressionStream
      const ds = new DecompressionStream('gzip');
      const stream = new Blob([gz]).stream().pipeThrough(ds);
      return this.rozbalCast(new Uint8Array(await new Response(stream).arrayBuffer()), cast, tile);
    }

    rozbalCast(u8, cast, tile) {                        // už rozbalený (gunzip) obsah jedné části
      tile = tile || { n: 0, sek: {}, casti: {}, exc: { n: 0, bunka: new Uint16Array(0), t: new Uint16Array(0), val: new Uint16Array(0) } };
      const dv = new DataView(u8.buffer, u8.byteOffset, u8.byteLength);
      const n = dv.getUint32(0, true), nrec = dv.getUint32(4, true);
      let pos = 8;
      for (const s of this.M.casti[cast]) {
        const len = n * s.m, a = u8.slice(pos, pos + len); pos += len;
        if (s.lag > 0) {
          for (let r = 0; r < n; r++) {
            const b = r * s.m;
            for (let i = s.lag; i < s.m; i++) a[b + i] = (a[b + i] + a[b + i - s.lag]) & 255;
          }
        }
        tile.sek[s.jmeno] = a;
      }
      if (cast === 'model') {
        const e = u8.subarray(pos, pos + 6 * nrec);
        const exc = { n: nrec, bunka: new Uint16Array(nrec), t: new Uint16Array(nrec), val: new Uint16Array(nrec) };
        for (let i = 0; i < nrec; i++) {
          exc.bunka[i] = e[i] | (e[nrec + i] << 8); exc.t[i] = e[2 * nrec + i] | (e[3 * nrec + i] << 8); exc.val[i] = e[4 * nrec + i] | (e[5 * nrec + i] << 8);
        }
        tile.exc = exc;
      }
      tile.n = n; tile.casti[cast] = true;
      return tile;
    }

    // ---------- přístup k uloženým hodnotám ----------
    _idx(tile, k, j) {                                  // index měsíce (0..K_ANOM-1) pro rok*12+měsíc = j
      const NY12 = this.NY * 12; let v = tile.sek.idx_lo[k * NY12 + j];
      const mb = Math.ceil(NY12 / 8);
      for (let p = 0; p < this.M.nb_idx - 8; p++) v |= ((tile.sek['idx_hi' + p][k * mb + (j >> 3)] >> (7 - (j & 7))) & 1) << (8 + p);
      return v;
    }
    _suma(tile, k, nm, y) {
      const NY = this.NY, mb = Math.ceil(NY / 8);
      return tile.sek['sum_' + nm + '_lo'][k * NY + y] | (((tile.sek['sum_' + nm + '_hi'][k * mb + (y >> 3)] >> (7 - (y & 7))) & 1) << 8);
    }
    _r(tile, k, t) {                                    // zbytek TNn v krocích (zigzag); bez části 'tnn' 0
      if (!tile.sek.r_lo) return 0;
      const zz = tile.sek.r_lo[k * this.NT + t] | (tile.sek.r_hi[k * this.NT + t] << 8);
      return (zz >>> 1) ^ -(zz & 1);
    }

    rocni(tile, k, y0, y1) {                            // přesné roční hodnoty (bez dekódování dekád); y0..y1 = jen tyto roky
      const NY = this.NY, o = { soucty: {}, tg: new Float32Array(NY), txx: new Float32Array(NY),
        posledniJarniMraz: new Uint8Array(NY), prvniPodzimniMraz: new Uint8Array(NY) };
      y0 = y0 === undefined ? 0 : y0; y1 = y1 === undefined ? NY - 1 : y1;
      for (const f of this.M.COUNTS) { const a = new Uint16Array(NY); for (let y = y0; y <= y1; y++) a[y] = this._suma(tile, k, f, y); o.soucty[f] = a; }
      for (let y = y0; y <= y1; y++) {
        o.tg[y] = tile.sek.tgq[k * NY + y] * 0.5 - 75; o.txx[y] = tile.sek.txxq[k * NY + y] * 0.5 - 75;
        o.posledniJarniMraz[y] = tile.sek.ls[k * NY + y]; o.prvniPodzimniMraz[y] = tile.sek.fa[k * NY + y];
      }
      return o;                                         // ls/fa: číslo dekády 1..18 (0 = žádný mráz), jižní polokoule od července
    }

    _unwarp(y, j) { return fr(interp(y / this.wp[j], this.wy[j], this.wx[j])); }

    // dekódované spojité hodnoty (před vynucením ročních hodnot): P a základ Bf, oba Float32Array(NY*36*NF)
    _predikce(tile, k, y0, y1, fs) {                    // fs (volitelné): jen tyto veličiny (indexy), bez základu Bf
      const { NF, NY, nper, M } = this, K2 = M.K2;
      const fl = fs || Array.from({ length: NF }, (_, i) => i);
      y0 = y0 === undefined ? 0 : y0; y1 = y1 === undefined ? NY - 1 : y1;
      const tl = (y1 - y0 + 1) * 36 * NF;
      const P = new Float32Array(tl), Bf = new Float32Array(tl);
      const sq = new Float32Array(nper * 12);
      for (let i = 0; i < nper * 12; i++) sq[i] = fr(M.SIG_MIN * Math.pow(2, tile.sek.sc[k * nper * 12 + i] / 8));
      const Brec = new Float32Array(nper * 36 * NF);
      for (let i = 0; i < nper * 36; i++) {
        const a = tile.sek.p1[k * nper * 36 + i], b = tile.sek.p2[k * nper * 36 + i];
        for (const f of fl) Brec[i * NF + f] = fr(this.C1[a * NF + f] + this.C2[(a * K2 + b) * NF + f]);
      }
      for (let y = y0; y <= y1; y++) {
        const pi = this.yp[y];
        for (let mo = 0; mo < 12; mo++) {
          const ci = this._idx(tile, k, y * 12 + mo) * 42, s = sq[pi * 12 + mo];
          for (let dd = 0; dd < 3; dd++) {
            const d = mo * 3 + dd, o = ((y - y0) * 36 + d) * NF;
            for (const f of fl) {
              const br = Brec[(pi * 36 + d) * NF + f];
              const pre = fr(br + fr(this.CA[ci + dd * NF + f] * s));
              P[o + f] = this._unwarp(pre, f); if (!fs) Bf[o + f] = this._unwarp(br, f);
            }
          }
        }
      }
      return { P, Bf };
    }

    // rychlé dekódování pro mapu: jen vybrané veličiny (např. [IX.tg, IX.tnn]) bez vynucování ročních součtů.
    // Pokud je načtená část 'tnn', jsou dekádní TNn přesné (do stropu). Výsledek: Float32Array((y1-y0+1)*36*NF).
    dekodujRychle(tile, k, y0, y1, fs) {
      const { NF, IX, M } = this, kr = M.krok_tnn, qc = this.qc;
      const { P } = this._predikce(tile, k, y0, y1, fs);
      if ((!fs || fs.indexOf(IX.tnn) >= 0) && tile.sek.r_lo) {
        for (let y = y0; y <= y1; y++) for (let d = 0; d < 36; d++) {
          const o = ((y - y0) * 36 + d) * NF + IX.tnn, ppt = P[o];
          let qp = rint(ppt / kr); if (qc !== null) qp = Math.min(qp, qc);
          const qt = qp + this._r(tile, k, y * 36 + d);
          P[o] = fr((qc !== null && qt >= qc) ? Math.max(ppt, qc * kr) : qt * kr);
        }
      }
      return P;
    }

    // rychlý odhad jedné dekády (bez vynucení ročních hodnot) – např. pro mapy; Float32Array(NF)
    dekodujDekadu(tile, k, rok, dekada) {
      const { P } = this._predikce(tile, k, rok, rok), o = dekada * this.NF;
      return P.slice(o, o + this.NF);
    }

    // plné dekódování jedné buňky včetně přesných ročních hodnot a výjimek SDII: Float32Array(NY*36*NF)
    // y0, y1 (volitelné): dekóduje jen roky y0..y1; výsledek pak má Float32Array((y1-y0+1)*36*NF) s rokem y0 na začátku
    dekodujBunku(tile, k, jih, y0, y1) {
      const { NF, NY, NT, IX, M } = this, COUNTS = M.COUNTS;
      y0 = y0 === undefined ? 0 : y0; y1 = y1 === undefined ? NY - 1 : y1;
      const { P, Bf } = this._predikce(tile, k, y0, y1);
      const R = new Float32Array(P);
      const c11 = new Float64Array(36).fill(11), ones = new Uint8Array(36).fill(1), nuly = new Float64Array(36);
      const kr = M.krok_tnn, qc = this.qc;
      const ser = (A, y, j, mx) => { const a = new Float64Array(36); for (let d = 0; d < 36; d++) { const v = A[((y - y0) * 36 + d) * NF + j]; a[d] = mx ? Math.max(v, 0) : v; } return a; };
      const tmp = new Float64Array(36), rr1 = new Float64Array(36), cap = new Float64Array(36);
      const put = (y, j, a) => { for (let d = 0; d < 36; d++) R[((y - y0) * 36 + d) * NF + j] = a[d]; };
      for (let y = y0; y <= y1; y++) {
        const S = {}; for (const f of COUNTS) S[f] = this._suma(tile, k, f, y);
        const ls = tile.sek.ls[k * NY + y], fa = tile.sek.fa[k * NY + y];
        const lsi = ls - 1, fai = fa > 0 ? fa + 17 : 36;
        const W = new Uint8Array(36), baz = new Float64Array(36);
        for (let d = 0; d < 36; d++) {
          const s = (d - (jih ? 18 : 0) + 36) % 36;
          const nula = (s <= 17 && s > lsi) || (s >= 18 && s < fai);
          const hr = (s === lsi && lsi >= 0) || (s === fai && fai < 36);
          W[d] = nula ? 0 : 1; baz[d] = hr ? 1 : 0;
        }
        skaluj(ser(P, y, IX.fd, true), S.fd, W, baz, c11, ser(Bf, y, IX.fd, true), tmp); put(y, IX.fd, tmp);
        skaluj(ser(P, y, IX.rr1, true), S.rr1, ones, nuly, c11, ser(Bf, y, IX.rr1, true), rr1); put(y, IX.rr1, rr1);
        skaluj(ser(P, y, IX.r20, true), S.r20, ones, nuly, rr1, ser(Bf, y, IX.r20, true), tmp); put(y, IX.r20, tmp);
        skaluj(ser(P, y, IX.tr, true), S.tr, ones, nuly, c11, ser(Bf, y, IX.tr, true), tmp); put(y, IX.tr, tmp);
        cap.set(c11);
        for (const f of ['su20', 'su25', 'su30', 'su35']) {
          skaluj(ser(P, y, IX[f], true), S[f], ones, nuly, cap, ser(Bf, y, IX[f], true), tmp); put(y, IX[f], tmp); cap.set(tmp);
        }
        // teploty
        const tgm = tile.sek.tgq[k * NY + y] * 0.5 - 75, txxm = tile.sek.txxq[k * NY + y] * 0.5 - 75;
        let st = 0; for (let d = 0; d < 36; d++) st += P[((y - y0) * 36 + d) * NF + IX.tg];
        const delta = tgm - st / 36;
        const tnn = new Float64Array(36), txx = new Float64Array(36), tn = new Float64Array(36), tx = new Float64Array(36), tg = new Float64Array(36);
        let amax = 0;
        for (let d = 0; d < 36; d++) {
          const o = ((y - y0) * 36 + d) * NF;
          const ppt = P[o + IX.tnn];
          let qp = rint(ppt / kr); if (qc !== null) qp = Math.min(qp, qc);
          const qt = qp + this._r(tile, k, y * 36 + d);
          tnn[d] = fr((qc !== null && qt >= qc) ? Math.max(ppt, qc * kr) : qt * kr);
          txx[d] = fr(Math.min(fr(P[o + IX.txx] + delta), txxm));
          tn[d] = fr(P[o + IX.tn] + delta); tx[d] = fr(P[o + IX.tx] + delta); tg[d] = fr(P[o + IX.tg] + delta);
          if (txx[d] > txx[amax]) amax = d;
        }
        txx[amax] = txxm;
        for (let d = 0; d < 36; d++) {
          const o = ((y - y0) * 36 + d) * NF;
          const tn2 = Math.max(tn[d], tnn[d]), tx2 = Math.min(tx[d], txx[d]);
          R[o + IX.tnn] = tnn[d]; R[o + IX.txx] = txx[d]; R[o + IX.tn] = tn2; R[o + IX.tx] = tx2;
          R[o + IX.tg] = Math.min(Math.max(tg[d], tn2), tx2);
          R[o + IX.sdii] = R[o + IX.rr1] <= 0 ? 0 : Math.max(P[o + IX.sdii], 1.0);
        }
      }
      const ex = tile.exc;                               // výjimky SDII (skutečné extrémy bez ořezu)
      for (let i = 0; i < ex.n; i++) if (ex.bunka[i] === k && ex.t[i] >= y0 * 36 && ex.t[i] < (y1 + 1) * 36) R[(ex.t[i] - y0 * 36) * NF + IX.sdii] = f16(ex.val[i]);
      return R;
    }
  }

  const API = { KlimaDekoder, f16 };
  if (typeof module !== 'undefined' && module.exports) module.exports = API; else root.KlimaDekoder = KlimaDekoder;
})(typeof window !== 'undefined' ? window : globalThis);
