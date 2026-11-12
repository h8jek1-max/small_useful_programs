#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ODHAD VQ3 – varianta bezExt + TNn po dekádách, chyby v jednotkách, percentily a bias

Jedním během (data se čtou jen jednou) porovná tři scénáře nad stejným indexem (varianta bezExt):

  S0  bezExt jako dosud: roční TNn, TXx, TG a součty počtů přesně, dekádní TNn z indexu
  S1  S0 + TNn PO DEKÁDÁCH uložené přesně (v zadaném kroku); roční TNn se z něj dopočítá, proud se vypustí
  S2  S1 + vypuštěné vybrané roční ukazatele (--vypust), jejich roční hodnota se pak jen odhadne z indexu

Výpis: pro každý scénář a každý indikátor v jeho jednotkách (dny, °C, mm/den):
RMSE, MAE, bias (se znaménkem), P95, P99 a maximum absolutní chyby (test na buňkách mimo trénink),
dále rozdělení podle zeměpisné šířky (tropy / mírné pásmo / chladné oblasti), bias extrémů podle
třídy skutečné hodnoty (nejchladnější TNn, nejteplejší TXx), roční chyby a velikosti po gzipu.

Vyžaduje odhad_vq2.py a odhad_indikatoru.py ve stejné složce. Nic nezapisuje kromě odhad_vq3.txt.

    python odhad_vq3.py                                   výchozí: krok 0.1 °C, strop +10 °C, vypustit su20,su25,txx
    python odhad_vq3.py --krok 0.5 --strop 10             TNn po 0,5 °C
    python odhad_vq3.py --strop zadny                     TNn v celém rozsahu
    python odhad_vq3.py --vypust su20,su25,su30,txx,tg    jiný výběr vypuštěných ročních ukazatelů
    python odhad_vq3.py --roky 10 --bunek 600 --trenink 2000    rychlá zkouška
Povolené hodnoty --vypust: su20 su25 su30 su35 tr r20 rr1 tg txx (FD a jeho termíny mrazů zůstávají vždy)
"""
import sys, time, math
import numpy as np
from odhad_vq2 import *          # sdílené funkce a nastavení (FEAT, IX, warp, kmeans, vynut, skaluj, ...)
import odhad_vq2 as v2
import odhad_indikatoru as oi

# =============================== NASTAVENÍ ===============================
VAR = arg("--varianta", "bezExtremu")
if VAR not in CONFIGS: sys.exit("Neznámá varianta: " + VAR)
KROK = float(arg("--krok", 0.1))
_s = arg("--strop", "10")
STROP = None if _s.lower() in ("zadny", "none", "-") else float(_s)
VYP = [x for x in arg("--vypust", "su20,su25,txx").split(",") if x]
POVOLENO = {"su20", "su25", "su30", "su35", "tr", "r20", "rr1", "tg", "txx"}
for x in VYP:
    if x not in POVOLENO: sys.exit(f"--vypust: nepovolená hodnota '{x}'")
ROCNI_JMENO = {"su20": "SU20", "su25": "SU25", "su30": "SU30", "su35": "SU35", "tr": "TR", "r20": "R20mm",
               "rr1": "RR1", "tg": "TG", "txx": "TXx"}
LOG = []
def p(s=""):
    print(s); LOG.append(s)

SC = [("S0", False, ()), ("S1", True, ()), ("S2", True, tuple(VYP))]
BW = 0.05; NBIN = 2001
PASMA = ["tropy (|šířka| < 23,5°)", "mírné pásmo (23,5–50°)", "chladné oblasti (> 50°)"]
KLS = {"tnn": ([-30, -15, 0, 10], ["< -30 °C", "-30 až -15", "-15 až 0", "0 až 10", "> 10 °C"]),
       "txx": ([30, 40], ["< 30 °C", "30 až 40", "> 40 °C"])}
ANNITEMS = list(ANN_EX) + [("TN roční průměr", "tn", "mean"), ("TX roční průměr", "tx", "mean"),
                           ("dekád s mrazem", "fd", "frost")]
JEDN = {"TG": "°C", "TNn": "°C", "TXx": "°C", "TN roční průměr": "°C", "TX roční průměr": "°C"}


def pasmo(lat):
    a = np.abs(lat)
    return np.where(a < 23.5, 0, np.where(a < 50, 1, 2)).astype(np.int32)


# =============================== VYNUCENÍ (rozšířené) ===============================
def vynut2(P, Bf, Pt, jih, tnn_dek=None, vyp=()):
    """Jako vynut z odhad_vq2, navíc: tnn_dek (n,NY,36) = dekódované dekádní TNn (nahradí roční vynucení),
    vyp = roční ukazatele, které se nevynucují (použije se součet/průměr/maximum z dekódovaných hodnot)."""
    n, ny_, _, _ = P.shape
    R = P.copy(); Bf = np.maximum(Bf, 0)
    ones = np.ones((n, ny_, 36), bool); nuly = np.zeros((n, ny_, 36)); c11 = np.full((n, ny_, 36), 11.0)
    g = lambda f: np.maximum(P[..., IX[f]], 0)
    S_ = lambda f: np.rint(g(f).sum(2)) if f in vyp else Pt[..., IX[f]].sum(2)
    ls, fa, ordi = mraz_kody(Pt[..., IX["fd"]], jih)
    inv = np.argsort(ordi, axis=1)
    s = np.arange(36)[None, None, :]
    lsi = ls.astype(int)[..., None] - 1
    fai = np.where(fa > 0, fa.astype(int) + 17, 36)[..., None]
    nula = ((s <= 17) & (s > lsi)) | ((s >= 18) & (s < fai))
    hr = ((s == lsi) & (lsi >= 0)) | ((s == fai) & (fai < 36))
    pick = lambda m: np.take_along_axis(m, inv[:, None, :], axis=2)
    R[..., IX["fd"]] = skaluj(g("fd"), S_("fd"), ~pick(nula), pick(hr).astype(np.float64), c11, Bf[..., IX["fd"]])
    R[..., IX["rr1"]] = skaluj(g("rr1"), S_("rr1"), ones, nuly, c11, Bf[..., IX["rr1"]])
    R[..., IX["r20"]] = skaluj(g("r20"), S_("r20"), ones, nuly, R[..., IX["rr1"]], Bf[..., IX["r20"]])
    R[..., IX["tr"]] = skaluj(g("tr"), S_("tr"), ones, nuly, c11, Bf[..., IX["tr"]])
    cap = c11
    for f in ("su20", "su25", "su30", "su35"):
        R[..., IX[f]] = skaluj(g(f), S_(f), ones, nuly, cap, Bf[..., IX[f]])
        cap = R[..., IX[f]]
    # --- teploty
    if "tg" in vyp: tg_m = P[..., IX["tg"]].mean(2)
    else: tg_m = kv05(Pt[..., IX["tg"]].mean(2))
    delta = (tg_m - P[..., IX["tg"]].mean(2))[..., None]
    for f in ("tg", "tn", "tx", "tnn", "txx"): R[..., IX[f]] = P[..., IX[f]] + delta
    if tnn_dek is not None:
        tnn = np.asarray(tnn_dek, np.float32).copy()
    else:
        tnn_m = kv05(Pt[..., IX["tnn"]].min(2))
        tnn = np.maximum(R[..., IX["tnn"]], tnn_m[..., None])
        np.put_along_axis(tnn, tnn.argmin(2)[..., None], tnn_m[..., None], axis=2)
    if "txx" in vyp:
        txx = R[..., IX["txx"]]
    else:
        txx_m = kv05(Pt[..., IX["txx"]].max(2))
        txx = np.minimum(R[..., IX["txx"]], txx_m[..., None])
        np.put_along_axis(txx, txx.argmax(2)[..., None], txx_m[..., None], axis=2)
    R[..., IX["tnn"]] = tnn; R[..., IX["txx"]] = txx
    tn = np.maximum(R[..., IX["tn"]], tnn); tx = np.minimum(R[..., IX["tx"]], txx)
    R[..., IX["tn"]] = tn; R[..., IX["tx"]] = tx
    R[..., IX["tg"]] = np.minimum(np.maximum(R[..., IX["tg"]], tn), tx)
    R[..., IX["sdii"]] = np.where(R[..., IX["rr1"]] <= 0, 0, np.maximum(P[..., IX["sdii"]], 1.0))
    return R


# =============================== TNn PO DEKÁDÁCH ===============================
def tnn_kvant(Pt_t, Pp_t, krok, strop):
    """Vrací (qt, qp, qc): kvantizovaná pravda, kvantizovaná predikce (obojí oříznuté stropem), strop v krocích."""
    qt = np.rint(Pt_t / krok); qp = np.rint(Pp_t / krok); qc = None
    if strop is not None:
        qc = math.floor(strop / krok + 1e-9)
        qt = np.minimum(qt, qc); qp = np.minimum(qp, qc)
    return qt.astype(np.int32), qp.astype(np.int32), qc

def tnn_dekod(qt, Pp_t, krok, qc):
    v = qt.astype(np.float32) * np.float32(krok)
    if qc is not None:
        v = np.where(qt >= qc, np.maximum(Pp_t, np.float32(qc * krok)), v)
    return v.astype(np.float32)

def zz_bajty(r):
    z = ((r << 1) ^ (r >> 31)).astype(np.uint32)
    return (z & 255).astype(np.uint8), ((z >> 8) & 255).astype(np.uint8)


# =============================== STATISTIKY ===============================
class Stat:
    def __init__(self):
        self.n = np.zeros(3); self.s = np.zeros((3, NF)); self.sq = np.zeros((3, NF)); self.sa = np.zeros((3, NF))
        self.mx = np.zeros((3, NF)); self.h = np.zeros((3, NF, NBIN), np.int64)
        self.kc = {k: np.zeros((3, len(v[1]))) for k, v in KLS.items()}   # počet, součet chyb, součet |chyb|
        self.ks = {k: np.zeros((3, len(v[1]))) for k, v in KLS.items()}
        self.ka = {k: np.zeros((3, len(v[1]))) for k, v in KLS.items()}

    def add(self, E, Pt, band):
        for b in range(3):
            m = band == b
            if not m.any(): continue
            Em = E[m]; Tm = Pt[m]
            self.n[b] += Em.shape[0] * Em.shape[1] * Em.shape[2]
            for j in range(NF):
                e = Em[..., j].ravel()
                a = np.abs(e)
                self.s[b, j] += float(e.sum(dtype=np.float64)); self.sq[b, j] += float((e * e).sum(dtype=np.float64))
                self.sa[b, j] += float(a.sum(dtype=np.float64)); self.mx[b, j] = max(self.mx[b, j], float(a.max()))
                self.h[b, j] += np.bincount(np.minimum((a / BW + 0.5).astype(np.int32), NBIN - 1), minlength=NBIN)
            for k, (hr, lab) in KLS.items():
                j = IX[k]
                cls = np.digitize(Tm[..., j].ravel(), hr)
                e = Em[..., j].ravel()
                self.kc[k][b] += np.bincount(cls, minlength=len(lab))
                self.ks[k][b] += np.bincount(cls, weights=e, minlength=len(lab))
                self.ka[k][b] += np.bincount(cls, weights=np.abs(e), minlength=len(lab))

    def vysl(self, bands=(0, 1, 2)):
        b = list(bands); n = max(self.n[b].sum(), 1)
        cum = np.cumsum(self.h[b].sum(0), axis=1)
        q = lambda f: np.minimum((cum < f * n).sum(1), NBIN - 1) * BW
        return {"rmse": np.sqrt(self.sq[b].sum(0) / n), "mae": self.sa[b].sum(0) / n, "bias": self.s[b].sum(0) / n,
                "p95": q(0.95), "p99": q(0.99), "max": self.mx[b].max(0), "n": n}

    def trida(self, k):
        c = self.kc[k].sum(0); cc = np.maximum(c, 1)
        return c / max(c.sum(), 1), self.ks[k].sum(0) / cc, self.ka[k].sum(0) / cc


def ann_hodnota(P, j, mode):
    if mode == "frost": return (P[..., j] > 0.5).sum(2)
    return rocni(P, j, mode)


# =============================== VYHODNOCENÍ JEDNOHO BLOKU TESTOVACÍCH BUNĚK ===============================
def vyhodnot(Pp, Bfz, Pt, jih, band, ST, ANN, TNN):
    """Pp: dekódované spojité hodnoty před vynucením (n,NY,36,NF), Bfz: základ, Pt: pravda."""
    n, ny_, _, _ = Pt.shape
    j = IX["tnn"]
    for (krok, strop), d in TNN.items():
        qt, qp, _qc = tnn_kvant(Pt[..., j], Pp[..., j], krok, strop)
        lo, hi = zz_bajty((qt - qp).astype(np.int32))
        d["lo"].append(lo.reshape(n, ny_ * 36)); d["hi"].append(hi.reshape(n, ny_ * 36))
    qt, qp, qc = tnn_kvant(Pt[..., j], Pp[..., j], KROK, STROP)
    tnn_dek = tnn_dekod(qt, Pp[..., j], KROK, qc)
    for k, (nm, dek, vyp) in enumerate(SC):
        R = vynut2(Pp, Bfz, Pt, jih, tnn_dek if dek else None, vyp)
        E = R - Pt
        ST[k].add(E, Pt, band)
        for nm_, f, mode in ANNITEMS:
            a_t = ann_hodnota(Pt, IX[f], mode); a_r = ann_hodnota(R, IX[f], mode)
            ANN[k].setdefault(nm_, []).append((a_r - a_t).astype(np.float32))
        del R, E


# =============================== HLAVNÍ PROGRAM ===============================
def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    soubory = oi.najdi_soubory(oi.SCENAR)
    nper = len(soubory["fd"])
    p("=" * 78); p(f"ODHAD VQ3 – varianta '{VAR}' + TNn po dekádách (K = {K_ANOM})"); p("=" * 78)

    def otevri(i): return {k: Dataset(soubory[k][i]) for k in oi.PREFIXY}
    def zavri(ds):
        for d in ds.values(): d.close()

    ds = otevri(0)
    v0 = oi.najdi_promennou(ds["tg"], "TG"); o0 = oi.osy(v0)
    lat = np.asarray(ds["tg"].variables[v0.dimensions[o0["lat"]]][:], dtype=np.float64)
    lon = np.asarray(ds["tg"].variables[v0.dimensions[o0["lon"]]][:], dtype=np.float64)
    flip = lat[0] < lat[-1]
    lat_nu = lat[::-1] if flip else lat
    NCOLS = len(lon); radku = int(np.sum(lat_nu >= oi.LAT_MIN - 1e-6))
    tg0 = oi.cti_blok(v0, o0, 0, 1, flip, radku, False)[0]
    zavri(ds)
    maska = np.isfinite(tg0)
    T = oi.TILE
    rows_l, cols_l = [], []
    for ty in range(math.ceil(radku / T)):
        for tx in range(math.ceil(NCOLS / T)):
            blok = maska[ty * T:(ty + 1) * T, tx * T:(tx + 1) * T]
            r_, c_ = np.nonzero(blok)
            rows_l.append(ty * T + r_); cols_l.append(tx * T + c_)
    rows_all = np.concatenate(rows_l); cols_all = np.concatenate(cols_l); NC = len(rows_all)

    nt_soubor, lim = [], []
    for i in range(nper):
        ds = otevri(i)
        nt = len(ds["fd"].dimensions["time"])
        for k, d in ds.items():
            if len(d.dimensions["time"]) != nt: sys.exit("Chyba: různý počet kroků mezi indikátory")
        nt_soubor.append(nt); lim.append((min(nt, ROKY * 36) // 36) * 36); zavri(ds)
    ny = [l // 36 for l in lim]; NY = sum(ny); NT = NY * 36
    NM_PLNE = sum(nt_soubor) // 3; NT_PLNE = sum(nt_soubor); NY_PLNE = NT_PLNE // 36

    ns = min(NS, NC // 3); ntr = min(NTR, NC - ns)
    perm = rng.permutation(NC)
    S = perm[:ns]; TR = perm[ns:ns + ntr]
    use = np.concatenate([S, TR]); n_use = len(use)
    jihS = lat_nu[rows_all[S]] < 0
    bandS = pasmo(lat_nu[rows_all[S]])
    rows_use, cols_use = rows_all[use], cols_all[use]
    ns_tr = int(ns * 0.7); n_te = ns - ns_tr
    p(f"  buněk pevniny: {NC}; vzorek {ns} (trénink {ns_tr}, test {n_te}); buňky pro trénink základu: {ntr}")
    p(f"  čteno: {NY} let ({NT} dekád); varianta: {VAR}")
    p(f"  TNn po dekádách: krok {KROK:g} °C, strop {'žádný' if STROP is None else '%g °C' % STROP}; vypouštěné roční ukazatele (S2): {', '.join(VYP) if VYP else 'žádné'}")

    # ---------- čtení ----------
    p(""); p("[1/4] Čtení dat")
    print()
    Xs = np.zeros((ns, NT, NF), np.float16)
    base = np.zeros((n_use, nper, 36, NF), np.float32)
    g = 0; kel = {}; tr0 = time.time()
    for i in range(nper):
        ds = otevri(i)
        var = {k: oi.najdi_promennou(ds[k], oi.PREFIXY[k]) for k in oi.PREFIXY}
        o = {k: oi.osy(var[k]) for k in oi.PREFIXY}
        por = oi.poradi_su(ds["su"], var["su"], o["su"], nt_soubor[i], flip, radku)
        for a in range(0, lim[i], oi.BLOK):
            b = min(a + oi.BLOK, lim[i]); nb = b - a
            def rd(k, nez): return oi.cti_blok(var[k], o[k], a, b, flip, radku, nez)[..., rows_use, cols_use]
            F = np.empty((nb, n_use, NF), np.float32)
            for j, k in enumerate(("tg", "tn", "tx", "tnn", "txx")):
                x = rd(k, False)
                if (i, k) not in kel: kel[(i, k)] = oi.je_kelvin(var[k], x)
                if kel[(i, k)]: x = x - np.float32(273.15)
                F[:, :, j] = warp(np.nan_to_num(x, nan=0.0), k)
            for k in ("fd", "tr", "rr1", "r20"):
                F[:, :, IX[k]] = warp(np.nan_to_num(rd(k, True), nan=0.0), k)
            su = rd("su", True)[:, por, :]
            for j, k in enumerate(("su20", "su25", "su30", "su35")):
                F[:, :, IX[k]] = warp(np.nan_to_num(su[:, j, :], nan=0.0), k)
            F[:, :, IX["sdii"]] = warp(np.nan_to_num(rd("sdii", True), nan=0.0), "sdii")
            for t in range(nb): base[:, i, (a + t) % 36, :] += F[t]
            Xs[:, g:g + nb, :] = F[:, :ns, :].transpose(1, 0, 2).astype(np.float16)
            g += nb
            oi.progres(g, NT, tr0, "dekády")
        zavri(ds)
        base[:, i] /= ny[i]
    Btr = base[ns:]; Bs = base[:ns]

    # ---------- základ ----------
    p(""); p("[2/4] Trénink kódové knihy základu")
    V = Btr.reshape(-1, NF)
    vz = V[rng.choice(len(V), size=min(TRAIN_VZORKU_ZAKL1, len(V)), replace=False)]
    C1 = kmeans(vz, K1, KMEANS_ITER, rng)
    lab1 = prirad(V, C1)[0]; res = V - C1[lab1]
    C2 = np.zeros((K1, K2, NF), np.float32)
    for k in range(K1):
        m = lab1 == k; n = int(m.sum())
        if n: C2[k] = doplnit(kmeans(res[m], min(K2, max(1, n // 3)), 8, rng), K2)
    def kod_zaklad(Vx):
        k1 = prirad(Vx, C1)[0]; rs = Vx - C1[k1]; k2 = np.zeros(len(Vx), np.int32)
        for k in np.unique(k1):
            m = k1 == k; k2[m] = prirad(rs[m], C2[k])[0]
        return k1.astype(np.uint8), k2.astype(np.uint8), C1[k1] + C2[k1, k2]
    p1, p2, Vrec = kod_zaklad(Bs.reshape(-1, NF))
    Brec = Vrec.reshape(ns, nper, 36, NF)
    p1 = p1.reshape(ns, nper * 36); p2 = p2.reshape(ns, nper * 36)
    print(f"  hotovo za {oi.cas(time.time() - t0)}")

    def bexp(Barr):
        n = Barr.shape[0]
        return np.concatenate([np.tile(Barr[:, i][:, None], (1, ny[i], 1, 1)).reshape(n, ny[i] * 36, NF) for i in range(nper)], axis=1)
    def anomalie(c0, c1, cols):
        n = c1 - c0
        X = Xs[c0:c1].astype(np.float32).reshape(n, NY, 36, NF)
        A = np.empty_like(X); off = 0
        for i in range(nper):
            A[:, off:off + ny[i]] = X[:, off:off + ny[i]] - Bs[c0:c1, i][:, None]; off += ny[i]
        A = A.reshape(n, NY, 12, 42)
        Ac = A[..., cols]
        sig = np.empty((n, nper, 12), np.float32); off = 0
        for i in range(nper):
            sig[:, i] = np.sqrt((Ac[:, off:off + ny[i]] ** 2).mean(axis=(1, 3))); off += ny[i]
        sc = sig_kod(sig); sq = sig_dek(sc)
        sqx = np.concatenate([np.repeat(sq[:, i][:, None, :], ny[i], axis=1) for i in range(nper)], axis=1)[..., None]
        return A / sqx, sc, sqx

    # ---------- roční přesné hodnoty a termíny mrazů (pro velikosti) ----------
    p(""); p("[3/4] Roční přesné hodnoty (pro odhad velikosti)")
    ann_kod = {}; m_ls = []; m_fa = []
    tbs = time.time()
    chunks = [(c, min(c + CHUNK_BUNEK, ns)) for c in range(0, ns, CHUNK_BUNEK)]
    for ci, (c0, c1) in enumerate(chunks):
        n = c1 - c0
        Wt = Xs[c0:c1].astype(np.float32).reshape(n, NY, 36, NF)
        Pa = fyz(Wt); Pa[..., 5:13] = np.clip(np.rint(Pa[..., 5:13]), 0, 11)
        for nm_, f, mode in ANN_EX: ann_kod.setdefault(nm_, []).append(rocni(Pa, IX[f], mode))
        ls_, fa_, _o = mraz_kody(Pa[..., IX["fd"]], jihS[c0:c1]); m_ls.append(ls_); m_fa.append(fa_)
        oi.progres(ci + 1, len(chunks), tbs, "buňky")

    # ---------- varianta: kódová kniha indexu + vyhodnocení scénářů ----------
    nb_idx = max(1, math.ceil(math.log2(K_ANOM)))
    feats = CONFIGS[VAR]
    cols = [d * NF + IX[f] for d in range(3) for f in feats]
    p(""); p(f"[4/4] Varianta '{VAR}': kódováno {len(feats)} hodnot ({len(cols)} čísel), vyhodnocení 3 scénářů")
    tv = time.time()
    trch = [(c, min(c + CHUNK_BUNEK, ns_tr)) for c in range(0, ns_tr, CHUNK_BUNEK)]
    kvota = math.ceil(TRAIN_VZORKU_ANOM / max(len(trch), 1))
    sber = []
    for c0, c1 in trch:
        Zn, _sc, _sq = anomalie(c0, c1, cols)
        Zf = Zn.reshape(-1, 42)
        sber.append(Zf[rng.choice(len(Zf), size=min(kvota, len(Zf)), replace=False)])
    Zs = np.vstack(sber); del sber
    Csub = kmeans(Zs[:, cols], K_ANOM, KMEANS_ITER, rng)
    lab = prirad(Zs[:, cols], Csub)[0]
    cnt = np.bincount(lab, minlength=K_ANOM)
    CA = np.zeros((K_ANOM, 42), np.float32)
    for j in range(42): CA[:, j] = np.bincount(lab, weights=Zs[:, j], minlength=K_ANOM) / np.maximum(cnt, 1)
    for k in np.nonzero(cnt == 0)[0]: CA[k, cols] = Csub[k]
    del Zs
    idx_all = np.zeros((ns, NY, 12), np.uint16); sig_all = np.zeros((ns, nper, 12), np.uint8)

    ST = [Stat() for _ in SC]; ANN = [dict() for _ in SC]
    kombinace = [(0.1, None), (0.1, 10.0), (0.5, None), (0.5, 10.0), (1.0, 10.0), (KROK, STROP)]
    TNN = {}
    for kb in kombinace: TNN.setdefault(kb, {"lo": [], "hi": []})
    for ci, (c0, c1) in enumerate(chunks):
        Zn, sc, sqx = anomalie(c0, c1, cols)
        n = c1 - c0
        idx = prirad(Zn.reshape(-1, 42)[:, cols], Csub)[0].reshape(n, NY, 12)
        idx_all[c0:c1] = idx; sig_all[c0:c1] = sc
        te0 = max(c0, ns_tr)
        if te0 >= c1:
            oi.progres(ci + 1, len(chunks), tv, "buňky"); continue
        s0 = te0 - c0
        Arec = (CA[idx[s0:]] * sqx[s0:]).reshape(-1, NY, 36, NF)
        Br = bexp(Brec[te0:c1]).reshape(-1, NY, 36, NF)
        Wte = Xs[te0:c1].astype(np.float32).reshape(-1, NY, 36, NF)
        Pt = fyz(Wte); Pt[..., 5:13] = np.clip(np.rint(Pt[..., 5:13]), 0, 11)
        Pp = fyz(Br + Arec); Bfz = fyz(Br)
        vyhodnot(Pp, Bfz, Pt, jihS[te0:c1], bandS[te0:c1], ST, ANN, TNN)
        del Arec, Br, Wte, Pt, Pp, Bfz
        oi.progres(ci + 1, len(chunks), tv, "buňky")

    ix2 = idx_all.reshape(ns, NY * 12).astype(np.int32)
    by = gz_chunky((ix2 & 255).astype(np.uint8))[0]
    for bit in range(8, nb_idx):
        by += gz_chunky(np.packbits(((ix2 >> bit) & 1).astype(np.uint8), axis=1), (None,))[0]
    mb_idx = by / (ns * NY * 12) * NC * NM_PLNE / 1e6
    bs_ = gz_chunky(sig_all.reshape(ns, nper * 12), (None,))[0]
    mb_sig = bs_ / (ns * nper * 12) * NC * nper * 12 / 1e6
    H = entropie(idx_all, K_ANOM)

    # ---------- velikosti ----------
    def sz(a2, lagy): return gz_chunky(a2, lagy)[0]
    z1 = sz(p1, (None,)); z2 = sz(p2, (None,))
    mb_zak = (z1 + z2) / (2 * ns * nper * 36) * NC * nper * 36 * 2 / 1e6
    mb_knz = (K1 * NF * 2 + K1 * K2 * NF * 2) / 1e6
    mb_ann = {}
    for nm_, f, mode in ANN_EX:
        A_ = np.vstack(ann_kod[nm_])
        if mode in ("min", "max", "mean"):
            by = sz(np.clip(np.rint((A_ + 75.0) / 0.5), 0, 254).astype(np.uint8), (None, 1))
        else:
            v = np.clip(np.rint(A_), 0, 511).astype(np.int32)
            by = sz((v & 255).astype(np.uint8), (None, 1)) + sz(np.packbits(((v >> 8) & 1).astype(np.uint8), axis=1), (None,))
        mb_ann[nm_] = by / (ns * NY) * NC * NY_PLNE / 1e6
    mb_mraz = (sz(np.vstack(m_ls), (None, 1)) + sz(np.vstack(m_fa), (None, 1))) / (ns * NY) * NC * NY_PLNE / 1e6
    mb_tnn = {}
    for kb, d in TNN.items():
        lo = np.vstack(d["lo"]); hi = np.vstack(d["hi"])
        mb_tnn[kb] = (sz(lo, (None,)) + sz(hi, (None,))) / (n_te * NT) * NC * NT_PLNE / 1e6
    ann_tot = sum(mb_ann.values())
    mb_s0 = mb_zak + ann_tot + mb_mraz + mb_knz + mb_idx + mb_sig
    mb_s1 = mb_s0 - mb_ann["TNn"] + mb_tnn[(KROK, STROP)]
    mb_s2 = mb_s1 - sum(mb_ann[ROCNI_JMENO[x]] for x in VYP)

    # ---------- výpis ----------
    jn = [f"{nm} – {txt}" for nm, txt in zip((s[0] for s in SC), (
        "bezExt jako dosud (roční TNn přesně)",
        f"+ TNn po dekádách (krok {KROK:g} °C, strop {'žádný' if STROP is None else '%g °C' % STROP}), roční TNn vypuštěn",
        "S1 + vypuštěno ročně: " + (", ".join(VYP) if VYP else "nic")))]
    p(""); p("=" * 78); p("DEKÁDNÍ CHYBY V JEDNOTKÁCH (dny / °C / mm za den), test na %d buňkách mimo trénink, celý svět" % n_te); p("=" * 78)
    p("  bias = průměr (odhad − skutečnost): kladný = nadhodnocuje, záporný = podhodnocuje")
    p(f"  P95, P99 = chyba, kterou překročí jen 5 % resp. 1 % dekád (přesnost {BW:.2f}); max = největší chyba ve vzorku")
    for k, (nm, dek, vyp) in enumerate(SC):
        r = ST[k].vysl()
        p(""); p(f"  {jn[k]}")
        p(f"  {'indikátor':<11}{'RMSE':>8}{'MAE':>8}{'bias':>8}{'P95':>8}{'P99':>8}{'max':>9}")
        for j, f in enumerate(FEAT):
            p(f"  {POPIS[f]:<11}{r['rmse'][j]:>8.2f}{r['mae'][j]:>8.2f}{r['bias'][j]:>+8.2f}{r['p95'][j]:>8.2f}{r['p99'][j]:>8.2f}{r['max'][j]:>9.2f}")

    p(""); p("=" * 78); p("PODLE ZEMĚPISNÉ ŠÍŘKY – scénář S2 (RMSE / bias / P95)"); p("=" * 78)
    rb = [ST[2].vysl((b,)) for b in range(3)]
    p(f"  {'indikátor':<11}" + "".join(f"{x:>26}" for x in ("tropy", "mírné pásmo", "chladné oblasti")))
    for j, f in enumerate(FEAT):
        p(f"  {POPIS[f]:<11}" + "".join(f"{rb[b]['rmse'][j]:>9.2f}/{rb[b]['bias'][j]:+6.2f}/{rb[b]['p95'][j]:>6.2f}" for b in range(3)))
    p(""); p("  TNn a TXx podle pásma, všechny scénáře (RMSE / bias / P95):")
    for f in ("tnn", "txx"):
        for k, (nm, dek, vyp) in enumerate(SC):
            rr = [ST[k].vysl((b,)) for b in range(3)]; j = IX[f]
            p(f"  {POPIS[f]:<9}{nm}  " + "  ".join(f"{rr[b]['rmse'][j]:5.2f}/{rr[b]['bias'][j]:+5.2f}/{rr[b]['p95'][j]:5.2f}" for b in range(3)))
    p(f"  (pásma: {PASMA[0]} | {PASMA[1]} | {PASMA[2]}; počet dekád v pásmech: " + ", ".join(f"{int(ST[0].n[b])}" for b in range(3)) + ")")

    p(""); p("=" * 78); p("BIAS EXTRÉMŮ PODLE SKUTEČNÉ HODNOTY (záporný bias u TNn = odhad chladnější než skutečnost)"); p("=" * 78)
    if STROP is not None: p(f"  (S1/S2: TNn nad {STROP:g} °C se neukládá, tam se použije odhad z indexu)")
    for f in ("tnn", "txx"):
        p(f"  {POPIS[f]}: třída skutečné hodnoty | podíl dekád | bias S0 / S1 / S2 | MAE S0 / S1 / S2")
        sh = None; bi = []; ma = []
        for k in range(len(SC)):
            sh, b_, m_ = ST[k].trida(f); bi.append(b_); ma.append(m_)
        for c, lab in enumerate(KLS[f][1]):
            p(f"    {lab:<12}{sh[c] * 100:>6.1f} %   " + " / ".join(f"{bi[k][c]:+6.2f}" for k in range(3)) + "    " + " / ".join(f"{ma[k][c]:6.2f}" for k in range(3)))

    p(""); p("=" * 78); p("ROČNÍ UKAZATELE – chyba v jednotkách (MAE / P95 / bias), za rok a buňku"); p("=" * 78)
    p(f"  {'ukazatel':<17}{'jedn.':<7}" + "".join(f"{s[0]:>24}" for s in SC))
    for nm_, f, mode in ANNITEMS:
        cells = ""
        for k in range(len(SC)):
            e = np.concatenate(ANN[k][nm_], axis=0).ravel()
            cells += f"{np.abs(e).mean():>8.2f}/{np.percentile(np.abs(e), 95):>6.2f}/{e.mean():+6.2f}"
        p(f"  {nm_:<17}{JEDN.get(nm_, 'dny'):<7}" + cells)

    p(""); p("=" * 78); p("VELIKOSTI PO GZIPU (MB, odhad pro všech %d buněk)" % NC); p("=" * 78)
    p(f"  základ {mb_zak:6.1f} | termíny mrazů {mb_mraz:4.1f} | kódová kniha základu {mb_knz:3.1f} | index {mb_idx:6.1f} MB (entropie {H:.2f} z {nb_idx} bitů) + sigma {mb_sig:.1f}")
    p("  roční přesné proudy (MB):  " + ", ".join(f"{k} {v:.1f}" for k, v in sorted(mb_ann.items(), key=lambda x: -x[1])) + f"   celkem {ann_tot:.1f}")
    p("  TNn po dekádách (odhad; zbytek proti predikci z indexu):")
    for (krok, strop), v in sorted(mb_tnn.items(), key=lambda x: (x[0][0], x[0][1] is not None, x[0][1] or 0)):
        p(f"    krok {krok:<4g} °C, strop {'žádný' if strop is None else '%g °C' % strop:<8}: {v:6.1f} MB  → celkem S1 s tímto TNn {mb_s0 - mb_ann['TNn'] + v:6.1f} MB")
    p("")
    p(f"{'SCÉNÁŘ':<70}{'MB':>8}")
    p(f"{'S0 bezExt (roční TNn přesně)':<70}{mb_s0:>8.1f}")
    p(f"{'S1 + TNn po dekádách, roční TNn vypuštěn':<70}{mb_s1:>8.1f}")
    p(f"{'S2 + vypuštěno ročně: ' + (', '.join(VYP) if VYP else 'nic'):<70}{mb_s2:>8.1f}")
    vyp_all = sum(mb_ann[ROCNI_JMENO[x]] for x in POVOLENO)
    p(f"  rozdíl S2 proti S0: {mb_s2 - mb_s0:+.1f} MB (všechny povolené roční proudy dohromady: {vyp_all:.1f} MB)")
    p(f"  celkový čas {oi.cas(time.time() - t0)}")
    open("odhad_vq3.txt", "w", encoding="utf-8").write("\n".join(LOG) + "\n")
    print("  (výpis uložen do odhad_vq3.txt)")


if __name__ == "__main__":
    main()
