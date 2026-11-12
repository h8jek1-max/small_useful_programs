#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ODHAD VQ2 – které indikátory dát do indexu měsíce a co dopočítat; velikost "dělení" s bezeztrátovými počty

Jedním během (data se čtou jen jednou) změří několik VARIANT, jak sestavit index měsíce (K vzorů, 1 index na měsíc):

  vse         index vybírají všechny 14 hodnot (jako dosud, VQ+v)
  uzivatel6   index jen z TG, TN, TX, TNn, TXx, SDII  (VAŠE dělení; počty FD, TR, SU, RR1, R20mm se ukládají zvlášť)
  TX+TNn      index bez TN a TXx (TN a TXx se dopočítají z TX a TNn)
  bezExtremu  index bez TNn a TXx (dopočítají se z TN a TX; roční minimum/maximum jsou uložené přesně)
  jadro       index jen z TN, TX, RR1, R20mm (všechno ostatní se dopočítá)

Princip: index se vybírá jen podle "kódovaných" hodnot, ale kódová kniha uchovává průměr VŠECH 42 hodnot
v dané skupině, takže nekódované hodnoty se dopočítají jako podmíněný průměr (bez dalších bitů).
Roční přesné hodnoty (součty počtů, TG průměr, TNn minimum, TXx maximum, termíny mrazů) se pak VYNUTÍ.

Dále změří, kolik by zabraly počty FD, TR, RR1, R20mm, SU20..SU35 uložené BEZEZTRÁTOVĚ po dekádách a po měsících,
a z toho spočte velikost vašeho dělení (S1 = bezeztrátové dekády, S2 = bezeztrátové měsíční součty).

Vyžaduje odhad_indikatoru.py a odhad_vq.py ve stejné složce. Nic nezapisuje kromě odhad_vq2.txt.

    python odhad_vq2.py                         všechny varianty, K = 512
    python odhad_vq2.py --k 256
    python odhad_vq2.py --varianty vse,jadro    jen vybrané
    python odhad_vq2.py --roky 10               rychlejší zkouška
    python odhad_vq2.py --bunek 2000 --trenink 8000
"""
# (samostatná verze: nepotřebuje odhad_vq.py, jen odhad_indikatoru.py)
import sys, os, re, time, math
import numpy as np
from netCDF4 import Dataset
import odhad_indikatoru as oi

# =============================== NASTAVENÍ ===============================
NF = 14
FEAT = ["tg", "tn", "tx", "tnn", "txx", "fd", "tr", "rr1", "r20", "su20", "su25", "su30", "su35", "sdii"]
IX = {n: i for i, n in enumerate(FEAT)}
POPIS = {"tg": "TG °C", "tn": "TN °C", "tx": "TX °C", "tnn": "TNn °C", "txx": "TXx °C", "fd": "FD dny",
         "tr": "TR dny", "rr1": "RR1 dny", "r20": "R20mm dny", "su20": "SU20 dny", "su25": "SU25 dny",
         "su30": "SU30 dny", "su35": "SU35 dny", "sdii": "SDII mm/d"}

# pokřivení: (hranice, sklony) – sklon = kolik "jednotek důležitosti" připadá na °C resp. den resp. mm
WARP = {
    "tg":   ([-100, -10, 0, 5, 15, 100],      [0.6, 1.2, 2.0, 1.5, 0.6]),
    "tn":   ([-100, -12, -3, 3, 100],          [0.5, 1.0, 3.0, 0.7]),
    "tx":   ([-100, 20, 28, 38, 100],          [0.5, 1.0, 2.5, 1.0]),
    "tnn":  ([-100, -25, -12, -3, 3, 100],     [0.4, 1.0, 2.0, 3.0, 0.5]),
    "txx":  ([-100, 25, 33, 42, 100],          [0.4, 1.0, 2.5, 1.0]),
    "fd":   ([-1, 0, 1, 3, 400],               [1.0, 3.0, 2.0, 0.7]),
    "tr":   ([-1, 0, 1, 3, 400],               [1.0, 3.0, 2.0, 0.7]),
    "rr1":  ([-1, 400],                        [1.0]),
    "r20":  ([-1, 0, 1, 3, 400],               [1.0, 3.0, 2.0, 1.0]),
    "su20": ([-1, 0, 1, 5, 400],               [1.0, 1.5, 1.0, 0.6]),
    "su25": ([-1, 0, 1, 5, 400],               [1.0, 1.5, 1.0, 0.6]),
    "su30": ([-1, 0, 1, 5, 400],               [1.0, 2.5, 1.5, 0.8]),
    "su35": ([-1, 0, 1, 3, 400],               [1.0, 3.0, 2.0, 1.0]),
    "sdii": ([-1, 1, 5, 20, 400],              [1.0, 0.8, 0.4, 0.15]),
}
PRIORITA = {n: 1.0 for n in FEAT}     # např. PRIORITA["tnn"] = 1.5 zvýrazní zimovzdornost

K_ANOM = 256
K1, K2 = 256, 256
SIG_MIN = 0.05
TRAIN_VZORKU_ANOM = 300000
TRAIN_VZORKU_ZAKL1 = 300000
KMEANS_ITER = 12
CHUNK_BUNEK = 300
GZ_UROVEN = 9

def arg(n, v):
    if n in sys.argv:
        i = sys.argv.index(n)
        if i + 1 < len(sys.argv): return sys.argv[i + 1]
    return v
ROKY = int(arg("--roky", 30)); NS = int(arg("--bunek", 3000)); NTR = int(arg("--trenink", 10000))
SEED = int(arg("--seed", 42))
K_ANOM = int(arg("--k", K_ANOM))
LOG = []
def p(s=""):
    print(s); LOG.append(s)


# =============================== POKŘIVENÍ ===============================
_W = {}
for _n, (_h, _s) in WARP.items():
    _x = np.array(_h, float); _y = [_x[0]]
    for _k in range(len(_s)): _y.append(_y[-1] + _s[_k] * (_x[_k + 1] - _x[_k]))
    _W[_n] = (_x, np.array(_y))

def warp(x, n): return (np.interp(x, _W[n][0], _W[n][1]) * PRIORITA[n]).astype(np.float32)
def unwarp(y, n): return np.interp(np.asarray(y, np.float32) / PRIORITA[n], _W[n][1], _W[n][0]).astype(np.float32)
def fyz(W):
    P = np.empty(W.shape, np.float32)
    for j, n in enumerate(FEAT): P[..., j] = unwarp(W[..., j], n)
    return P


# =============================== K-MEANS ===============================
def prirad(X, C, chunk=40000):
    cn = (C * C).sum(1)
    lab = np.empty(len(X), np.int32); dist = np.empty(len(X), np.float32)
    for a in range(0, len(X), chunk):
        x = X[a:a + chunk]
        d = (x * x).sum(1)[:, None] - 2 * (x @ C.T) + cn[None, :]
        l = d.argmin(1)
        lab[a:a + chunk] = l
        dist[a:a + chunk] = np.maximum(d[np.arange(len(x)), l], 0)
    return lab, dist

def kmeans(X, K, iters, rng):
    X = np.ascontiguousarray(X, dtype=np.float32)
    N, D = X.shape
    K = max(1, min(K, N))
    C = np.empty((K, D), np.float32)
    C[0] = X[rng.integers(N)]
    d = X - C[0]; d2 = np.einsum("ij,ij->i", d, d)
    for k in range(1, K):
        s = d2.sum()
        i = rng.choice(N, p=d2 / s) if s > 0 else rng.integers(N)
        C[k] = X[i]
        d = X - C[k]; d2 = np.minimum(d2, np.einsum("ij,ij->i", d, d))
    for _ in range(iters):
        lab, _d = prirad(X, C)
        cnt = np.bincount(lab, minlength=K)
        for j in range(D):
            s = np.bincount(lab, weights=X[:, j], minlength=K)
            C[:, j] = np.where(cnt > 0, s / np.maximum(cnt, 1), C[:, j])
        prazdne = np.nonzero(cnt == 0)[0]
        if len(prazdne): C[prazdne] = X[rng.integers(N, size=len(prazdne))]
    return C

def doplnit(C, K):
    if len(C) >= K: return C
    return np.vstack([C, np.repeat(C[:1], K - len(C), axis=0)])


# =============================== SMĚRODATNÁ ODCHYLKA ===============================
def sig_kod(s): return np.clip(np.rint(np.log2(np.maximum(s, SIG_MIN) / SIG_MIN) * 8), 0, 255).astype(np.uint8)
def sig_dek(c): return (SIG_MIN * 2.0 ** (c.astype(np.float32) / 8)).astype(np.float32)


# =============================== ODHAD GZIP ===============================
def gz_chunky(a2, lagy=(None, 1, 12)):
    """a2: (buňky, kroky) uint8. Vrací (nejlepší bajty, název filtru) při dělení na bloky po 450 buňkách."""
    nejl = None
    for lg in lagy:
        s = 0
        for c0 in range(0, a2.shape[0], 450):
            x = np.ascontiguousarray(a2[c0:c0 + 450])
            s += len(oi.gz(np.ascontiguousarray(oi.filtr(x, None, lg, False)).tobytes(), GZ_UROVEN))
        if nejl is None or s < nejl[0]: nejl = (s, "raw" if lg is None else f"delta{lg}")
    return nejl

def entropie(idx, K=256):
    c = np.bincount(idx.ravel(), minlength=K).astype(np.float64); c = c[c > 0] / c.sum()
    return float(-(c * np.log2(c)).sum())

def uprav(P):
    """Vnutí fyzikální souvislosti a celočíselnost po dekódování."""
    P = P.copy()
    P[..., IX["tn"]] = np.minimum(P[..., IX["tn"]], P[..., IX["tg"]])
    P[..., IX["tx"]] = np.maximum(P[..., IX["tx"]], P[..., IX["tg"]])
    P[..., IX["tnn"]] = np.minimum(P[..., IX["tnn"]], P[..., IX["tn"]])
    P[..., IX["txx"]] = np.maximum(P[..., IX["txx"]], P[..., IX["tx"]])
    P[..., 5:13] = np.clip(np.rint(P[..., 5:13]), 0, 11)
    P[..., IX["r20"]] = np.minimum(P[..., IX["r20"]], P[..., IX["rr1"]])
    P[..., IX["su25"]] = np.minimum(P[..., IX["su25"]], P[..., IX["su20"]])
    P[..., IX["su30"]] = np.minimum(P[..., IX["su30"]], P[..., IX["su25"]])
    P[..., IX["su35"]] = np.minimum(P[..., IX["su35"]], P[..., IX["su30"]])
    P[..., IX["sdii"]] = np.where(P[..., IX["rr1"]] <= 0, 0, np.maximum(P[..., IX["sdii"]], 1.0))
    return P



# =============================== VYNUCENÍ PŘESNÝCH ROČNÍCH HODNOT ===============================
def mraz_kody(fd, jih):
    """fd (n,NY,36) = dny s mrazem v dekádě; jih (n,) bool. Vrací kódy (ls, fa) 0..18 a pořadí dekád.
    Rok se pro jižní polokouli posune o půl roku (dekáda 18 = začátek července)."""
    ordi = (np.arange(36)[None, :] + np.where(jih, 18, 0)[:, None]) % 36
    fs = np.take_along_axis(fd, ordi[:, None, :], axis=2)
    h1 = fs[..., :18] > 0; h2 = fs[..., 18:] > 0
    ls = np.where(h1.any(2), 18 - np.argmax(h1[..., ::-1], axis=2), 0)    # poslední jarní mráz: 1..18, 0 = žádný
    fa = np.where(h2.any(2), np.argmax(h2, axis=2) + 1, 0)                # první podzimní mráz: 1..18, 0 = žádný
    return ls.astype(np.uint8), fa.astype(np.uint8), ordi

def skaluj(Q, S, W, bazovy, cap, Bf):
    """Rozdělí roční součet S do 36 dekád podle vah Q (dekódovaný vzor), v povolených dekádách W,
    nad minimem bazovy a pod stropem cap; výsledek jsou celá čísla se součtem S (pokud to strop dovolí)."""
    zb = np.maximum(S.astype(np.float64) - bazovy.sum(2), 0)
    w = Q * W
    sw = w.sum(2)
    w = np.where((sw <= 0)[..., None], Bf * W, w); sw = w.sum(2)
    w = np.where((sw <= 0)[..., None], W.astype(np.float64), w); sw = w.sum(2)
    x = np.where(sw[..., None] > 0, w / np.maximum(sw, 1e-12)[..., None] * zb[..., None], 0.0)
    hlava = np.maximum(cap - bazovy, 0)
    for _ in range(3):
        x = np.minimum(x, hlava)
        vol = np.maximum(hlava - x, 0) * W
        sv = vol.sum(2)
        manko = np.minimum(zb - x.sum(2), sv)
        x = x + np.where(sv[..., None] > 0, vol / np.maximum(sv, 1e-12)[..., None] * manko[..., None], 0.0)
    x = np.minimum(x, hlava)
    fl = np.floor(x)
    rem = np.maximum(np.rint(zb - fl.sum(2)), 0)
    frac = np.where((fl < hlava) & W, x - fl, -1.0)
    rank = np.argsort(np.argsort(-frac, axis=2), axis=2)
    add = ((rank < rem[..., None]) & (frac >= 0)).astype(np.float64)
    return bazovy + fl + add

def kv05(x):
    """Roční teploty se ukládají s krokem 0,5 °C od -75 °C."""
    return np.clip(np.rint((x + 75.0) / 0.5), 0, 254) * 0.5 - 75.0

def vynut(P, Bf, Pt, jih):
    """P: dekódované spojité hodnoty (n,NY,36,NF), Bf: základní vzor, Pt: pravda (počty zaokrouhlené), jih: (n,) bool.
    Vynutí přesné roční součty počtů, roční průměr TG, roční minimum TNn, maximum TXx a termíny mrazů
    (v reálné databázi by se četly z uložených ročních hodnot)."""
    n, ny_, _, _ = P.shape
    R = P.copy(); Bf = np.maximum(Bf, 0)
    ones = np.ones((n, ny_, 36), bool); nuly = np.zeros((n, ny_, 36)); c11 = np.full((n, ny_, 36), 11.0)
    S_ = lambda f: Pt[..., IX[f]].sum(2)
    g = lambda f: np.maximum(P[..., IX[f]], 0)
    # --- FD: termíny mrazů
    ls, fa, ordi = mraz_kody(Pt[..., IX["fd"]], jih)
    inv = np.argsort(ordi, axis=1)
    s = np.arange(36)[None, None, :]
    lsi = ls.astype(int)[..., None] - 1
    fai = np.where(fa > 0, fa.astype(int) + 17, 36)[..., None]
    nula = ((s <= 17) & (s > lsi)) | ((s >= 18) & (s < fai))
    hr = ((s == lsi) & (lsi >= 0)) | ((s == fai) & (fai < 36))
    pick = lambda m: np.take_along_axis(m, inv[:, None, :], axis=2)
    R[..., IX["fd"]] = skaluj(g("fd"), S_("fd"), ~pick(nula), pick(hr).astype(np.float64), c11, Bf[..., IX["fd"]])
    # --- ostatní počty (řetězově kvůli SU20 >= SU25 >= SU30 >= SU35 a R20mm <= RR1)
    R[..., IX["rr1"]] = skaluj(g("rr1"), S_("rr1"), ones, nuly, c11, Bf[..., IX["rr1"]])
    R[..., IX["r20"]] = skaluj(g("r20"), S_("r20"), ones, nuly, R[..., IX["rr1"]], Bf[..., IX["r20"]])
    R[..., IX["tr"]] = skaluj(g("tr"), S_("tr"), ones, nuly, c11, Bf[..., IX["tr"]])
    cap = c11
    for f in ("su20", "su25", "su30", "su35"):
        R[..., IX[f]] = skaluj(g(f), S_(f), ones, nuly, cap, Bf[..., IX[f]])
        cap = R[..., IX[f]]
    # --- teploty
    tg_m = kv05(Pt[..., IX["tg"]].mean(2)); tnn_m = kv05(Pt[..., IX["tnn"]].min(2)); txx_m = kv05(Pt[..., IX["txx"]].max(2))
    delta = (tg_m - P[..., IX["tg"]].mean(2))[..., None]
    for f in ("tg", "tn", "tx", "tnn", "txx"): R[..., IX[f]] = P[..., IX[f]] + delta
    tnn = np.maximum(R[..., IX["tnn"]], tnn_m[..., None])
    np.put_along_axis(tnn, tnn.argmin(2)[..., None], tnn_m[..., None], axis=2)
    txx = np.minimum(R[..., IX["txx"]], txx_m[..., None])
    np.put_along_axis(txx, txx.argmax(2)[..., None], txx_m[..., None], axis=2)
    R[..., IX["tnn"]] = tnn; R[..., IX["txx"]] = txx
    tn = np.maximum(R[..., IX["tn"]], tnn); tx = np.minimum(R[..., IX["tx"]], txx)
    R[..., IX["tn"]] = tn; R[..., IX["tx"]] = tx
    R[..., IX["tg"]] = np.minimum(np.maximum(R[..., IX["tg"]], tn), tx)
    R[..., IX["sdii"]] = np.where(R[..., IX["rr1"]] <= 0, 0, np.maximum(P[..., IX["sdii"]], 1.0))
    return R

def rocni(P, j, mode):
    if mode == "sum": return P[..., j].sum(2)
    if mode == "min": return P[..., j].min(2)
    if mode == "max": return P[..., j].max(2)
    return P[..., j].mean(2)



import sys, os, re, time, math
import numpy as np
from netCDF4 import Dataset
import odhad_indikatoru as oi

CONFIGS = {
    "vse": list(FEAT),
    "uzivatel6": ["tg", "tn", "tx", "tnn", "txx", "sdii"],
    "TX+TNn": [f for f in FEAT if f not in ("tn", "txx")],
    "bezExtremu": [f for f in FEAT if f not in ("tnn", "txx")],
    "jadro": ["tn", "tx", "rr1", "r20"],
}
KRATKE = {"vse": "vse", "uzivatel6": "uziv6", "TX+TNn": "TX+TNn", "bezExtremu": "bezExt", "jadro": "jadro"}
COUNTS = ["fd", "tr", "rr1", "r20", "su20", "su25", "su30", "su35"]
ANN_EX = [("FD", "fd", "sum"), ("TR", "tr", "sum"), ("RR1", "rr1", "sum"), ("R20mm", "r20", "sum"),
          ("SU20", "su20", "sum"), ("SU25", "su25", "sum"), ("SU30", "su30", "sum"), ("SU35", "su35", "sum"),
          ("TNn", "tnn", "min"), ("TXx", "txx", "max"), ("TG", "tg", "mean")]
EXTRA = [("TN roční průměr", "tn"), ("TX roční průměr", "tx"), ("dekád s mrazem", "fd")]

K_ANOM = int(arg("--k", 512))   # výchozí 512 vzorů (9 bitů na měsíc)
VARIANTY = [x for x in arg("--varianty", ",".join(CONFIGS)).split(",") if x in CONFIGS]
LOG = []
def p(s=""):
    print(s); LOG.append(s)


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    soubory = oi.najdi_soubory(oi.SCENAR)
    nper = len(soubory["fd"])
    p("=" * 78); p(f"ODHAD VQ2 – které indikátory kódovat v indexu (K = {K_ANOM}), velikost dělení"); p("=" * 78)

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
    rows_use, cols_use = rows_all[use], cols_all[use]
    ns_tr = int(ns * 0.7)
    p(f"  buněk pevniny: {NC}; vzorek {ns} (trénink {ns_tr}, test {ns - ns_tr}); buňky pro trénink základu: {ntr}")
    p(f"  čteno: {NY} let ({NT} dekád); varianty: {', '.join(VARIANTY)}")

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
    def bexp_okna(Xw, okno):
        out = np.empty_like(Xw); off = 0
        for i in range(nper):
            y0 = off
            while y0 < off + ny[i]:
                y1 = min(y0 + okno, off + ny[i]); out[:, y0:y1] = Xw[:, y0:y1].mean(axis=1, keepdims=True); y0 = y1
            off += ny[i]
        return out
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

    # ---------- společný průchod vzorkem: bezeztrátové počty, roční přesné, referenční varianty ----------
    p(""); p("[3/4] Společná měření (bezeztrátové počty, roční přesné hodnoty, referenční varianty)")
    cnt_d = {f: [] for f in COUNTS}; cnt_m = {f: [] for f in COUNTS}
    ann_kod = {}; m_ls = []; m_fa = []
    REF = ["klim30", "Z", "Z10"]
    sse_ref = {k: np.zeros(NF) for k in REF}
    sx_t = np.zeros(NF); n_m = 0
    ext_ref = {k: {n: [0.0, 0] for n, _ in EXTRA} for k in ("Z", "Z10")}
    tbs = time.time()
    chunks = [(c, min(c + CHUNK_BUNEK, ns)) for c in range(0, ns, CHUNK_BUNEK)]
    for ci, (c0, c1) in enumerate(chunks):
        n = c1 - c0
        Wt = Xs[c0:c1].astype(np.float32).reshape(n, NY, 36, NF)
        Pa = fyz(Wt); Pa[..., 5:13] = np.clip(np.rint(Pa[..., 5:13]), 0, 11)
        for f in COUNTS:
            a = Pa[..., IX[f]]
            cnt_d[f].append(a.reshape(n, NT).astype(np.uint8))
            cnt_m[f].append(a.reshape(n, NY, 12, 3).sum(3).reshape(n, NY * 12).astype(np.uint8))
        for nm_, f, mode in ANN_EX:
            ann_kod.setdefault(nm_, []).append(rocni(Pa, IX[f], mode))
        ls_, fa_, _o = mraz_kody(Pa[..., IX["fd"]], jihS[c0:c1]); m_ls.append(ls_); m_fa.append(fa_)
        te0 = max(c0, ns_tr)
        if te0 < c1:
            s0 = te0 - c0
            Pt = Pa[s0:]; Wte = Wt[s0:]; jih = jihS[te0:c1]
            Bt = bexp(Bs[te0:c1]).reshape(-1, NY, 36, NF); Br = bexp(Brec[te0:c1]).reshape(-1, NY, 36, NF)
            Pbr = fyz(Br)
            out = {"klim30": fyz(Bt), "Z": vynut(Pbr, Pbr, Pt, jih)}
            P10 = fyz(bexp_okna(Wte, 10)); out["Z10"] = vynut(P10, P10, Pt, jih)
            for k in REF: sse_ref[k] += ((out[k] - Pt) ** 2).sum(axis=(0, 1, 2))
            sx_t += Pt.sum(axis=(0, 1, 2)); n_m += Pt.shape[0] * NY * 36
            for k in ("Z", "Z10"):
                for nm_, f in EXTRA:
                    j = IX[f]
                    if f == "fd": a_t = (Pt[..., j] > 0.5).sum(2); a_r = (out[k][..., j] > 0.5).sum(2)
                    else: a_t = Pt[..., j].mean(2); a_r = out[k][..., j].mean(2)
                    ext_ref[k][nm_][0] += float(np.abs(a_r - a_t).sum()); ext_ref[k][nm_][1] += a_t.size
        oi.progres(ci + 1, len(chunks), tbs, "buňky")

    # ---------- varianty indexu ----------
    res = {}
    nb_idx = max(1, math.ceil(math.log2(K_ANOM)))
    for vi, name in enumerate(VARIANTY):
        feats = CONFIGS[name]
        cols = [d * NF + IX[f] for d in range(3) for f in feats]
        p(""); p(f"[4/4] Varianta '{name}' ({vi + 1}/{len(VARIANTY)}): kódováno {len(feats)} hodnot ({len(cols)} čísel)")
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
        sse = np.zeros(NF); ext = {n: [0.0, 0] for n, _ in EXTRA}
        for ci, (c0, c1) in enumerate(chunks):
            Zn, sc, sqx = anomalie(c0, c1, cols)
            n = c1 - c0
            idx = prirad(Zn.reshape(-1, 42)[:, cols], Csub)[0].reshape(n, NY, 12)
            idx_all[c0:c1] = idx; sig_all[c0:c1] = sc
            te0 = max(c0, ns_tr)
            if te0 >= c1: continue
            s0 = te0 - c0
            Arec = (CA[idx[s0:]] * sqx[s0:]).reshape(-1, NY, 36, NF)
            Br = bexp(Brec[te0:c1]).reshape(-1, NY, 36, NF)
            Wte = Xs[te0:c1].astype(np.float32).reshape(-1, NY, 36, NF)
            Pt = fyz(Wte); Pt[..., 5:13] = np.clip(np.rint(Pt[..., 5:13]), 0, 11)
            R = vynut(fyz(Br + Arec), fyz(Br), Pt, jihS[te0:c1])
            sse += ((R - Pt) ** 2).sum(axis=(0, 1, 2))
            for nm_, f in EXTRA:
                j = IX[f]
                if f == "fd": a_t = (Pt[..., j] > 0.5).sum(2); a_r = (R[..., j] > 0.5).sum(2)
                else: a_t = Pt[..., j].mean(2); a_r = R[..., j].mean(2)
                ext[nm_][0] += float(np.abs(a_r - a_t).sum()); ext[nm_][1] += a_t.size
            oi.progres(ci + 1, len(chunks), tv, "buňky")
        ix2 = idx_all.reshape(ns, NY * 12).astype(np.int32)
        by = gz_chunky((ix2 & 255).astype(np.uint8))[0]
        for bit in range(8, nb_idx):
            by += gz_chunky(np.packbits(((ix2 >> bit) & 1).astype(np.uint8), axis=1), (None,))[0]
        mb_idx = by / (ns * NY * 12) * NC * NM_PLNE / 1e6
        bs_ = gz_chunky(sig_all.reshape(ns, nper * 12), (None,))[0]
        mb_sig = bs_ / (ns * nper * 12) * NC * nper * 12 / 1e6
        res[name] = {"rmse": np.sqrt(sse / n_m), "ext": {n: e[0] / max(e[1], 1) for n, e in ext.items()},
                     "mb_idx": mb_idx, "mb_sig": mb_sig, "H": entropie(idx_all, K_ANOM)}

    # ---------- velikosti společných částí ----------
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
    mb_cd = {f: sz(np.vstack(cnt_d[f]), (None, 1, 36)) / (ns * NT) * NC * NT_PLNE / 1e6 for f in COUNTS}
    mb_cm = {f: sz(np.vstack(cnt_m[f]), (None, 1, 12)) / (ns * NY * 12) * NC * NM_PLNE / 1e6 for f in COUNTS}

    # ---------- výpis ----------
    jm = [KRATKE[n] for n in VARIANTY]
    p(""); p("=" * 78); p("DEKÁDNÍ HODNOTY – RMSE ve fyzikálních jednotkách (menší = lépe), test na buňkách mimo trénink"); p("=" * 78)
    p("  klim = jen 30letý průměr dekády; Z = základ + přesné roční (bez indexu); ostatní sloupce = index měsíce + přesné roční")
    p("  hvězdička = u varianty uziv6 se tato hodnota ukládá bezeztrátově zvlášť (chyba 0)")
    hl = f"{'indikátor':<11}{'klim':>7}{'Z':>7}" + "".join(f"{x:>9}" for x in jm)
    p(hl)
    for j, f in enumerate(FEAT):
        rk = math.sqrt(sse_ref["klim30"][j] / n_m); rz = math.sqrt(sse_ref["Z"][j] / n_m)
        cells = ""
        for n in VARIANTY:
            if n == "uzivatel6" and f in COUNTS: cells += f"{'0.00*':>9}"
            else: cells += f"{res[n]['rmse'][j]:>9.2f}"
        p(f"{POPIS[f]:<11}{rk:>7.2f}{rz:>7.2f}{cells}")
    p("")
    p("Roční ukazatele, které se přesně neukládají (MAE za rok):")
    p(f"{'ukazatel':<17}{'Z':>7}{'Z10':>7}" + "".join(f"{x:>9}" for x in jm))
    for nm_, _ in EXTRA:
        v = {k: ext_ref[k][nm_][0] / max(ext_ref[k][nm_][1], 1) for k in ext_ref}
        cells = "".join(f"{res[n]['ext'][nm_]:>9.2f}" if not (n == "uzivatel6" and nm_ == "dekád s mrazem") else f"{'0.00*':>9}" for n in VARIANTY)
        p(f"{nm_:<17}{v['Z']:>7.2f}{v['Z10']:>7.2f}{cells}")

    p(""); p("=" * 78); p("VELIKOSTI PO GZIPU (MB, odhad pro všech %d buněk)" % NC); p("=" * 78)
    ann_tot = sum(mb_ann.values()); ann_temp = mb_ann["TG"] + mb_ann["TNn"] + mb_ann["TXx"]
    p(f"  základ (2 B / dekádu roku / období) {mb_zak:6.1f} | přesné roční (11 proudů) {ann_tot:5.1f} | termíny mrazů {mb_mraz:4.1f} | kódová kniha základu {mb_knz:3.1f}")
    p("  index měsíce a směrodatné odchylky podle varianty:")
    for n in VARIANTY: p(f"    {n:<11} index {res[n]['mb_idx']:6.1f} MB (entropie {res[n]['H']:.2f} z {nb_idx} bitů) + sigma {res[n]['mb_sig']:.1f} MB")
    p("  počty BEZEZTRÁTOVĚ, po dekádách / po měsících (součty):")
    for f in COUNTS: p(f"    {f:<5} {mb_cd[f]:6.1f} MB / {mb_cm[f]:6.1f} MB")
    cd_tot = sum(mb_cd.values()); cm_tot = sum(mb_cm.values())
    p(f"    celkem {cd_tot:6.1f} MB / {cm_tot:6.1f} MB    (jen FD dekádně: {mb_cd['fd']:.1f} MB)")
    p("")
    p(f"{'VARIANTA':<66}{'MB':>8}")
    for n in VARIANTY:
        if n == "uzivatel6": continue
        h = mb_zak + ann_tot + mb_mraz + mb_knz + res[n]["mb_idx"] + res[n]["mb_sig"]
        p(f"{('H  index z ' + n + ' + přesné roční + mrazy'):<66}{h:>8.1f}")
        p(f"{('H+FD dekádně přesně (' + n + '; roční FD a mrazy lze vypustit)'):<66}{h + mb_cd['fd'] - mb_mraz - mb_ann['FD']:>8.1f}")
    if "uzivatel6" in res:
        r = res["uzivatel6"]; zk = mb_zak + mb_knz + r["mb_idx"] + r["mb_sig"] + ann_temp
        p(f"{'S1 vaše dělení: index z 6 hodnot + počty bezeztrátově po dekádách':<66}{zk + cd_tot:>8.1f}")
        p(f"{'S2 vaše dělení: index z 6 hodnot + počty bezeztrátově po měsících':<66}{zk + cm_tot:>8.1f}")
        p("  (u S1/S2 je základ započten celý, tj. mírně pesimisticky; roční TG, TNn, TXx přesně)")
    p(f"  celkový čas {oi.cas(time.time() - t0)}")
    open("odhad_vq2.txt", "w", encoding="utf-8").write("\n".join(LOG) + "\n")
    print("  (výpis uložen do odhad_vq2.txt)")


if __name__ == "__main__":
    main()
