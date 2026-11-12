#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STAVBA DATABÁZE agroklimatických indikátorů (varianta bezExt + TNn po dekádách + výjimky SDII)

Navazuje na odhad_vq3.py. Pracuje PO DÁVKÁCH, vypisuje procenta průběhu a po přerušení (Ctrl+C, vypnutí,
--max-minut) PŘI DALŠÍM SPUŠTĚNÍ POKRAČUJE tam, kde skončila. Stačí pořád spouštět stejný příkaz.

Fáze:
  A  jednorázové přečtení netCDF souborů do kompaktního dočasného úložiště (cca 4,3 GB na disku,
     v paměti jen jeden časový blok)
  B  výcvik kódových knih (základ + index měsíce) na vzorku buněk
  C  kódování všech buněk po dávkách (bloky po 450 buňkách) a zápis streamů do složky databáze
  D  zápis meta.json; volitelně smazání dočasného úložiště

Dekádní SDII se nikde neořezává. Skutečné hodnoty SDII >= --sdii-prah (výchozí 50 mm/den) se ukládají
přesně jako výjimky (pár tisíc záznamů), zbytek SDII je z indexu stejně jako v odhadu.
Teploty se ukládají jako TNn po dekádách (krok --krok, strop --strop), roční TNn se z něj dopočítá.
Nic se nevypouští: roční součty počtů, TG, TXx a termíny mrazů zůstávají přesné.

Vyžaduje odhad_vq2.py a odhad_indikatoru.py ve stejné složce. Výstup: složka databaze/ (--vystup).

    python stavba_db.py                      výchozí: krok 0.5 °C, strop 10 °C, dávka 2 bloky (900 buněk)
    python stavba_db.py --max-minut 60       po 60 minutách bezpečně skončí; příště pokračuje
    python stavba_db.py --davka 1            menší paměť (450 buněk najednou)
    python stavba_db.py --overit             po stavbě: dekóduje uložené kontrolní buňky a vypíše chyby
    python stavba_db.py --smazat-temp        po dokončení smaže dočasné úložiště (4,3 GB)
    python stavba_db.py --znovu              začne úplně od začátku (smaže vše ve složce databáze)
    python stavba_db.py --znovu-kodovani     zopakuje jen fázi C (např. po změně --sdii-prah)
    python stavba_db.py --roky 6 --blok 50   rychlá zkouška na pár letech
"""
import sys, os, time, math, json, zlib, types, shutil
import numpy as np
from netCDF4 import Dataset
import odhad_indikatoru as oi
from odhad_vq2 import (FEAT, IX, NF, POPIS, COUNTS, WARP, PRIORITA, SIG_MIN, K1, K2, K_ANOM, KMEANS_ITER,
                       TRAIN_VZORKU_ANOM, TRAIN_VZORKU_ZAKL1, SEED, CONFIGS,
                       arg, warp, fyz, prirad, doplnit, sig_kod, sig_dek, mraz_kody, skaluj, kv05, rocni)

# =============================== NASTAVENÍ ===============================
VYSTUP = arg("--vystup", "databaze")
CB = int(arg("--blok", 450))                 # buněk v jednom bloku (jednotka čtení v PHP)
DAVKA = int(arg("--davka", 2))               # bloků zpracovaných najednou
CH = int(arg("--chunk", 150))                # buněk zpracovaných najednou uvnitř bloku (paměť)
KROK = float(arg("--krok", 0.5))
_s = arg("--strop", "10")
STROP = None if _s.lower() in ("zadny", "none", "-") else float(_s)
SDII_PRAH = float(arg("--sdii-prah", 50))
NS_TR = int(arg("--vyber", 2100))            # buněk pro výcvik indexu
NTR_B = int(arg("--trenink", 10000))         # buněk pro výcvik základu
ROKY = int(arg("--roky", 1000))
MAXMIN = float(arg("--max-minut", 0))
MAXDAVEK = int(arg("--max-davek", 0))
GZ = int(arg("--gz", 9))
VARIANTA = "bezExtremu"
FLAG = lambda n: n in sys.argv
T0 = time.time()
NB_IDX = max(1, math.ceil(math.log2(K_ANOM)))
RANGE_SDII_MAX = 60000.0


class Pauza(Exception):
    pass


def cas(s):
    s = int(max(s, 0)); h, r = divmod(s, 3600); m, ss = divmod(r, 60)
    return f"{h}:{m:02d}:{ss:02d}" if h else f"{m}:{ss:02d}"

def p(s=""):
    print(s)
    try:
        with open(os.path.join(VYSTUP, "stavba_log.txt"), "a", encoding="utf-8") as f: f.write(s + "\n")
    except OSError:
        pass

def prog(h, c, h0, t0, popis):
    pct = 100.0 * h / c if c else 100.0
    el = time.time() - t0
    rate = (h - h0) / el if el > 0 and h > h0 else 0
    zb = (c - h) / rate if rate > 0 else 0
    sys.stdout.write(f"\r  {popis}: {pct:5.1f} %   uplynulo {cas(el)}, zbývá cca {cas(zb)}      ")
    sys.stdout.flush()
    if h >= c: sys.stdout.write("\n")

def cas_vyprsel():
    return MAXMIN > 0 and (time.time() - T0) / 60.0 >= MAXMIN

def cesta(n): return os.path.join(VYSTUP, n)


# =============================== STAV (pokračování) ===============================
def nacti_stav():
    if os.path.exists(cesta("stav.json")):
        with open(cesta("stav.json"), encoding="utf-8") as f: return json.load(f)
    return {}

def uloz_stav(stav):
    tmp = cesta("stav.json.tmp")
    with open(tmp, "w", encoding="utf-8") as f: json.dump(stav, f)
    os.replace(tmp, cesta("stav.json"))


# =============================== PŘÍPRAVA (mřížka, buňky, časová osa) ===============================
def priprava():
    G = types.SimpleNamespace()
    soubory = oi.najdi_soubory(oi.SCENAR)
    nper = len(soubory["fd"])
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
    G.rows = np.concatenate(rows_l); G.cols = np.concatenate(cols_l); G.NC = len(G.rows)
    nt_soubor, lim = [], []
    for i in range(nper):
        ds = otevri(i)
        nt = len(ds["fd"].dimensions["time"])
        for k, d in ds.items():
            if len(d.dimensions["time"]) != nt: sys.exit("Chyba: různý počet kroků mezi indikátory")
        nt_soubor.append(nt); lim.append((min(nt, ROKY * 36) // 36) * 36); zavri(ds)
    G.nper = nper; G.nt_soubor = nt_soubor; G.lim = lim
    G.ny = [l // 36 for l in lim]; G.NY = sum(G.ny); G.NT = G.NY * 36
    G.lat = lat_nu; G.lon = lon; G.flip = flip; G.radku = radku
    G.otevri = otevri; G.zavri = zavri
    G.jih = lat_nu[G.rows] < 0
    G.nblok = math.ceil(G.NC / CB)
    return G


def F_z_uloziste(Tm, Cm, Sm):
    """Dočasné úložiště -> pokřivené hodnoty (nb, n, NF)."""
    nb, n = Sm.shape
    F = np.empty((nb, n, NF), np.float32)
    for j, k in enumerate(("tg", "tn", "tx", "tnn", "txx")):
        F[:, :, j] = warp(Tm[..., j].astype(np.float32) * np.float32(0.01), k)
    for j, k in enumerate(COUNTS):
        F[:, :, IX[k]] = warp(Cm[..., j].astype(np.float32), k)
    F[:, :, IX["sdii"]] = warp(Sm.astype(np.float32), "sdii")
    return F

def otevri_uloziste(G, mode):
    NT, NC = G.NT, G.NC
    return (np.memmap(cesta("tmp_T.i16"), dtype=np.int16, mode=mode, shape=(NT, NC, 5)),
            np.memmap(cesta("tmp_C.u8"), dtype=np.uint8, mode=mode, shape=(NT, NC, 8)),
            np.memmap(cesta("tmp_S.f16"), dtype=np.float16, mode=mode, shape=(NT, NC)))


# =============================== FÁZE A: netCDF -> dočasné úložiště ===============================
def faze_A(G, stav):
    if stav.get("A_hotovo"): return
    p(""); p("[A] Čtení netCDF souborů do dočasného úložiště (jednorázově)")
    NT, NC = G.NT, G.NC
    novy = stav.get("A_t", 0) == 0
    mmT, mmC, mmS = otevri_uloziste(G, "w+" if novy else "r+")
    sd = stav.setdefault("sdii", {"max": 0.0, "n30": 0, "n50": 0, "n100": 0, "n200": 0, "n400": 0, "nfill": 0})
    g_hotovo = stav.get("A_t", 0); h0 = g_hotovo; t0 = time.time()
    off = 0; kel = {}
    for i in range(G.nper):
        if off + G.lim[i] <= g_hotovo:
            off += G.lim[i]; continue
        ds = G.otevri(i)
        var = {k: oi.najdi_promennou(ds[k], oi.PREFIXY[k]) for k in oi.PREFIXY}
        o = {k: oi.osy(var[k]) for k in oi.PREFIXY}
        por = oi.poradi_su(ds["su"], var["su"], o["su"], G.nt_soubor[i], G.flip, G.radku)
        for a in range(0, G.lim[i], oi.BLOK):
            b = min(a + oi.BLOK, G.lim[i]); nb = b - a; gg = off + a
            if gg + nb <= g_hotovo: continue
            def rd(k, nez): return oi.cti_blok(var[k], o[k], a, b, G.flip, G.radku, nez)[..., G.rows, G.cols]
            Tb = np.empty((nb, NC, 5), np.int16); Cb = np.empty((nb, NC, 8), np.uint8)
            for j, k in enumerate(("tg", "tn", "tx", "tnn", "txx")):
                x = rd(k, False)
                if (i, k) not in kel: kel[(i, k)] = oi.je_kelvin(var[k], x)
                if kel[(i, k)]: x = x - np.float32(273.15)
                x = np.nan_to_num(x, nan=0.0); x = np.where(np.abs(x) > 1e4, 0, x)
                Tb[:, :, j] = np.clip(np.rint(x * 100.0), -32000, 32000)
            for j, k in enumerate(("fd", "tr", "rr1", "r20")):
                x = np.nan_to_num(rd(k, True), nan=0.0); x = np.where((x > 1e4) | (x < 0), 0, x)
                Cb[:, :, j] = np.clip(np.rint(x), 0, 255)
            su = rd("su", True)[:, por, :]
            for j in range(4):
                x = np.nan_to_num(su[:, j, :], nan=0.0); x = np.where((x > 1e4) | (x < 0), 0, x)
                Cb[:, :, 4 + j] = np.clip(np.rint(x), 0, 255)
            sx = np.nan_to_num(rd("sdii", True), nan=0.0)
            bad = (sx > 1e5) | (sx < 0); sd["nfill"] += int(bad.sum()); sx = np.where(bad, 0, sx)
            sd["max"] = max(sd["max"], float(sx.max()))
            for nm, th in (("n30", 30), ("n50", 50), ("n100", 100), ("n200", 200), ("n400", 400)):
                sd[nm] += int((sx >= th).sum())
            mmT[gg:gg + nb] = Tb; mmC[gg:gg + nb] = Cb; mmS[gg:gg + nb] = np.clip(sx, 0, RANGE_SDII_MAX).astype(np.float16)
            for m_ in (mmT, mmC, mmS): m_.flush()
            g_hotovo = gg + nb; stav["A_t"] = g_hotovo; uloz_stav(stav)
            prog(g_hotovo, NT, h0, t0, "dekády")
            if cas_vyprsel() and g_hotovo < NT: G.zavri(ds); raise Pauza()
        G.zavri(ds); off += G.lim[i]
    stav["A_hotovo"] = True; uloz_stav(stav)
    p(f"  SDII ve zdroji: maximum {sd['max']:.1f} mm/den; počty dekád (buňka × dekáda) ≥30: {sd['n30']}, ≥50: {sd['n50']}, "
      f"≥100: {sd['n100']}, ≥200: {sd['n200']}, ≥400: {sd['n400']}; neplatných (výplňových) hodnot: {sd['nfill']}")


# =============================== K-MEANS S PROCENTY ===============================
def kmeans_p(X, K, iters, rng, popis):
    X = np.ascontiguousarray(X, dtype=np.float32)
    N, D = X.shape
    K = max(1, min(K, N))
    C = np.empty((K, D), np.float32)
    C[0] = X[rng.integers(N)]
    d = X - C[0]; d2 = np.einsum("ij,ij->i", d, d)
    t0 = time.time(); krok = max(1, K // 40)
    for k in range(1, K):
        s = d2.sum()
        i = rng.choice(N, p=d2 / s) if s > 0 else rng.integers(N)
        C[k] = X[i]
        d = X - C[k]; d2 = np.minimum(d2, np.einsum("ij,ij->i", d, d))
        if popis and k % krok == 0: prog(k, K, 0, t0, popis + " (start)")
    t0 = time.time()
    for it in range(iters):
        lab, _d = prirad(X, C, 20000)
        cnt = np.bincount(lab, minlength=K)
        for j in range(D):
            s = np.bincount(lab, weights=X[:, j], minlength=K)
            C[:, j] = np.where(cnt > 0, s / np.maximum(cnt, 1), C[:, j])
        prazdne = np.nonzero(cnt == 0)[0]
        if len(prazdne): C[prazdne] = X[rng.integers(N, size=len(prazdne))]
        if popis: prog(it + 1, iters, 0, t0, popis + " (iterace)")
    return C

def f16(a): return np.asarray(a, np.float32).astype(np.float16).astype(np.float32)


# =============================== FÁZE B: výcvik kódových knih ===============================
def faze_B(G, stav):
    if stav.get("B_hotovo") and os.path.exists(cesta("kody.npz")): return
    p(""); p("[B] Výcvik kódových knih (základ + index měsíce)")
    NT, NC, nper, NY, ny = G.NT, G.NC, G.nper, G.NY, G.ny
    rng = np.random.default_rng(SEED)
    perm = rng.permutation(NC)
    ns = min(NS_TR, NC // 3); ntr = min(NTR_B, NC - ns)
    S = perm[:ns]; TRc = perm[ns:ns + ntr]
    zaklad_hotov = os.path.exists(cesta("kody_zaklad.npz"))
    use = S if zaklad_hotov else np.concatenate([S, TRc])
    n_use = len(use)
    mmT, mmC, mmS = otevri_uloziste(G, "r")
    Xs = np.zeros((ns, NT, NF), np.float16)
    base = np.zeros((n_use, nper, 36, NF), np.float32)
    t0 = time.time(); off = 0
    p("  čtení vzorku z úložiště")
    for i in range(nper):
        for a in range(0, G.lim[i], oi.BLOK):
            b = min(a + oi.BLOK, G.lim[i]); nb = b - a; gg = off + a
            Tm = np.asarray(mmT[gg:gg + nb])[:, use, :]; Cm = np.asarray(mmC[gg:gg + nb])[:, use, :]; Sm = np.asarray(mmS[gg:gg + nb])[:, use]
            F = F_z_uloziste(Tm, Cm, Sm)
            for t in range(nb): base[:, i, (a + t) % 36, :] += F[t]
            Xs[:, gg:gg + nb, :] = F[:, :ns, :].transpose(1, 0, 2).astype(np.float16)
            prog(gg + nb, NT, 0, t0, "vzorek")
        off += G.lim[i]
        base[:, i] /= ny[i]
    Bs = base[:ns]
    # ---- základ
    if zaklad_hotov:
        z = np.load(cesta("kody_zaklad.npz")); C1 = z["C1"]; C2 = z["C2"]
    else:
        Btr = base[ns:]
        V = Btr.reshape(-1, NF)
        vz = V[rng.choice(len(V), size=min(TRAIN_VZORKU_ZAKL1, len(V)), replace=False)]
        C1 = kmeans_p(vz, K1, KMEANS_ITER, rng, "základ 1")
        lab1 = prirad(V, C1, 20000)[0]; res = V - C1[lab1]
        C2 = np.zeros((K1, K2, NF), np.float32); t1 = time.time()
        for k in range(K1):
            m = lab1 == k; n = int(m.sum())
            if n: C2[k] = doplnit(kmeans_p(res[m], min(K2, max(1, n // 3)), 8, rng, None), K2)
            prog(k + 1, K1, 0, t1, "základ 2")
        C1 = f16(C1); C2 = f16(C2)
        np.savez(cesta("kody_zaklad.npz"), C1=C1, C2=C2)
        del Btr, V, vz, res
    # ---- index měsíce (anomálie)
    feats = CONFIGS[VARIANTA]
    cols = [d * NF + IX[f] for d in range(3) for f in feats]
    def anomalie_tr(c0, c1):
        n = c1 - c0
        X = Xs[c0:c1].astype(np.float32).reshape(n, NY, 36, NF)
        A = np.empty_like(X); o2 = 0
        for i in range(nper):
            A[:, o2:o2 + ny[i]] = X[:, o2:o2 + ny[i]] - Bs[c0:c1, i][:, None]; o2 += ny[i]
        A = A.reshape(n, NY, 12, 42)
        Ac = A[..., cols]
        sig = np.empty((n, nper, 12), np.float32); o2 = 0
        for i in range(nper):
            sig[:, i] = np.sqrt((Ac[:, o2:o2 + ny[i]] ** 2).mean(axis=(1, 3))); o2 += ny[i]
        sq = sig_dek(sig_kod(sig))
        sqx = np.concatenate([np.repeat(sq[:, i][:, None, :], ny[i], axis=1) for i in range(nper)], axis=1)[..., None]
        return A / sqx
    trch = [(c, min(c + 300, ns)) for c in range(0, ns, 300)]
    kvota = math.ceil(TRAIN_VZORKU_ANOM / max(len(trch), 1)); sber = []; t2 = time.time()
    for ci, (c0, c1) in enumerate(trch):
        Zf = anomalie_tr(c0, c1).reshape(-1, 42)
        sber.append(Zf[rng.choice(len(Zf), size=min(kvota, len(Zf)), replace=False)])
        prog(ci + 1, len(trch), 0, t2, "anomálie vzorku")
    Zs = np.vstack(sber); del sber, Xs, base
    Csub = kmeans_p(Zs[:, cols], K_ANOM, KMEANS_ITER, rng, "index měsíce")
    lab = prirad(Zs[:, cols], Csub, 20000)[0]
    cnt = np.bincount(lab, minlength=K_ANOM)
    CA = np.zeros((K_ANOM, 42), np.float32)
    for j in range(42): CA[:, j] = np.bincount(lab, weights=Zs[:, j], minlength=K_ANOM) / np.maximum(cnt, 1)
    for k in np.nonzero(cnt == 0)[0]: CA[k, cols] = Csub[k]
    CA = f16(CA)
    np.savez(cesta("kody.npz"), C1=C1, C2=C2, CA=CA, Csub=Csub.astype(np.float32), cols=np.array(cols))
    with open(cesta("kody.bin"), "wb") as f:      # pro PHP: C1, C2, CA jako float16 little-endian
        for a_ in (C1, C2, CA): f.write(np.ascontiguousarray(a_.astype("<f2")).tobytes())
    stav["B_hotovo"] = True; uloz_stav(stav)
    p("  kódové knihy uloženy (kody.npz, kody.bin)")


def nacti_kody():
    z = np.load(cesta("kody.npz"))
    return types.SimpleNamespace(C1=z["C1"], C2=z["C2"], CA=z["CA"], Csub=z["Csub"], cols=list(z["cols"]))


# =============================== KÓDOVÁNÍ A DEKÓDOVÁNÍ JEDNOHO BLOKU ===============================
def tnn_kvant(Pt_t, Pp_t, krok, strop):
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

def sig_expand(sc, G, n):
    sq = sig_dek(sc.reshape(n, G.nper, 12))
    return np.concatenate([np.repeat(sq[:, i][:, None, :], G.ny[i], axis=1) for i in range(G.nper)], axis=1)[..., None]

def bexp(Brec, G):
    n = Brec.shape[0]
    return np.concatenate([np.tile(Brec[:, i][:, None], (1, G.ny[i], 1, 1)).reshape(n, G.ny[i] * 36, NF) for i in range(G.nper)], axis=1)

def rekonstrukce(G, KOD, p1, p2, idx, sc):
    """Společné pro kódování i dekódování: (Pp, Bf) = dekódované spojité hodnoty a základ ve fyzikálních jednotkách."""
    n = p1.shape[0]
    Brec = (KOD.C1[p1.reshape(-1)] + KOD.C2[p1.reshape(-1), p2.reshape(-1)]).reshape(n, G.nper, 36, NF)
    sqx = sig_expand(sc, G, n)
    Arec = (KOD.CA[idx.reshape(n, G.NY, 12)] * sqx).reshape(n, G.NY, 36, NF)
    Br = bexp(Brec, G).reshape(n, G.NY, 36, NF)
    return fyz(Br + Arec), fyz(Br)

def kod_zaklad(KOD, Vx):
    k1 = prirad(Vx, KOD.C1, 20000)[0]; rs = Vx - KOD.C1[k1]; k2 = np.zeros(len(Vx), np.int32)
    for k in np.unique(k1):
        m = k1 == k; k2[m] = prirad(rs[m], KOD.C2[k], 20000)[0]
    return k1.astype(np.uint8), k2.astype(np.uint8)

def kodovat(G, KOD, Xw, Sraw, jih):
    """Xw (n,NY,36,NF) pokřivená data; Sraw (n,NT) surové SDII float16; jih (n,) bool."""
    n = Xw.shape[0]
    Pa = fyz(Xw); Pa[..., 5:13] = np.clip(np.rint(Pa[..., 5:13]), 0, 11)
    d = {}
    d["sums"] = np.clip(np.rint(np.stack([Pa[..., IX[f]].sum(2) for f in COUNTS], axis=1)), 0, 511).astype(np.uint16)
    d["tgq"] = np.clip(np.rint((Pa[..., IX["tg"]].mean(2) + 75.0) / 0.5), 0, 254).astype(np.uint8)
    d["txxq"] = np.clip(np.rint((Pa[..., IX["txx"]].max(2) + 75.0) / 0.5), 0, 254).astype(np.uint8)
    d["ls"], d["fa"], _o = mraz_kody(Pa[..., IX["fd"]], jih)
    # základ
    Bs = np.empty((n, G.nper, 36, NF), np.float32); o2 = 0
    for i in range(G.nper):
        Bs[:, i] = Xw[:, o2:o2 + G.ny[i]].mean(axis=1); o2 += G.ny[i]
    k1, k2 = kod_zaklad(KOD, Bs.reshape(-1, NF))
    d["p1"] = k1.reshape(n, G.nper * 36); d["p2"] = k2.reshape(n, G.nper * 36)
    # anomálie a index
    A = np.empty_like(Xw); o2 = 0
    for i in range(G.nper):
        A[:, o2:o2 + G.ny[i]] = Xw[:, o2:o2 + G.ny[i]] - Bs[:, i][:, None]; o2 += G.ny[i]
    A = A.reshape(n, G.NY, 12, 42)
    Ac = A[..., KOD.cols]
    sig = np.empty((n, G.nper, 12), np.float32); o2 = 0
    for i in range(G.nper):
        sig[:, i] = np.sqrt((Ac[:, o2:o2 + G.ny[i]] ** 2).mean(axis=(1, 3))); o2 += G.ny[i]
    sc = sig_kod(sig)
    sqx = sig_expand(sc.reshape(n, -1), G, n)
    Z = A / sqx
    idx = prirad(Z.reshape(-1, 42)[:, KOD.cols], KOD.Csub, 20000)[0].reshape(n, G.NY * 12)
    d["idx"] = idx.astype(np.uint16); d["sc"] = sc.reshape(n, G.nper * 12)
    del A, Ac, Z
    # TNn po dekádách (zbytek proti predikci z indexu)
    Pp, _Bf = rekonstrukce(G, KOD, d["p1"], d["p2"], d["idx"], d["sc"])
    qt, qp, _qc = tnn_kvant(Pa[..., IX["tnn"]], Pp[..., IX["tnn"]], KROK, STROP)
    d["r"] = (qt - qp).reshape(n, G.NT).astype(np.int32)
    # výjimky SDII (skutečné extrémy, bez ořezu)
    ci, ti = np.nonzero(Sraw >= np.float16(SDII_PRAH))
    d["exc"] = np.stack([ci.astype(np.uint16), ti.astype(np.uint16), Sraw[ci, ti].view(np.uint16)], axis=0)
    return d

def vynut_db(P, Bf, S, ls, fa, tg_m, txx_m, jih, tnn_dek):
    """Jako vynut2 z odhad_vq3 (bez vypouštění), ale z uložených ročních hodnot místo skutečných dat."""
    n, ny_, _, _ = P.shape
    R = P.copy(); Bf = np.maximum(Bf, 0)
    ones = np.ones((n, ny_, 36), bool); nuly = np.zeros((n, ny_, 36)); c11 = np.full((n, ny_, 36), 11.0)
    g = lambda f: np.maximum(P[..., IX[f]], 0)
    Sf = lambda f: S[:, COUNTS.index(f), :].astype(np.float32)
    ordi = (np.arange(36)[None, :] + np.where(jih, 18, 0)[:, None]) % 36
    inv = np.argsort(ordi, axis=1)
    s = np.arange(36)[None, None, :]
    lsi = ls.astype(int)[..., None] - 1
    fai = np.where(fa > 0, fa.astype(int) + 17, 36)[..., None]
    nula = ((s <= 17) & (s > lsi)) | ((s >= 18) & (s < fai))
    hr = ((s == lsi) & (lsi >= 0)) | ((s == fai) & (fai < 36))
    pick = lambda m: np.take_along_axis(m, inv[:, None, :], axis=2)
    R[..., IX["fd"]] = skaluj(g("fd"), Sf("fd"), ~pick(nula), pick(hr).astype(np.float64), c11, Bf[..., IX["fd"]])
    R[..., IX["rr1"]] = skaluj(g("rr1"), Sf("rr1"), ones, nuly, c11, Bf[..., IX["rr1"]])
    R[..., IX["r20"]] = skaluj(g("r20"), Sf("r20"), ones, nuly, R[..., IX["rr1"]], Bf[..., IX["r20"]])
    R[..., IX["tr"]] = skaluj(g("tr"), Sf("tr"), ones, nuly, c11, Bf[..., IX["tr"]])
    cap = c11
    for f in ("su20", "su25", "su30", "su35"):
        R[..., IX[f]] = skaluj(g(f), Sf(f), ones, nuly, cap, Bf[..., IX[f]])
        cap = R[..., IX[f]]
    delta = (tg_m - P[..., IX["tg"]].mean(2))[..., None]
    for f in ("tg", "tn", "tx", "tnn", "txx"): R[..., IX[f]] = P[..., IX[f]] + delta
    tnn = np.asarray(tnn_dek, np.float32).copy()
    txx = np.minimum(R[..., IX["txx"]], txx_m[..., None])
    np.put_along_axis(txx, txx.argmax(2)[..., None], txx_m[..., None], axis=2)
    R[..., IX["tnn"]] = tnn; R[..., IX["txx"]] = txx
    tn = np.maximum(R[..., IX["tn"]], tnn); tx = np.minimum(R[..., IX["tx"]], txx)
    R[..., IX["tn"]] = tn; R[..., IX["tx"]] = tx
    R[..., IX["tg"]] = np.minimum(np.maximum(R[..., IX["tg"]], tn), tx)
    R[..., IX["sdii"]] = np.where(R[..., IX["rr1"]] <= 0, 0, np.maximum(P[..., IX["sdii"]], 1.0))
    return R

def dekodovat(G, KOD, d, jih):
    """d: slovník polí pro n buněk (jako z kodovat/načtení bloku) -> R (n,NY,36,NF) ve fyzikálních jednotkách."""
    n = d["p1"].shape[0]
    Pp, Bf = rekonstrukce(G, KOD, d["p1"], d["p2"], d["idx"], d["sc"])
    qt, qp, qc = tnn_kvant(Pp[..., IX["tnn"]], Pp[..., IX["tnn"]], KROK, STROP)   # qp z predikce (qt se přepíše)
    qt = qp + d["r"].reshape(n, G.NY, 36)
    tnn_dek = tnn_dekod(qt, Pp[..., IX["tnn"]], KROK, qc)
    tg_m = d["tgq"].astype(np.float32) * np.float32(0.5) - np.float32(75.0)
    txx_m = d["txxq"].astype(np.float32) * np.float32(0.5) - np.float32(75.0)
    R = vynut_db(Pp, Bf, d["sums"], d["ls"], d["fa"], tg_m, txx_m, jih, tnn_dek)
    ex = d["exc"]
    if ex.shape[1]:
        ci = ex[0].astype(np.int64); ti = ex[1].astype(np.int64)
        R[ci, ti // 36, ti % 36, IX["sdii"]] = ex[2].view(np.float16).astype(np.float32)
    return R


# =============================== ZABALENÍ DO STREAMŮ ===============================
def filtr(a, lg):
    if lg == 0: return a
    y = a.copy(); y[:, lg:] = a[:, lg:] - a[:, :-lg]; return y

def inv_filtr(y, lg):
    if lg == 0: return y
    n, m = y.shape
    mm = -(-m // lg) * lg
    pad = np.zeros((n, mm), np.uint8); pad[:, :m] = y
    return np.cumsum(pad.reshape(n, mm // lg, lg), axis=1, dtype=np.uint8).reshape(n, mm)[:, :m]

def zabal(a2, lagy):
    best = None
    for lg in lagy:
        z = zlib.compress(np.ascontiguousarray(filtr(a2, lg)).tobytes(), GZ)
        if best is None or len(z) < len(best[0]): best = (z, lg)
    return best

def rozbal(data, lg, n, m):
    a = np.frombuffer(zlib.decompress(data), np.uint8).reshape(n, m)
    return inv_filtr(a, lg)

def roviny(v, m):
    """9bitové hodnoty -> (spodní bajt, [bitové roviny nad 8. bitem zabalené po osmi])."""
    v = v.astype(np.int32)
    return (v & 255).astype(np.uint8), [np.packbits(((v >> b) & 1).astype(np.uint8), axis=1) for b in range(8, NB_IDX)]

def casti_bloku(G, d):
    """Vrací {sklad: [(jméno, bajty, lag), ...]} pro jeden blok buněk."""
    n = d["p1"].shape[0]
    S = {}
    c = []
    for nm in ("p1", "p2", "sc"):
        z, lg = zabal(d[nm], (0,)); c.append((nm, z, lg))
    S["zaklad"] = c
    c = []
    lo, hi = roviny(d["idx"], G.NY * 12)
    z, lg = zabal(lo, (0, 1, 12)); c.append(("idx_lo", z, lg))
    for k, pl in enumerate(hi):
        z, lg = zabal(pl, (0,)); c.append((f"idx_hi{k}", z, lg))
    S["index"] = c
    c = []
    for j, f in enumerate(COUNTS):
        lo, hi = roviny(d["sums"][:, j, :], G.NY)
        z, lg = zabal(lo, (0, 1)); c.append((f"sum_{f}_lo", z, lg))
        z, lg = zabal(np.packbits((d["sums"][:, j, :] >> 8).astype(np.uint8) & 1, axis=1), (0,)); c.append((f"sum_{f}_hi", z, lg))
    for nm in ("tgq", "txxq", "ls", "fa"):
        z, lg = zabal(d[nm], (0, 1)); c.append((nm, z, lg))
    S["rocni"] = c
    zz = ((d["r"] << 1) ^ (d["r"] >> 31)).astype(np.uint32)
    c = []
    z, lg = zabal((zz & 255).astype(np.uint8), (0,)); c.append(("r_lo", z, lg))
    z, lg = zabal(((zz >> 8) & 255).astype(np.uint8), (0,)); c.append(("r_hi", z, lg))
    S["tnn"] = c
    ex = d["exc"]
    raw = np.ascontiguousarray(np.stack([(ex[0] & 255), (ex[0] >> 8), (ex[1] & 255), (ex[1] >> 8), (ex[2] & 255), (ex[2] >> 8)]).astype(np.uint8))
    S["sdii"] = [("exc", zlib.compress(raw.tobytes(), GZ), 0)]
    return S

def nacti_blok(G, stav, blok, n):
    """Načte a dekóduje streamy jednoho bloku -> slovník polí (n buněk). Používá se při ověření."""
    tab = stav["sklady"]
    def cast(sklad, jmeno):
        for (nm, off, ln, lg) in tab[sklad]["tabulka"][blok]:
            if nm == jmeno:
                with open(cesta(sklad + ".dat"), "rb") as f:
                    f.seek(off); return f.read(ln), lg
        raise KeyError(jmeno)
    def a2(sklad, jmeno, m):
        dat, lg = cast(sklad, jmeno); return rozbal(dat, lg, n, m)
    def bit9(sklad, lo_name, hi_names, m):
        v = a2(sklad, lo_name, m).astype(np.int32)
        for k, hn in enumerate(hi_names):
            pl = rozbal(cast(sklad, hn)[0], 0, n, -(-m // 8))
            v |= (np.unpackbits(pl, axis=1)[:, :m].astype(np.int32) << (8 + k))
        return v
    d = {}
    d["p1"] = a2("zaklad", "p1", G.nper * 36); d["p2"] = a2("zaklad", "p2", G.nper * 36); d["sc"] = a2("zaklad", "sc", G.nper * 12)
    d["idx"] = bit9("index", "idx_lo", [f"idx_hi{k}" for k in range(NB_IDX - 8)], G.NY * 12).astype(np.uint16)
    sums = np.empty((n, 8, G.NY), np.uint16)
    for j, f in enumerate(COUNTS):
        lo = a2("rocni", f"sum_{f}_lo", G.NY).astype(np.int32)
        hi = np.unpackbits(rozbal(cast("rocni", f"sum_{f}_hi")[0], 0, n, -(-G.NY // 8)), axis=1)[:, :G.NY].astype(np.int32)
        sums[:, j, :] = lo | (hi << 8)
    d["sums"] = sums
    for nm in ("tgq", "txxq", "ls", "fa"): d[nm] = a2("rocni", nm, G.NY)
    zz = a2("tnn", "r_lo", G.NT).astype(np.uint32) | (a2("tnn", "r_hi", G.NT).astype(np.uint32) << 8)
    zz = zz.astype(np.int64); d["r"] = ((zz >> 1) ^ -(zz & 1)).astype(np.int32)
    dat, _lg = cast("sdii", "exc"); raw = np.frombuffer(zlib.decompress(dat), np.uint8)
    nrec = len(raw) // 6; raw = raw.reshape(6, nrec).astype(np.uint16)
    d["exc"] = np.stack([raw[0] | (raw[1] << 8), raw[2] | (raw[3] << 8), raw[4] | (raw[5] << 8)], axis=0)
    return d


# =============================== FÁZE C: kódování po dávkách ===============================
def faze_C(G, stav):
    if stav.get("C_hotovo"): return
    p(""); p(f"[C] Kódování buněk po dávkách (celkem {G.nblok} bloků po {CB} buňkách, dávka {DAVKA} bloků)")
    KOD = nacti_kody()
    NT, NC, NY = G.NT, G.NC, G.NY
    mmT, mmC, mmS = otevri_uloziste(G, "r")
    skl = stav.setdefault("sklady", {})
    for jm in ("zaklad", "index", "rocni", "tnn", "sdii"):
        s = skl.setdefault(jm, {"velikost": 0, "tabulka": []})
        s["tabulka"] = s["tabulka"][:stav.get("C_bloku", 0)]
        fp = cesta(jm + ".dat")
        if os.path.exists(fp):
            with open(fp, "r+b") as f: f.truncate(s["velikost"])
        else:
            open(fp, "wb").close()
    hotovo = stav.get("C_bloku", 0)
    ov = None
    if os.path.exists(cesta("overeni.npz")):
        z = np.load(cesta("overeni.npz")); ov = {"bunka": list(z["bunka"]), "data": z["data"]}
    else:
        ov = {"bunka": [], "data": np.zeros((0, NY, 36, NF), np.float16)}
    celkem_exc = stav.get("C_exc", 0)
    h0 = min(hotovo * CB, NC); t0 = time.time(); davek = 0
    for bt0 in range(hotovo, G.nblok, DAVKA):
        bt1 = min(bt0 + DAVKA, G.nblok); c0 = bt0 * CB; c1 = min(bt1 * CB, NC); ncb = c1 - c0
        # --- načtení dávky z úložiště
        Xs = np.zeros((ncb, NT, NF), np.float16); Sx = np.zeros((ncb, NT), np.float16)
        for g0 in range(0, NT, 72):
            g1 = min(g0 + 72, NT)
            Tm = np.asarray(mmT[g0:g1, c0:c1, :]); Cm = np.asarray(mmC[g0:g1, c0:c1, :]); Sm = np.asarray(mmS[g0:g1, c0:c1])
            F = F_z_uloziste(Tm, Cm, Sm)
            Xs[:, g0:g1, :] = F.transpose(1, 0, 2).astype(np.float16); Sx[:, g0:g1] = Sm.T
        for bl in range(bt0, bt1):
            b0 = bl * CB - c0; b1 = min((bl + 1) * CB, NC) - c0; nbl = b1 - b0
            ds_ = []
            for s0 in range(b0, b1, CH):
                s1 = min(s0 + CH, b1)
                Xw = Xs[s0:s1].astype(np.float32).reshape(s1 - s0, NY, 36, NF)
                jih = G.jih[c0 + s0:c0 + s1]
                d = kodovat(G, KOD, Xw, Sx[s0:s1], jih)
                d["exc"][0] = d["exc"][0] + np.uint16(s0 - b0)
                if bl % 3 == 0 and s0 == b0 and (c0 + s0) not in ov["bunka"] and len(ov["bunka"]) < 80:
                    Pa = fyz(Xw[:1]); Pa[..., 5:13] = np.clip(np.rint(Pa[..., 5:13]), 0, 11)
                    Pa[0, ..., IX["sdii"]] = Sx[s0].astype(np.float32).reshape(NY, 36)   # skutečné SDII bez ořezu
                    ov["bunka"].append(c0 + s0); ov["data"] = np.concatenate([ov["data"], Pa.astype(np.float16)], axis=0)
                ds_.append(d)
                prog(min(c0 + s1, NC), NC, h0, t0, "buňky")
            d = {k: (np.concatenate([x[k] for x in ds_], axis=1 if k == "exc" else 0)) for k in ds_[0]}
            celkem_exc += d["exc"].shape[1]
            for jm, casti in casti_bloku(G, d).items():
                s = skl[jm]; zaz = []
                with open(cesta(jm + ".dat"), "ab") as f:
                    for (nm, z, lg) in casti:
                        zaz.append([nm, s["velikost"], len(z), lg]); f.write(z); s["velikost"] += len(z)
                s["tabulka"].append(zaz)
            stav["C_bloku"] = bl + 1
        stav["C_exc"] = celkem_exc
        np.savez(cesta("overeni.npz"), bunka=np.array(ov["bunka"]), data=ov["data"])
        uloz_stav(stav)
        del Xs, Sx
        davek += 1
        if (cas_vyprsel() or (MAXDAVEK and davek >= MAXDAVEK)) and bt1 < G.nblok: raise Pauza()
    stav["C_hotovo"] = True; uloz_stav(stav)
    p(f"  výjimek SDII (≥ {SDII_PRAH:g} mm/den) uloženo: {celkem_exc}")


# =============================== FÁZE D: závěr ===============================
def faze_D(G, stav):
    meta = {"verze": 1, "varianta": VARIANTA, "popis": "bezExt: základ(p1,p2)+index měsíce(9 bitů)+sigma; roční přesné součty počtů, TG, TXx, termíny mrazů; TNn po dekádách; výjimky SDII",
            "FEAT": FEAT, "COUNTS": COUNTS, "NF": NF, "nper": G.nper, "ny": G.ny, "NY": G.NY, "NT": G.NT, "nt_soubor": G.nt_soubor,
            "NC": G.NC, "blok_bunek": CB, "nblok": G.nblok, "K_ANOM": K_ANOM, "K1": K1, "K2": K2, "krok_tnn": KROK, "strop_tnn": STROP,
            "sdii_prah": SDII_PRAH, "SIG_MIN": SIG_MIN, "sigma_na_oktavu": 8,
            "roční_teploty": {"krok": 0.5, "od": -75.0, "max_kod": 254},
            "warp": {k: {"hranice": [float(x) for x in h], "sklony": [float(x) for x in s]} for k, (h, s) in WARP.items()},
            "priorita": {k: float(v) for k, v in PRIORITA.items()},
            "lat": [float(x) for x in G.lat], "lon": [float(x) for x in G.lon],
            "sklady": {jm: [[c[0] for c in tab[0]] if tab else [], s["velikost"]] for jm, s in stav["sklady"].items() for tab in [s["tabulka"]]}}
    with open(cesta("meta.json"), "w", encoding="utf-8") as f: json.dump(meta, f)
    with open(cesta("tabulky.json"), "w", encoding="utf-8") as f: json.dump({jm: s["tabulka"] for jm, s in stav["sklady"].items()}, f)
    with open(cesta("bunky.bin"), "wb") as f:
        f.write(np.stack([G.rows, G.cols], axis=1).astype("<u2").tobytes())
    p(""); p("[D] Hotovo. Velikosti souborů (MB):")
    tot = 0
    for jm in ("zaklad", "index", "rocni", "tnn", "sdii"):
        mb = os.path.getsize(cesta(jm + ".dat")) / 1e6; tot += mb; p(f"    {jm + '.dat':<12}{mb:9.1f}")
    kb = os.path.getsize(cesta("kody.bin")) / 1e6; tot += kb
    p(f"    {'kody.bin':<12}{kb:9.1f}\n    {'CELKEM':<12}{tot:9.1f}")
    stav["D_hotovo"] = True; uloz_stav(stav)
    if FLAG("--smazat-temp"):
        for n_ in ("tmp_T.i16", "tmp_C.u8", "tmp_S.f16"):
            if os.path.exists(cesta(n_)): os.remove(cesta(n_))
        p("  dočasné úložiště smazáno")
    else:
        p("  dočasné úložiště (tmp_*.*, cca 4 GB) lze smazat: python stavba_db.py --smazat-temp")


# =============================== OVĚŘENÍ ===============================
def overit(G, stav):
    p(""); p("OVĚŘENÍ: dekódování kontrolních buněk ze souborů databáze a porovnání se zdrojem")
    KOD = nacti_kody()
    z = np.load(cesta("overeni.npz")); buny = list(z["bunka"]); T = z["data"].astype(np.float32)
    E = []; 
    for k, c in enumerate(buny):
        bl = c // CB; n = min((bl + 1) * CB, G.NC) - bl * CB; r = c - bl * CB
        d = nacti_blok(G, stav, bl, n)
        dr = {a: (b[:, :] if a == "exc" else b[r:r + 1]) for a, b in d.items()}
        ex = d["exc"]; m = ex[0] == r
        dr["exc"] = np.stack([ex[0][m] * 0, ex[1][m], ex[2][m]], axis=0)
        R = dekodovat(G, KOD, dr, G.jih[c:c + 1])
        E.append(R[0] - T[k])
    E = np.stack(E)                                      # (nbuněk, NY, 36, NF)
    p(f"  kontrolních buněk: {len(buny)}; dekád na buňku: {G.NT}")
    p(f"  {'indikátor':<11}{'RMSE':>8}{'MAE':>8}{'bias':>8}{'P95':>8}{'P99':>8}{'max':>9}")
    for j, f in enumerate(FEAT):
        e = E[..., j].ravel(); a = np.abs(e)
        p(f"  {POPIS[f]:<11}{np.sqrt((e * e).mean()):>8.2f}{a.mean():>8.2f}{e.mean():>+8.2f}{np.percentile(a, 95):>8.2f}{np.percentile(a, 99):>8.2f}{a.max():>9.2f}")
    p("  (kontrolní buňky jsou vybrané z celého světa; hodnoty mají být blízké výsledkům odhad_vq3.py, scénář S1)")


# =============================== HLAVNÍ PROGRAM ===============================
def main():
    if FLAG("--znovu") and os.path.isdir(VYSTUP):
        shutil.rmtree(VYSTUP)
    os.makedirs(VYSTUP, exist_ok=True)
    stav = nacti_stav()
    if FLAG("--znovu-kodovani"):
        for k in ("C_hotovo", "C_bloku", "C_exc", "sklady", "D_hotovo"): stav.pop(k, None)
        for n_ in ("overeni.npz", "zaklad.dat", "index.dat", "rocni.dat", "tnn.dat", "sdii.dat"):
            if os.path.exists(cesta(n_)): os.remove(cesta(n_))
        uloz_stav(stav)
    p("=" * 78); p(f"STAVBA DATABÁZE – varianta {VARIANTA}, TNn krok {KROK:g} °C, strop {'žádný' if STROP is None else '%g °C' % STROP}, "
                   f"SDII výjimky od {SDII_PRAH:g} mm/den, K = {K_ANOM}"); p("=" * 78)
    G = priprava()
    p(f"  buněk pevniny: {G.NC}; let: {G.NY} (dekád {G.NT}); bloků po {CB} buňkách: {G.nblok}")
    sig = {"CB": CB, "K": K_ANOM, "krok": KROK, "strop": STROP, "NC": G.NC, "NT": G.NT}
    if stav.get("sig") and stav["sig"] != sig and not FLAG("--znovu"):
        sys.exit("Chyba: nastavení se liší od rozpracované stavby (stav.json). Použijte stejné parametry, nebo --znovu.")
    stav["sig"] = sig; uloz_stav(stav)
    try:
        if not FLAG("--overit") or not stav.get("D_hotovo"):
            faze_A(G, stav)
            faze_B(G, stav)
            faze_C(G, stav)
            if not stav.get("D_hotovo") or FLAG("--smazat-temp"): faze_D(G, stav)
        if FLAG("--overit"):
            overit(G, stav)
    except Pauza:
        p("\n  Pozastaveno (časový limit / počet dávek). Spusťte stejný příkaz znovu, stavba pokračuje.")
        return
    except KeyboardInterrupt:
        p("\n  Přerušeno. Spusťte stejný příkaz znovu, stavba pokračuje od poslední dokončené dávky.")
        return
    p(f"  čas tohoto běhu {cas(time.time() - T0)}")


if __name__ == "__main__":
    main()
