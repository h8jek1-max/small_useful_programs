#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ODHAD VELIKOSTI v1 – agroklimatické indikátory C3S: dekádní vs. měsíční uložení

Nic nezapisuje (kromě textového výpisu odhad_velikosti.txt). Změří, kolik by po kódování
a gzipu zabraly jednotlivé proudy dat, a sečte je do tří variant:

    A  vše dekádně
    B  jádro dekádně (TG, TN, TX, FD, RR1), extrémy a počty (TNn, TXx, TR, R20mm, SU) měsíčně
    C  vše měsíčně

Vstup: soubory FD_*.nc, TR_*.nc, RR1_*.nc, R20mm_*.nc, SU_*.nc, SDII_*.nc, TG_*.nc, TN_*.nc,
       TX_*.nc, TNn_*.nc, TXx_*.nc (ve složce databaze/ nebo vedle skriptu; po třech obdobích).

Spuštění:
    python odhad_indikatoru.py                 výchozí vzorek 108 kroků (3 roky) z každého souboru
    python odhad_indikatoru.py --kroky 72      menší vzorek (méně RAM i času)
    python odhad_indikatoru.py --scenar rcp6p0 když ve složce leží víc scénářů/modelů

Kódování (shodné s návrhem):
    TG            1 B, krok 0,4 °C od -55 °C
    TN, TX, TNn, TXx   1 B, kladný rozdíl od už zakódovaného TG/TN/TX, krok 0,25 °C (bez hromadění chyby)
    FD, TR, RR1, R20mm, SU20..35  počty dní (dekádně 0..15, měsíčně 0..255), případně po dvou do bajtu
    SDII          1 B, logaritmická stupnice (krok 5 %) – v návrhu se místo toho počítá RR/RR1

Filtry (modulo 256, bezeztrátové): raw, delta po kroku, sezónní delta (36 dekád / 12 měsíců),
a každé z nich ještě s prostorovou predikcí ze souseda. Sloupec "PHP" = jen časové filtry
(rychlý bodový dotaz), sloupec "JS" = s prostorovou predikcí (nutné rozbalit celou dlaždici).
Odhad se škáluje z vzorku na všech dekád; počítejte s nepřesností ±10–15 %.
"""
import sys, os, re, glob, gzip, time, math, calendar
import numpy as np
from netCDF4 import Dataset, num2date

# =============================== NASTAVENÍ ===============================
TILE        = 30
LAT_MIN     = -60.0
BLOK        = 12          # kroků načtených najednou (musí být násobek 3)
KROKU_NA_SOUBOR = 108     # vzorek z každého souboru (násobek 36 = celé roky)
TG_KROK, TG_OD = 0.4, -55.0
D_KROK      = 0.25
KOD_MAX     = 250
POC_MAX_D   = 15          # strop počtu dní v dekádě (4 bity)
LAG_D, LAG_M = 36, 12     # sezónní posun v dekádách / měsících
VZOREK_DLAZDIC = 8
GZ_VYBER, GZ_FINAL = 6, 9

PREFIXY = {"fd": "FD", "tr": "TR", "rr1": "RR1", "r20": "R20mm", "su": "SU", "sdii": "SDII",
           "tg": "TG", "tn": "TN", "tx": "TX", "tnn": "TNn", "txx": "TXx"}
TEPLOTY  = ["tg", "dtn", "dtx", "dtnn", "dtxx"]
POCTY    = ["fd", "tr", "rr1", "r20", "su20", "su25", "su30", "su35"]
JMENA    = TEPLOTY + POCTY + ["sdii"]
PARY     = [("fd", "rr1"), ("tr", "r20"), ("su20", "su25"), ("su30", "su35")]

def arg(nazev, vychozi):
    if nazev in sys.argv:
        i = sys.argv.index(nazev)
        if i + 1 < len(sys.argv): return sys.argv[i + 1]
    return vychozi

KROKU = int(arg("--kroky", KROKU_NA_SOUBOR))
SCENAR = arg("--scenar", None)
LOG = []
def p(s=""):
    print(s); LOG.append(s)


# =============================== POMOCNÉ ===============================
def cas(s):
    s = int(s); return f"{s // 60:02d}:{s % 60:02d}"

def progres(hotovo, celkem, t0, popis):
    pr = hotovo / celkem if celkem else 1
    el = time.time() - t0
    zbyva = el / pr - el if pr > 0 else 0
    sys.stdout.write(f"\r  {popis}: {hotovo}/{celkem} ({pr * 100:5.1f} %)  uplynulo {cas(el)}  zbývá ~{cas(zbyva)}   ")
    sys.stdout.flush()
    if hotovo >= celkem: print()

def najdi_soubory(scenar):
    vse = sorted(set(glob.glob("databaze/*.nc") + glob.glob("*.nc")))
    klice = {}   # indikátor -> {klíč scénáře -> [(perioda, cesta)]}
    for k, pref in PREFIXY.items():
        klice[k] = {}
        for pth in vse:
            b = os.path.basename(pth)
            if not b.startswith(pref + "_"): continue
            m = re.search(r"(\d{8})-(\d{8})", b)
            kl = re.sub(r"\d{8}-\d{8}", "", b[len(pref):])
            klice[k].setdefault(kl, []).append((m.group(1) if m else b, pth))
        if not klice[k]:
            sys.exit(f"Chyba: nenalezeny soubory {pref}_*.nc (hledáno v databaze/ a v aktuální složce)")
    spolecne = set.intersection(*[set(v) for v in klice.values()])
    if scenar: spolecne = {x for x in spolecne if scenar in x}
    if not spolecne: sys.exit("Chyba: indikátory nemají společný model/scénář (zkuste --scenar rcp6p0)")
    if len(spolecne) > 1:
        pref6 = [x for x in spolecne if "rcp6p0" in x]
        vyber = sorted(pref6 or spolecne)[0]
        print(f"  upozornění: nalezeno víc scénářů/modelů {sorted(spolecne)}, používám '{vyber}'")
    else:
        vyber = next(iter(spolecne))
    return {k: [pth for _, pth in sorted(klice[k][vyber])] for k in PREFIXY}

def osy(var):
    d = [x.lower() for x in var.dimensions]
    def hledej(pref):
        for i, x in enumerate(d):
            if x.startswith(pref): return i
        return None
    o = {"t": hledej("time"), "lat": hledej("lat"), "lon": hledej("lon")}
    if None in o.values(): o = {"t": 0, "lat": var.ndim - 2, "lon": var.ndim - 1}
    zbytek = [i for i in range(var.ndim) if i not in o.values()]
    if zbytek: o["l"] = zbytek[0]
    return o

def najdi_promennou(nc, jmeno):
    if jmeno in nc.variables: return nc.variables[jmeno]
    for n, v in nc.variables.items():
        if n.lower() == jmeno.lower(): return v
    for n, v in nc.variables.items():
        if v.ndim >= 3 and not n.lower().endswith(("bnds", "bounds")) and "lat" in " ".join(v.dimensions).lower():
            print(f"  (proměnná '{jmeno}' nenalezena, používám '{n}')"); return v
    sys.exit(f"V souboru {nc.filepath()} není proměnná {jmeno}. Nalezeno: {list(nc.variables)}")

def cti_blok(var, o, t0, t1, flip, radku, nezaporne):
    idx = [slice(None)] * var.ndim
    idx[o["t"]] = slice(t0, t1)
    a = var[tuple(idx)]
    if np.ma.isMaskedArray(a): a = a.filled(np.nan)
    a = np.asarray(a, dtype=np.float32)
    with np.errstate(invalid="ignore"):
        a[a > 1e19] = np.nan
        if nezaporne:
            a[a < -1e-3] = np.nan
            np.maximum(a, 0, out=a)
    poradi = [o["t"]] + ([o["l"]] if "l" in o else []) + [o["lat"], o["lon"]]
    a = np.transpose(a, poradi)
    if flip: a = a[..., ::-1, :]
    return a[..., :radku, :]

def dny_v_krocich(nc, nt):
    t = nc.variables["time"]
    b = getattr(t, "bounds", None)
    if b and b in nc.variables:
        bb = np.asarray(nc.variables[b][:], dtype=np.float64)
        n = np.rint(bb[:, 1] - bb[:, 0]).astype(int)
        if n.min() >= 7 and n.max() <= 12: return n
    kal = getattr(t, "calendar", "standard").lower()
    dts = num2date(t[:], t.units, kal)
    n = np.zeros(nt, dtype=int)
    for i, d in enumerate(dts):
        if kal == "360_day": dim = 30
        elif kal in ("noleap", "365_day"): dim = [31,28,31,30,31,30,31,31,30,31,30,31][d.month - 1]
        else: dim = calendar.monthrange(d.year, d.month)[1]
        n[i] = 10 if d.day <= 20 else dim - 20
    return n

def je_kelvin(var, blok):
    u = getattr(var, "units", None)
    if u: return u.strip().lower() in ("k", "kelvin", "degk", "deg_k")
    return float(np.nanmean(blok)) > 150

def poradi_su(d, v, o, nt, flip, radku):
    """Pořadí hladin SU vzestupně podle prahu (20,25,30,35)."""
    dn = v.dimensions[o["l"]]
    if dn in d.variables:
        h = np.asarray(d.variables[dn][:], dtype=float)
        if h.size == 4 and np.all(np.isfinite(h)): return np.argsort(h)
    a = cti_blok(v, o, min(18, nt - 1), min(24, nt), flip, radku, True)
    prum = np.nanmean(a, axis=(0, 2, 3))
    print("  upozornění: osa prahů SU nemá hodnoty, pořadí určeno z dat (červenec)")
    return np.argsort(-prum)


# =============================== KÓDOVÁNÍ ===============================
orez, celkem_k = {}, {}
def kvant(x, krok, nazev, od=0.0, maxk=KOD_MAX):
    k = np.rint((x - od) / krok)
    zle = ~np.isfinite(k)
    k[zle] = 0
    mimo = (k < 0) | (k > maxk)
    orez[nazev] = orez.get(nazev, 0) + int(mimo.sum()) + int(zle.sum())
    celkem_k[nazev] = celkem_k.get(nazev, 0) + k.size
    return np.clip(k, 0, maxk).astype(np.uint8)

def rek(k, krok, od=0.0):
    return od + k.astype(np.float32) * np.float32(krok)

def kod_teplot(tg, tn, tx, tnn, txx, s):
    tgk = kvant(tg, TG_KROK, "tg" + s, TG_OD)
    tgr = rek(tgk, TG_KROK, TG_OD)
    dtn = kvant(tgr - tn, D_KROK, "dtn" + s)
    tnr = tgr - rek(dtn, D_KROK)
    dtx = kvant(tx - tgr, D_KROK, "dtx" + s)
    txr = tgr + rek(dtx, D_KROK)
    return {"tg": tgk, "dtn": dtn, "dtx": dtx,
            "dtnn": kvant(tnr - tnn, D_KROK, "dtnn" + s),
            "dtxx": kvant(txx - txr, D_KROK, "dtxx" + s)}

def kod_poctu(x, nazev, maxk):
    return kvant(x, 1.0, nazev, 0.0, maxk)

def sestav_sdii_lut():
    v = [0.0, 1.0]
    while v[-1] < 150: v.append(v[-1] * 1.05)
    return np.array(v)
SDII_LUT = sestav_sdii_lut()

def kod_sdii(x):
    return np.clip(np.searchsorted(SDII_LUT, np.nan_to_num(x)), 0, 255).astype(np.uint8)

def mes_vazeny(x, n3):
    nm = n3.shape[0]
    r = x.reshape(nm, 3, -1)
    return (r * n3[:, :, None]).sum(1) / n3.sum(1)[:, None]

def mes_soucet(x):
    return x.reshape(x.shape[0] // 3, 3, -1).sum(1)


# =============================== FILTRY A GZIP ===============================
def gz(b, uroven):
    return gzip.compress(b, compresslevel=uroven, mtime=0)

def predikce_dlazdice(sub):
    rr_, cc_ = np.nonzero(sub)
    rank = np.full(sub.shape, -1, np.int32); rank[rr_, cc_] = np.arange(len(rr_), dtype=np.int32)
    left = np.full(len(rr_), -1, np.int32); m = cc_ > 0; left[m] = rank[rr_[m], cc_[m] - 1]
    up = np.full(len(rr_), -1, np.int32);   m = rr_ > 0; up[m] = rank[rr_[m] - 1, cc_[m]]
    return np.where(left >= 0, left, up)

def filtr(x, P, lag, prostor):
    if lag and x.shape[1] > lag:
        d = x.copy(); d[:, lag:] = x[:, lag:] - x[:, :-lag]; x = d
    if prostor:
        r = x.copy(); m = P >= 0; r[m] = x[m] - x[P[m]]; x = r
    return x

def mody_pro(L):
    return [("raw", None, False), ("delta1", 1, False), (f"sezon{L}", L, False),
            ("prostor", None, True), ("prostor+delta1", 1, True), (f"prostor+sezon{L}", L, True)]


# =============================== HLAVNÍ ===============================
def main():
    t_start = time.time()
    if KROKU % 3 or KROKU < 12: sys.exit("--kroky musí být násobek 3 (doporučeno 36, 72 nebo 108)")
    soubory = najdi_soubory(SCENAR)
    nper = len(soubory["fd"])
    for k, v in soubory.items():
        if len(v) != nper: sys.exit(f"Chyba: {PREFIXY[k]} má {len(v)} souborů, FD {nper}")
    p("=" * 72); p("ODHAD VELIKOSTI – agroklimatické indikátory (dekádně vs. měsíčně)"); p("=" * 72)
    for i in range(nper):
        p(f"  období {i + 1}: " + os.path.basename(soubory["fd"][i]).split("_dek_")[-1])

    def otevri(i): return {k: Dataset(soubory[k][i]) for k in PREFIXY}
    def zavri(ds):
        for d in ds.values(): d.close()

    # ---- geometrie z TG ----
    ds = otevri(0)
    v0 = najdi_promennou(ds["tg"], "TG"); o0 = osy(v0)
    lat = np.asarray(ds["tg"].variables[v0.dimensions[o0["lat"]]][:], dtype=np.float64)
    lon = np.asarray(ds["tg"].variables[v0.dimensions[o0["lon"]]][:], dtype=np.float64)
    flip = lat[0] < lat[-1]
    lat_nu = lat[::-1] if flip else lat
    NCOLS = len(lon); radku = int(np.sum(lat_nu >= LAT_MIN - 1e-6))
    PY, PX = math.ceil(radku / TILE), math.ceil(NCOLS / TILE)
    zavri(ds)

    # ---- délky, kontrola časových os ----
    nt_soubor, nd_all, lim = [], [], []
    for i in range(nper):
        ds = otevri(i)
        nt = len(ds["fd"].dimensions["time"]); tref = ds["fd"].variables["time"][:]
        for k, d in ds.items():
            if len(d.dimensions["time"]) != nt: sys.exit(f"Chyba: různý počet kroků v souborech období {i + 1} ({PREFIXY[k]})")
            if not np.allclose(d.variables["time"][:], tref): sys.exit(f"Chyba: časová osa {PREFIXY[k]} se liší od FD (období {i + 1})")
        t = ds["fd"].variables["time"]
        d1 = num2date(t[:1], t.units, getattr(t, "calendar", "standard"))[0]
        if d1.month != 1 or d1.day > 10:
            print(f"  upozornění: soubor {i + 1} nezačíná první dekádou ledna – měsíční slučování po třech krocích je pak nepřesné")
        nd_all.append(dny_v_krocich(ds["fd"], nt)); nt_soubor.append(nt)
        lim.append((min(nt, KROKU) // 3) * 3)
        zavri(ds)
    NT_PLNE = sum(nt_soubor); NT_D = sum(lim); NT_M = NT_D // 3
    nd_vse = np.concatenate([nd_all[i][:lim[i]] for i in range(nper)])

    # ---- maska pevniny z TG ----
    ds = otevri(0)
    vtg = najdi_promennou(ds["tg"], "TG")
    tg0 = cti_blok(vtg, osy(vtg), 0, 1, flip, radku, False)[0]
    zavri(ds)
    maska = np.isfinite(tg0)
    rows_l, cols_l, t_start_idx, t_n, sub_masky, pos = [], [], [], [], [], 0
    for ty in range(PY):
        for tx in range(PX):
            sub = np.zeros((TILE, TILE), bool)
            blok = maska[ty * TILE:(ty + 1) * TILE, tx * TILE:(tx + 1) * TILE]
            sub[:blok.shape[0], :blok.shape[1]] = blok
            r_, c_ = np.nonzero(sub)
            rows_l.append(ty * TILE + r_); cols_l.append(tx * TILE + c_)
            t_start_idx.append(pos); t_n.append(len(r_)); pos += len(r_); sub_masky.append(sub)
    rows_all = np.concatenate(rows_l); cols_all = np.concatenate(cols_l)
    NC = len(rows_all)
    p(f"  mřížka {NCOLS} x {len(lat)}, ukládá se {radku} řad (šířka >= {LAT_MIN}), dlaždice {TILE}x{TILE} -> {PX} x {PY}")
    p(f"  buněk pevniny: {NC}, neprázdných dlaždic: {sum(1 for n in t_n if n)}")
    p(f"  úplná data: {NT_PLNE} dekád ({NT_PLNE // 3} měsíců); vzorek: {NT_D} dekád ({NT_M} měsíců) = {NT_D / NT_PLNE * 100:.1f} %")
    ram = NC * (NT_D + NT_M) * len(JMENA) / 1e6
    p(f"  RAM pro kódy: ~{ram:.0f} MB")

    DEK = {n: np.zeros((NC, NT_D), np.uint8) for n in JMENA}
    MES = {n: np.zeros((NC, NT_M), np.uint8) for n in JMENA}

    # ---- čtení a kódování ----
    p(""); p("[1/3] Čtení .nc souborů a kódování")
    print()
    t0 = time.time(); g = 0
    kel = {}
    for i in range(nper):
        ds = otevri(i)
        var = {k: najdi_promennou(ds[k], PREFIXY[k]) for k in PREFIXY}
        o = {k: osy(var[k]) for k in PREFIXY}
        if "l" not in o["su"]: sys.exit("SU nemá osu prahů")
        por_su = poradi_su(ds["su"], var["su"], o["su"], nt_soubor[i], flip, radku)
        print(f"  soubor {i + 1}/{nper}: pořadí prahů SU (index -> °C): " +
              ", ".join(f"{int(ix)}->{20 + 5 * j}" for j, ix in enumerate(por_su)))
        for a in range(0, lim[i], BLOK):
            b = min(a + BLOK, lim[i]); nb = b - a; nm = nb // 3
            def rd(k, nez):
                return cti_blok(var[k], o[k], a, b, flip, radku, nez)[..., rows_all, cols_all]
            T = {}
            for k in ("tg", "tn", "tx", "tnn", "txx"):
                x = rd(k, False)
                if (i, k) not in kel: kel[(i, k)] = je_kelvin(var[k], x)
                if kel[(i, k)]: x = x - np.float32(273.15)
                T[k] = x
            C = {k: rd(k, True) for k in ("fd", "tr", "rr1", "r20")}
            su = rd("su", True)[:, por_su, :]
            sd = rd("sdii", True)
            n = nd_vse[g:g + nb].astype(np.float32); n3 = n.reshape(nm, 3)
            gm = g // 3

            kd = kod_teplot(T["tg"], T["tn"], T["tx"], T["tnn"], T["txx"], "/dek")
            for k, v in kd.items(): DEK[k][:, g:g + nb] = v.T
            km = kod_teplot(mes_vazeny(T["tg"], n3), mes_vazeny(T["tn"], n3), mes_vazeny(T["tx"], n3),
                            T["tnn"].reshape(nm, 3, -1).min(1), T["txx"].reshape(nm, 3, -1).max(1), "/mes")
            for k, v in km.items(): MES[k][:, gm:gm + nm] = v.T

            for k in ("fd", "tr", "rr1", "r20"):
                DEK[k][:, g:g + nb] = kod_poctu(C[k], k + "/dek", POC_MAX_D).T
                MES[k][:, gm:gm + nm] = kod_poctu(mes_soucet(C[k]), k + "/mes", 255).T
            for j, nm_ in enumerate(("su20", "su25", "su30", "su35")):
                DEK[nm_][:, g:g + nb] = kod_poctu(su[:, j, :], nm_ + "/dek", POC_MAX_D).T
                MES[nm_][:, gm:gm + nm] = kod_poctu(mes_soucet(su[:, j, :]), nm_ + "/mes", 255).T
            DEK["sdii"][:, g:g + nb] = kod_sdii(sd).T
            vazy = np.nan_to_num(C["rr1"]); cit = vazy.reshape(nm, 3, -1).sum(1)
            cit_s = (np.nan_to_num(sd) * vazy).reshape(nm, 3, -1).sum(1)
            MES["sdii"][:, gm:gm + nm] = kod_sdii(np.where(cit > 0, cit_s / np.maximum(cit, 1e-9), 0)).T
            g += nb
            progres(g, NT_D, t0, "dekády")
        zavri(ds)

    p("")
    zle = {k: v / celkem_k[k] * 100 for k, v in orez.items() if v / celkem_k[k] > 0.01}
    if zle:
        p("  upozornění: podíl hodnot oříznutých/chybějících při kódování (%): " +
          ", ".join(f"{k} {v:.2f}" for k, v in sorted(zle.items())))
    else:
        p("  kódování bez oříznutí (všechny proudy < 0,01 % mimo rozsah)")

    # ---- filtry a gzip ----
    p(""); p("[2/3] Výběr filtru (vzorek dlaždic) a gzip všech dlaždic")
    print()
    predikce = [predikce_dlazdice(s) for s in sub_masky]
    nepr = [i for i in range(len(t_n)) if t_n[i] > 0]
    vz = [nepr[j] for j in np.linspace(0, len(nepr) - 1, min(VZOREK_DLAZDIC, len(nepr))).astype(int)]
    kus = lambda pole, i: pole[t_start_idx[i]:t_start_idx[i] + t_n[i]]
    vysl = {}   # (druh, jméno) -> dict

    def zpracuj(druh, jmeno, pole, lag, scale, jen_raw=False):
        mody = [("raw", None, False)] if jen_raw else mody_pro(lag)
        vel = {}
        for nazev, lg, pr in mody:
            vel[nazev] = sum(len(gz(np.ascontiguousarray(filtr(kus(pole, i), predikce[i], lg, pr)).tobytes(), GZ_VYBER)) for i in vz)
        php = {m: v for m, v in vel.items() if "prostor" not in m}
        f_php = min(php, key=php.get); f_any = min(vel, key=vel.get)
        pom = vel[f_any] / vel[f_php] if vel[f_php] else 1.0
        spec = {m[0]: m for m in mody}[f_php]
        celk = 0
        for i in nepr:
            celk += len(gz(np.ascontiguousarray(filtr(kus(pole, i), predikce[i], spec[1], spec[2])).tobytes(), GZ_FINAL))
        vysl[(druh, jmeno)] = {"php": celk * scale / 1e6, "js": celk * scale * pom / 1e6, "f_php": f_php, "f_js": f_any}

    sc_d = NT_PLNE / NT_D; sc_m = (NT_PLNE // 3) / NT_M
    celkem_proudu = len(JMENA) * 2 + len(PARY)
    t0 = time.time(); hotovo = 0
    for jm in JMENA:
        zpracuj("dek", jm, DEK[jm], LAG_D, sc_d); hotovo += 1; progres(hotovo, celkem_proudu, t0, "proudy")
        zpracuj("mes", jm, MES[jm], LAG_M, sc_m); hotovo += 1; progres(hotovo, celkem_proudu, t0, "proudy")
    for a_, b_ in PARY:
        zpracuj("par", a_ + "+" + b_, (DEK[a_] << 4) | DEK[b_], LAG_D, sc_d, jen_raw=True)
        hotovo += 1; progres(hotovo, celkem_proudu, t0, "proudy")

    # ---- souhrn ----
    p(""); p("[3/3] Výsledky (odhad pro všech %d dekád, MB po gzipu)" % NT_PLNE); p("")
    p(f"{'proud':<8}{'dekádně PHP':>13}{'filtr':>16}{'dekádně JS':>12}{'měsíčně PHP':>13}{'filtr':>16}{'měsíčně JS':>12}")
    for jm in JMENA:
        d, m = vysl[("dek", jm)], vysl[("mes", jm)]
        p(f"{jm:<8}{d['php']:>13.1f}{d['f_php']:>16}{d['js']:>12.1f}{m['php']:>13.1f}{m['f_php']:>16}{m['js']:>12.1f}")
    p("")
    p("Dekádní počty po dvou do bajtu (nibble) vs. zvlášť:")
    for a_, b_ in PARY:
        pr = vysl[("par", a_ + "+" + b_)]
        zv = vysl[("dek", a_)]["php"] + vysl[("dek", b_)]["php"]
        p(f"  {a_}+{b_:<6} spolu {pr['php']:6.1f} MB   zvlášť {zv:6.1f} MB   -> " + ("spojit" if pr['php'] < zv else "nechat zvlášť"))

    def par_dek(a_, b_, kl):
        return min(vysl[("par", a_ + "+" + b_)][kl], vysl[("dek", a_)][kl] + vysl[("dek", b_)][kl])
    def soucet(druh, jm, kl): return sum(vysl[(druh, x)][kl] for x in jm)

    varianty = {}
    for kl in ("php", "js"):
        pocty_d = sum(par_dek(a_, b_, kl) for a_, b_ in PARY)
        A = soucet("dek", TEPLOTY, kl) + pocty_d
        B = (soucet("dek", ["tg", "dtn", "dtx"], kl) + par_dek("fd", "rr1", kl) +
             soucet("mes", ["dtnn", "dtxx", "tr", "r20", "su20", "su25", "su30", "su35"], kl))
        C = soucet("mes", TEPLOTY + POCTY, kl)
        sd = {"A": vysl[("dek", "sdii")][kl], "B": vysl[("dek", "sdii")][kl], "C": vysl[("mes", "sdii")][kl]}
        varianty[kl] = {"A": (A, sd["A"]), "B": (B, sd["B"]), "C": (C, sd["C"])}

    raw = {"A": NC * NT_PLNE * (5 + 4) / 1e6,
           "B": (NC * NT_PLNE * 4 + NC * (NT_PLNE // 3) * 8) / 1e6,
           "C": NC * (NT_PLNE // 3) * 13 / 1e6}
    popis = {"A": "vše dekádně", "B": "jádro dekádně, extrémy/počty měsíčně", "C": "vše měsíčně"}
    p(""); p("=" * 72)
    p(f"{'varianta':<42}{'surově':>8}{'PHP (časové filtry)':>22}{'JS (+prostor)':>16}")
    for v in "ABC":
        a1, s1 = varianty["php"][v]; a2, s2 = varianty["js"][v]
        p(f"{v}  {popis[v]:<39}{raw[v]:>7.0f}M{a1:>14.0f} MB (+{s1:.0f} SDII){a2:>8.0f} MB (+{s2:.0f})")
    p("=" * 72)
    p(f"  surově = bajty před gzipem; v A/B/C bez SDII (v návrhu se počítá RR/RR1); SDII je uvedeno zvlášť.")
    p(f"  k tomu: rank.bin ~ {NCOLS * radku * 2 / 1e3:.0f} kB, rejstřík dlaždic jednotky kB, TG/H1 vztah netestován.")
    p(f"  klimatologie pro mapu (3 období x 36 dekád, 13 proudů): ~{NC * 108 * 13 / 1e6:.0f} MB před gzipem.")
    p(f"  vzorek tvoří prvních {KROKU} kroků z každého souboru; odhad ±10–15 %, trendy v dalších letech jsou mimo vzorek.")
    p(f"  celkový čas {cas(time.time() - t_start)}")
    with open("odhad_velikosti.txt", "w", encoding="utf-8") as f: f.write("\n".join(LOG) + "\n")
    print("  (výpis uložen do odhad_velikosti.txt)")


if __name__ == "__main__":
    main()
