#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SLUČOVAČ v2 – srážky (RR) + BEDD (6 hladin) z C3S .nc souborů do jedné gzip databáze dlaždic.

Spuštění:
    python slucovac.py --odhad     rychlý test: zpracuje jen ~2 roky z každého souboru a odhadne
                                   výslednou velikost (nic nezapisuje)
    python slucovac.py             plný běh

Co se ukládá (na jednu buňku pevniny a jednu dekádu = 2 bajty):
    RR  – 1 bajt: srážky na stupnici s krokem max(0,25 mm; 5 % hodnoty), 255 = chybí
    H1  – 1 bajt: průměrná denní "teplota nad 0 °C" = BEDD(0 °C)/počet dní v dekádě, krok 0,16 °C/den
Hladiny 2/4/6/8/10 °C dopočítá JavaScript z H1 pomocí tabulky střední hodnoty (kalibrace v JSONu).

Paměť: ~ 2 x 220 MB pole kódů + bloky po BLOK krocích. Dočasné soubory na disku nejsou potřeba.
"""
import sys, os, re, glob, json, gzip, time, base64, calendar, math
import numpy as np
from netCDF4 import Dataset, num2date

# =============================== NASTAVENÍ ===============================
VYSTUP_DATA     = "celysvet_agric_c3s.gz.bin"
VYSTUP_REJSTRIK = "registrace_agric_c3s.json"

TILE        = 30        # velikost dlaždice v buňkách (720/30 = 24, 300/30 = 10 řad)
LAT_MIN     = -60.0     # pod touto šířkou (Antarktida) se nic neukládá
BLOK        = 12        # kolik časových kroků se čte z disku najednou (víc = rychlejší HDD, víc RAM)

# srážky: stupnice kódů 0..K (0 = sucho). Krok = max(RR_ABS, RR_REL * hodnota)
RR_ABS      = 0.25      # mm
RR_REL      = 0.05      # 5 % -> max. chyba ±2,5 %; 0.08 = menší soubor, hrubší srážky
RR_MAX      = 2000.0    # mm / dekádu (strop stupnice)

# teplota: H1 = BEDD(0 °C) / počet dní v dekádě  [°C/den]
H1_KROK     = 0.16
H1_MAX_KOD  = 250       # 250 * 0.16 = 40 °C/den
KOD_CHYBI   = 255

PASMA_LAT   = [23.5, 40.0, 55.0]   # hranice |šířky| pro samostatné kalibrační tabulky (4 pásma)
MIN_VZORKU  = 200                  # řádek kalibrace s méně vzorky se nahradí průměrem přes pásma
GZIP_UROVEN = 9
FILTR_RR    = None      # None = vybere se automaticky; 0 raw, 1 časová delta, 2 prostorová, 3 obojí
FILTR_H1    = None
ODHAD_KROKU = 72        # kolik kroků z každého souboru se zpracuje v režimu --odhad

ODHAD = "--odhad" in sys.argv
PRAHY = [0, 2, 4, 6, 8, 10]


# =============================== POMOCNÉ ===============================
def cas(s):
    s = int(s); return f"{s // 60:02d}:{s % 60:02d}"

def progres(hotovo, celkem, t0, popis):
    pr = hotovo / celkem
    el = time.time() - t0
    zbyva = el / pr - el if pr > 0 else 0
    sys.stdout.write(f"\r  {popis}: {hotovo}/{celkem} ({pr * 100:5.1f} %)  uplynulo {cas(el)}  zbývá ~{cas(zbyva)}   ")
    sys.stdout.flush()
    if hotovo >= celkem: print()

def najdi_soubory():
    rr = sorted(glob.glob("databaze/RR_*.nc") or glob.glob("RR_*.nc"))
    bd = sorted(glob.glob("databaze/BEDD_*.nc") or glob.glob("BEDD_*.nc"))
    if not rr or len(rr) != len(bd):
        sys.exit(f"Chyba: počet souborů nesouhlasí nebo žádné nejsou (RR: {len(rr)}, BEDD: {len(bd)})")
    def perioda(p):
        m = re.search(r"(\d{8})-(\d{8})", os.path.basename(p)); return m.group(0) if m else None
    pr, pb = [perioda(p) for p in rr], [perioda(p) for p in bd]
    if None not in pr and None not in pb:
        if sorted(pr) != sorted(pb): sys.exit(f"Chyba: období RR {pr} a BEDD {pb} si neodpovídají")
        rr = [p for _, p in sorted(zip(pr, rr))]
        bd = sorted(bd, key=perioda)
    return list(zip(rr, bd))

def osy(var):
    """Najde osy proměnné podle jmen; vrací dict t, lat, lon, (l = hladiny)."""
    d = [x.lower() for x in var.dimensions]
    def hledej(pref):
        for i, x in enumerate(d):
            if x.startswith(pref): return i
        return None
    o = {"t": hledej("time"), "lat": hledej("lat"), "lon": hledej("lon")}
    if None in o.values():                       # záloha: (time, [level], lat, lon)
        o = {"t": 0, "lat": var.ndim - 2, "lon": var.ndim - 1}
    zbytek = [i for i in range(var.ndim) if i not in o.values()]
    if zbytek: o["l"] = zbytek[0]
    return o

def najdi_promennou(nc, jmeno):
    if jmeno in nc.variables: return nc.variables[jmeno]
    for n, v in nc.variables.items():
        if v.ndim >= 3 and not n.lower().endswith(("bnds", "bounds")) and "lat" in " ".join(v.dimensions).lower():
            print(f"  (proměnná '{jmeno}' nenalezena, používám '{n}')"); return v
    sys.exit(f"V souboru {nc.filepath()} není proměnná {jmeno}. Nalezeno: {list(nc.variables)}")

def cti_blok(var, o, t0, t1, flip, radku):
    """Vrátí pole (t, [hladiny], lat severem nahoru oříznuté na radku, lon) ve float32 s NaN."""
    idx = [slice(None)] * var.ndim
    idx[o["t"]] = slice(t0, t1)
    a = var[tuple(idx)]
    if np.ma.isMaskedArray(a): a = a.filled(np.nan)
    a = np.asarray(a, dtype=np.float32)
    a[a > 1e19] = np.nan
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
        if kal in ("360_day",): dim = 30
        elif kal in ("noleap", "365_day"): dim = [31,28,31,30,31,30,31,31,30,31,30,31][d.month - 1]
        else: dim = calendar.monthrange(d.year, d.month)[1]
        n[i] = 10 if d.day <= 20 else dim - 20
    return n

def rok_mesic_start(nc):
    t = nc.variables["time"]
    d = num2date(t[:1], t.units, getattr(t, "calendar", "standard"))[0]
    return d.year


# =============================== STUPNICE SRÁŽEK ===============================
def sestav_rr_lut():
    v = [0.0]
    while v[-1] < RR_MAX:
        v.append(v[-1] + max(RR_ABS, v[-1] * RR_REL))
    lut = np.array(v, dtype=np.float64)
    if len(lut) > 250: sys.exit(f"RR stupnice má {len(lut)} kódů (max 250) – zvyšte RR_REL nebo RR_ABS")
    return lut

def koduj_rr(rr, stred):
    k = np.searchsorted(stred, rr).astype(np.uint8)
    k[~np.isfinite(rr)] = KOD_CHYBI
    return k


# =============================== FILTRY (před gzipem) ===============================
def predikce_dlazdice(sub):
    """Pro každou buňku pevniny v dlaždici index souseda: vlevo, jinak nahoře, jinak -1."""
    rr_, cc_ = np.nonzero(sub)
    rank = np.full(sub.shape, -1, np.int32); rank[rr_, cc_] = np.arange(len(rr_), dtype=np.int32)
    left = np.full(len(rr_), -1, np.int32); m = cc_ > 0; left[m] = rank[rr_[m], cc_[m] - 1]
    up = np.full(len(rr_), -1, np.int32);   m = rr_ > 0; up[m] = rank[rr_[m] - 1, cc_[m]]
    return np.where(left >= 0, left, up)

def filtr(x, P, mode):
    """x: uint8 (buněk, kroků). Aritmetika modulo 256 (bezeztrátová)."""
    if mode in (1, 3):
        d = x.copy(); d[:, 1:] = x[:, 1:] - x[:, :-1]; x = d
    if mode in (2, 3):
        r = x.copy(); m = P >= 0; r[m] = x[m] - x[P[m]]; x = r
    return x

def gz(b, uroven=GZIP_UROVEN):
    return gzip.compress(b, compresslevel=uroven, mtime=0)

def vyber_filtr(nazev, ziskej, vzorek_id, predikce, pevny):
    if pevny is not None:
        print(f"  filtr {nazev}: pevně {pevny}"); return pevny
    vel = {m: 0 for m in range(4)}
    for i in vzorek_id:
        x = ziskej(i)
        for m in vel: vel[m] += len(gz(np.ascontiguousarray(filtr(x, predikce[i], m)).tobytes(), 6))
    nejlepsi = min(vel, key=vel.get)
    popis = ", ".join(f"{m}: {v / 1e6:.2f} MB" for m, v in vel.items())
    print(f"  filtr {nazev}: vzorek {len(vzorek_id)} dlaždic -> {popis}  => vybrán {nejlepsi}")
    return nejlepsi


# =============================== KALIBRACE ===============================
def sestav_kalibraci(cnt, S1, NB, NK):
    cnt2 = cnt.reshape(NB, NK).astype(np.float64)
    pool = cnt2.sum(0)
    kody = np.arange(NK)
    h1c = kody * H1_KROK
    tab = np.zeros((5, NB, NK))
    for k in range(5):
        s1 = S1[k].reshape(NB, NK)
        with np.errstate(invalid="ignore", divide="ignore"):
            m = np.where(cnt2 >= MIN_VZORKU, s1 / cnt2, np.nan)
            mp = np.where(pool >= MIN_VZORKU, s1.sum(0) / pool, np.nan)
        m = np.where(np.isnan(m), mp[None, :], m)
        for b in range(NB):
            ok = ~np.isnan(m[b])
            if not ok.any():
                m[b] = np.maximum(h1c - 2.0 * (k + 1), 0.0); continue
            i0, i1 = kody[ok][0], kody[ok][-1]
            m[b] = np.interp(kody, kody[ok], m[b][ok])
            m[b, kody > i1] = m[b, i1] + (kody[kody > i1] - i1) * H1_KROK   # teplá zóna: hladina roste 1:1 s H1
        tab[k] = m
    prev = h1c[None, :]
    for k in range(5):                                # monotonie: vyšší práh nikdy nemá víc než nižší
        tab[k] = np.clip(np.minimum(tab[k], prev), 0, None); prev = tab[k]
    return tab


# =============================== HLAVNÍ ===============================
def main():
    t_start = time.time()
    dvojice = najdi_soubory()
    print("=" * 70)
    print("SLUČOVAČ v2" + ("  – REŽIM ODHADU (nic se nezapisuje)" if ODHAD else ""))
    print("=" * 70)
    for a, b in dvojice: print(f"  {os.path.basename(a)}  +  {os.path.basename(b)}")

    # ---- geometrie z prvního souboru ----
    nc_rr = Dataset(dvojice[0][0]); nc_bd = Dataset(dvojice[0][1])
    v_rr = najdi_promennou(nc_rr, "RR"); v_bd = najdi_promennou(nc_bd, "BEDD")
    o_rr, o_bd = osy(v_rr), osy(v_bd)
    N_RR, N_BD = v_rr.name, v_bd.name
    if "l" not in o_bd: sys.exit("BEDD nemá osu hladin")
    lat = np.asarray(nc_rr.variables[v_rr.dimensions[o_rr["lat"]]][:], dtype=np.float64)
    lon = np.asarray(nc_rr.variables[v_rr.dimensions[o_rr["lon"]]][:], dtype=np.float64)
    flip = lat[0] < lat[-1]
    lat_nu = lat[::-1] if flip else lat                      # sever nahoře (jako Leaflet)
    NCOLS = len(lon); radku = int(np.sum(lat_nu >= LAT_MIN - 1e-6))
    dx = float(abs(lon[1] - lon[0])); dy = float(abs(lat[1] - lat[0]))
    lon0, lat0 = float(lon[0]), float(lat_nu[0])
    PY, PX = math.ceil(radku / TILE), math.ceil(NCOLS / TILE)
    print(f"  mřížka {NCOLS} x {len(lat)}, ukládá se {radku} řad (šířka >= {LAT_MIN}), dlaždice {TILE}x{TILE} -> {PX} x {PY}")

    # ---- délky, dny v dekádách, kontrola shody času ----
    nt_soubor, nd_all, limit = [], [], []
    for a, b in dvojice:
        with Dataset(a) as r, Dataset(b) as d:
            nt = len(r.dimensions["time"])
            if nt != len(d.dimensions["time"]): sys.exit(f"Chyba: různý počet kroků v {a} a {b}")
            if not np.allclose(r.variables["time"][:], d.variables["time"][:]):
                sys.exit(f"Chyba: časové osy {a} a {b} se liší")
            nd = dny_v_krocich(r, nt)
        nt_soubor.append(nt); nd_all.append(nd)
        limit.append(min(nt, ODHAD_KROKU) if ODHAD else nt)
    NT_CELKEM, NT = sum(nt_soubor), sum(limit)
    rok0 = rok_mesic_start(nc_rr)
    nc_rr.close(); nc_bd.close()
    nd_vse = np.concatenate([nd_all[i][:limit[i]] for i in range(len(dvojice))])
    print(f"  celkem {NT_CELKEM} dekád" + (f" (v odhadu {NT})" if ODHAD else "") + f", počátek roku {rok0}")

    # ---- maska pevniny z prvního kroku ----
    with Dataset(dvojice[0][0]) as r, Dataset(dvojice[0][1]) as d:
        rr0 = cti_blok(r.variables[N_RR], o_rr, 0, 1, flip, radku)[0]
        bd0 = cti_blok(d.variables[N_BD], o_bd, 0, 1, flip, radku)[0]
    maska = np.isfinite(rr0) & np.isfinite(bd0).all(axis=0)
    del rr0, bd0
    rows_l, cols_l, t_start_idx, t_n = [], [], [], []
    sub_masky = []
    pos = 0
    for ty in range(PY):
        for tx in range(PX):
            sub = np.zeros((TILE, TILE), bool)
            blok = maska[ty * TILE:(ty + 1) * TILE, tx * TILE:(tx + 1) * TILE]
            sub[:blok.shape[0], :blok.shape[1]] = blok
            r_, c_ = np.nonzero(sub)
            rows_l.append(ty * TILE + r_); cols_l.append(tx * TILE + c_)
            t_start_idx.append(pos); t_n.append(len(r_)); pos += len(r_)
            sub_masky.append(sub)
    rows_all = np.concatenate(rows_l); cols_all = np.concatenate(cols_l)
    NC = len(rows_all)
    pasmo_bunky = np.searchsorted(PASMA_LAT, np.abs(lat0 - dy * rows_all), side="right").astype(np.int64)
    NB, NK = len(PASMA_LAT) + 1, H1_MAX_KOD + 1
    print(f"  buněk pevniny: {NC}  (zbytek = moře/Antarktida se neukládá), neprázdných dlaždic: {sum(1 for n in t_n if n)}")
    print(f"  RAM pro kódy: {2 * NC * NT / 1e6:.0f} MB")

    lut = sestav_rr_lut(); stred = (lut[1:] + lut[:-1]) / 2
    rr_data = np.empty((NC, NT), np.uint8)
    h1_data = np.empty((NC, NT), np.uint8)

    # ---- pořadí hladin: z dat, ne z předpokladu ----
    poradi = None
    cnt = np.zeros(NB * NK, np.int64)
    S1 = [np.zeros(NB * NK) for _ in range(5)]
    S2 = [np.zeros(NB * NK) for _ in range(5)]
    max_rozdil, poruseni, orez, n_vz = 0.0, 0, 0, 0

    print("\n[1/3] Čtení .nc souborů a kódování")
    t0 = time.time(); g = 0
    for i, (fa, fb) in enumerate(dvojice):
        r = Dataset(fa); d = Dataset(fb)
        vr, vb = r.variables[N_RR], d.variables[N_BD]
        print(f"  soubor {i + 1}/{len(dvojice)}: {os.path.basename(fb)}")
        for a in range(0, limit[i], BLOK):
            b = min(a + BLOK, limit[i])
            rr = cti_blok(vr, o_rr, a, b, flip, radku)[:, rows_all, cols_all]
            bd = cti_blok(vb, o_bd, a, b, flip, radku)[:, :, rows_all, cols_all]    # (nb, hladiny, buněk)
            if poradi is None:
                prum = np.nanmean(bd, axis=(0, 2))
                poradi = np.argsort(-prum)                 # nejvyšší průměr = práh 0 °C ... nejnižší = 10 °C
                souradnice = None
                nm = d.variables[vb.dimensions[o_bd["l"]]] if vb.dimensions[o_bd["l"]] in d.variables else None
                if nm is not None:
                    hod = np.asarray(nm[:], dtype=float)
                    if len(hod) == 6: souradnice = np.argsort(hod)
                print("  pořadí hladin v souboru (index -> práh °C):",
                      ", ".join(f"{int(ix)}->{p}" for p, ix in zip(PRAHY, poradi)))
                if souradnice is not None and not np.array_equal(souradnice, poradi):
                    print("  !! POZOR: souřadnice 'threshold' říká jiné pořadí než data; používám pořadí podle dat")
            bd = bd[:, poradi, :]
            n = nd_vse[g:g + (b - a)].astype(np.float32)[:, None, None]
            h = bd / n                                       # °C/den
            h1 = h[:, 0, :]
            kod = np.clip(np.rint(h1 / H1_KROK), 0, H1_MAX_KOD)
            orez += int(np.nansum(h1 > H1_MAX_KOD * H1_KROK))
            kod[~np.isfinite(h1)] = KOD_CHYBI
            kod = kod.astype(np.uint8)
            rr_data[:, g:g + (b - a)] = koduj_rr(rr, stred).T
            h1_data[:, g:g + (b - a)] = kod.T
            # kalibrace: střední hodnota hladin 2..10 pro každý (pásmo, kód H1)
            ok = np.isfinite(h).all(axis=1) & (kod != KOD_CHYBI)           # (nb, buněk)
            if ok.any():
                ix = (pasmo_bunky[None, :] * NK + kod.astype(np.int64))[ok]
                cnt += np.bincount(ix, minlength=NB * NK)
                for k in range(5):
                    w = h[:, k + 1, :][ok].astype(np.float64)
                    S1[k] += np.bincount(ix, weights=w, minlength=NB * NK)
                    S2[k] += np.bincount(ix, weights=w * w, minlength=NB * NK)
                max_rozdil = max(max_rozdil, float(np.nanmax(h[:, 0, :] - h[:, 1, :])))
                poruseni += int(np.sum(h[:, 1:, :] > h[:, :-1, :] + 0.01))
                n_vz += int(h[:, 1:, :].size)
            g += b - a
            progres(g, NT, t0, "dekády")
        r.close(); d.close()
    del rr, bd, h, h1, kod

    print(f"\n  kontrola: max (BEDD0 - BEDD2)/počet dní = {max_rozdil:.2f}  (očekává se ≈ 2.00; "
          f"výrazně jiná hodnota = špatně určená délka dekády)")
    if poruseni > 0.001 * max(n_vz, 1):
        print(f"  !! {poruseni} případů, kdy vyšší práh má víc než nižší – zkontrolujte pořadí hladin")
    if orez: print(f"  upozornění: {orez} hodnot H1 nad {H1_MAX_KOD * H1_KROK:.0f} °C/den oříznuto")

    # ---- kalibrace ----
    print("\n[2/3] Kalibrační tabulka hladin 2–10 °C")
    tab = sestav_kalibraci(cnt, S1, NB, NK)
    Nvz = cnt.sum()
    nbar = float(nd_vse.mean())
    print("  typická chyba rekonstrukce (RMS, BEDD za dekádu):", end=" ")
    for k in range(5):
        with np.errstate(invalid="ignore", divide="ignore"):
            var = np.where(cnt > 0, S2[k] - S1[k] ** 2 / np.maximum(cnt, 1), 0).sum() / max(Nvz, 1)
        print(f"{PRAHY[k + 1]} °C: ±{math.sqrt(max(var, 0)) * nbar:.1f}", end="   ")
    print("\n  (chyba je náhodná se střední hodnotou ≈ 0, takže sezónní součty zůstávají nezkreslené)")

    # ---- výběr filtru a balení ----
    print("\n[3/3] Filtry a gzip dlaždic")
    predikce = [predikce_dlazdice(s) for s in sub_masky]
    nepr = [i for i in range(len(t_n)) if t_n[i] > 0]
    vz = [nepr[j] for j in np.linspace(0, len(nepr) - 1, min(6, len(nepr))).astype(int)]
    kus = lambda pole, i: pole[t_start_idx[i]:t_start_idx[i] + t_n[i]]
    f_rr = vyber_filtr("RR", lambda i: kus(rr_data, i), vz, predikce, FILTR_RR)
    f_h1 = vyber_filtr("H1", lambda i: kus(h1_data, i), vz, predikce, FILTR_H1)

    rejstrik, offset, celkem_gz = [], 0, 0
    fout = None if ODHAD else open(VYSTUP_DATA, "wb")
    t0 = time.time()
    for i in range(len(t_n)):
        if t_n[i] == 0:
            rejstrik.append({"offset": offset, "delka": 0})
        else:
            payload = (np.ascontiguousarray(filtr(kus(rr_data, i), predikce[i], f_rr)).tobytes() +
                       np.ascontiguousarray(filtr(kus(h1_data, i), predikce[i], f_h1)).tobytes())
            blok = gz(payload)
            if fout: fout.write(blok)
            rejstrik.append({"offset": offset, "delka": len(blok)})
            offset += len(blok); celkem_gz += len(blok)
        if i % 5 == 0 or i == len(t_n) - 1: progres(i + 1, len(t_n), t0, "dlaždice")
    if fout: fout.close()

    if ODHAD:
        odhad = celkem_gz * NT_CELKEM / NT / 1e6
        print("\n" + "=" * 70)
        print(f"ODHAD výsledné databáze pro všech {NT_CELKEM} dekád: ~{odhad:.0f} MB (±10–15 %)")
        print(f"  (vzorek {NT} dekád zabral {celkem_gz / 1e6:.1f} MB)")
        print("  Chcete menší? RR_REL = 0.08, H1_KROK = 0.25 nebo LAT_MIN = -56.")
        print("=" * 70); return

    maska_b64 = base64.b64encode(np.packbits(maska.ravel()).tobytes()).decode("ascii")
    reg = {
        "verze": 2,
        "mrizka": {"ncols": NCOLS, "radku": radku, "lon0": lon0, "lat0": lat0, "dx": dx, "dy": dy,
                   "tile": TILE, "dlazdic_x": PX, "dlazdic_y": PY},
        "pocet_kroku": NT, "rok_od": rok0,
        "dny": [int(x) for x in nd_vse],
        "maska_pevniny_b64": maska_b64,
        "rr_lut": [round(float(x), 3) for x in lut],
        "kod_chybi": KOD_CHYBI,
        "h1_krok": H1_KROK, "h1_max_kod": H1_MAX_KOD,
        "pasma_lat": PASMA_LAT,
        "praha_hladin": PRAHY[1:],
        "kalibrace": [[[round(float(x), 3) for x in tab[k, b]] for b in range(NB)] for k in range(5)],
        "filtr_rr": int(f_rr), "filtr_h1": int(f_h1),
        "dlazdice": rejstrik,
    }
    with open(VYSTUP_REJSTRIK, "w", encoding="utf-8") as f: json.dump(reg, f, ensure_ascii=False, separators=(",", ":"))
    print("\n" + "=" * 70)
    print("HOTOVO")
    print(f"  data:     {VYSTUP_DATA}  {os.path.getsize(VYSTUP_DATA) / 1e6:.1f} MB")
    print(f"  rejstřík: {VYSTUP_REJSTRIK}  {os.path.getsize(VYSTUP_REJSTRIK) / 1e3:.0f} kB")
    print(f"  celkový čas {cas(time.time() - t_start)}")
    print("=" * 70)


if __name__ == "__main__":
    main()
